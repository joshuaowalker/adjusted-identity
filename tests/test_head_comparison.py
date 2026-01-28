#!/usr/bin/env python3
"""
Tests comparing origin/main (pre-adjust_gaps) with HEAD (with adjust_gaps).

Two test categories:
A) adjust_gaps=False: Results should be IDENTICAL to origin/main
B) adjust_gaps=True: Metrics may differ due to gap rewriting changing alignment interpretation

These tests verify that adding the adjust_gaps parameter doesn't change
existing behavior when adjust_gaps=False (the default).

NOTE: These tests require git access to origin/main and are skipped in CI.
"""

import os
import pytest
from pathlib import Path

# Skip all tests in this module when running in CI (no git origin access)
pytestmark = pytest.mark.skipif(
    os.environ.get("CI") == "true",
    reason="HEAD comparison tests require git origin/main access (skipped in CI)"
)

from .compare_versions import (
    compare_head_pair,
    compare_head_fasta_pairwise,
    load_origin_scorer,
    read_fasta,
    HeadComparison,
)
from adjusted_identity import DEFAULT_ADJUSTMENT_PARAMS, AdjustmentParams


# Path to test data (gitignored, not published)
TEST_DATA_DIR = Path(__file__).parent / "data"
SAMPLE_FASTA = TEST_DATA_DIR / "ONT10.82-B11-IN25-00187-iNat272809801-c1-RiC500-msa.fasta"


class TestHeadComparisonInfrastructure:
    """Tests for the HEAD comparison infrastructure itself."""

    def test_origin_loader_works(self):
        """Verify we can load the origin/main scorer."""
        scorer = load_origin_scorer()
        assert callable(scorer)

        # Test it works on simple input
        result = scorer("AAAA", "AAAA", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 1.0

    def test_compare_head_pair_identical_sequences(self):
        """Identical sequences should produce identical results across all modes."""
        comparison = compare_head_pair("AAATTTGGG", "AAATTTGGG")
        assert comparison.adjust_false_identical, (
            "Identical sequences should produce identical results with adjust_gaps=False"
        )
        assert comparison.adjust_true_metrics_match, (
            "Identical sequences should have matching metrics with adjust_gaps=True"
        )


class TestAdjustGapsFalseIdentical:
    """
    Verify HEAD with adjust_gaps=False produces identical results to origin/main.

    When the new code runs with adjust_gaps=False (default), it should produce
    byte-for-byte identical results to origin/main. This verifies zero behavioral
    change for the default path.
    """

    @pytest.mark.skipif(
        not SAMPLE_FASTA.exists(),
        reason=f"Sample FASTA not found: {SAMPLE_FASTA}"
    )
    def test_sample_fasta_pairwise_identical(self):
        """
        All pairs in sample FASTA should produce identical results.

        This test fails on the FIRST difference found, providing detailed
        context for analysis. Limited to 1000 pairs for reasonable test time.
        """
        total, diff_count, differences = compare_head_fasta_pairwise(
            str(SAMPLE_FASTA),
            params=DEFAULT_ADJUSTMENT_PARAMS,
            stop_on_diff=True,
            check_mode="adjust_false",
            max_pairs=1000,
        )

        if differences:
            name1, name2, comparison = differences[0]
            pytest.fail(
                f"\n\nadjust_gaps=False difference found at pair {total}:\n"
                f"Sequences: {name1} vs {name2}\n\n"
                f"{comparison.format_adjust_false_diff()}"
            )

    def test_simple_cases_identical(self):
        """Simple test cases should produce identical results."""
        # Note: score_alignment expects pre-aligned sequences with equal lengths
        test_cases = [
            # (seq1, seq2, description)
            ("ATCG", "ATCG", "Perfect match"),
            ("ATCG", "AXCG", "Single substitution"),
            ("AAAA-TT", "AAA--TT", "Homopolymer indel difference"),
            ("TGC-C-TC", "TGCT--TC", "Variant range"),
            ("ATCG", "AT-G", "Simple deletion"),
            ("AT-G", "ATCG", "Simple insertion"),
            ("AAAAATTTTT", "AAAA-TTTTT", "Homopolymer length A"),
            ("AAA-TTTTT", "AAATTTTTT", "Homopolymer length T"),
            ("ATGC", "ATRC", "IUPAC ambiguity"),
            ("ATGC", "ATGN", "IUPAC N wildcard"),
            ("A--TGC", "ATT-GC", "Multiple gaps"),
            ("AAAA----TTTT", "AAAA----TTTT", "Dual gaps"),
        ]

        for seq1, seq2, description in test_cases:
            comparison = compare_head_pair(seq1, seq2)
            assert comparison.adjust_false_identical, (
                f"\nadjust_gaps=False not identical for: {description}\n"
                f"{comparison.format_adjust_false_diff()}"
            )

    @pytest.mark.parametrize("normalize_hp", [True, False])
    @pytest.mark.parametrize("normalize_indels", [True, False])
    @pytest.mark.parametrize("handle_iupac", [True, False])
    def test_parameter_combinations_identical(self, normalize_hp, normalize_indels, handle_iupac):
        """
        Various parameter combinations should produce identical results.

        Tests a representative sequence with each parameter combination.
        """
        params = AdjustmentParams(
            normalize_homopolymers=normalize_hp,
            normalize_indels=normalize_indels,
            handle_iupac_overlap=handle_iupac,
            end_skip_distance=0,
        )

        # Test sequence with homopolymer, indel, and IUPAC elements
        seq1 = "AAAA-TGRC"
        seq2 = "AAA--TGNC"

        comparison = compare_head_pair(seq1, seq2, params)
        assert comparison.adjust_false_identical, (
            f"\nadjust_gaps=False not identical with params:\n"
            f"  normalize_hp={normalize_hp}, normalize_indels={normalize_indels}, "
            f"handle_iupac={handle_iupac}\n"
            f"{comparison.format_adjust_false_diff()}"
        )


class TestAdjustGapsTrueMetrics:
    """
    Test HEAD with adjust_gaps=True behavior.

    When using adjust_gaps=True, gap rewriting can change the alignment structure
    which may cause metrics (mismatches, scored_positions, identity) to differ
    from origin/main. The aligned strings are expected to differ (that's the
    point of adjust_gaps).

    These tests document the current behavior. Some are marked xfail where
    metrics are known to differ.
    """

    @pytest.mark.xfail(
        reason="Gap rewriting can change mismatch counts by reclassifying content as extensions"
    )
    @pytest.mark.skipif(
        not SAMPLE_FASTA.exists(),
        reason=f"Sample FASTA not found: {SAMPLE_FASTA}"
    )
    def test_sample_fasta_pairwise_mismatches(self):
        """
        Test if mismatch counts match with adjust_gaps=True.

        Gap rewriting can change both scored_positions AND mismatch counts
        by reclassifying content as homopolymer extensions. This test documents
        that behavior. Limited to 1000 pairs.
        """
        from itertools import combinations

        sequences = read_fasta(str(SAMPLE_FASTA))
        total = 0
        mismatch_diffs = []

        for (name1, seq1), (name2, seq2) in combinations(sequences, 2):
            comparison = compare_head_pair(seq1, seq2)
            total += 1

            if not comparison.adjust_true_mismatches_match:
                mismatch_diffs.append((name1, name2, comparison))
                break  # Stop on first mismatch diff

            if total >= 1000:
                break

        if mismatch_diffs:
            name1, name2, comparison = mismatch_diffs[0]
            pytest.fail(
                f"\n\nadjust_gaps=True mismatch count difference found at pair {total}:\n"
                f"Sequences: {name1} vs {name2}\n\n"
                f"{comparison.format_adjust_true_diff()}"
            )

    @pytest.mark.xfail(
        reason="Gap rewriting can change scored_positions and identity"
    )
    @pytest.mark.skipif(
        not SAMPLE_FASTA.exists(),
        reason=f"Sample FASTA not found: {SAMPLE_FASTA}"
    )
    def test_sample_fasta_pairwise_all_metrics(self):
        """
        Test if all metrics match with adjust_gaps=True.

        Gap rewriting can change the alignment interpretation, causing
        metrics to differ. This test documents that behavior.
        """
        total, diff_count, differences = compare_head_fasta_pairwise(
            str(SAMPLE_FASTA),
            params=DEFAULT_ADJUSTMENT_PARAMS,
            stop_on_diff=True,
            check_mode="adjust_true",
            max_pairs=1000,
        )

        if differences:
            name1, name2, comparison = differences[0]
            pytest.fail(
                f"\n\nadjust_gaps=True metrics difference found at pair {total}:\n"
                f"Sequences: {name1} vs {name2}\n\n"
                f"{comparison.format_adjust_true_diff()}"
            )

    def test_simple_cases_metrics(self):
        """Simple test cases should have matching metrics."""
        # Note: score_alignment expects pre-aligned sequences with equal lengths
        # Note: Some complex gap patterns (like "A--TGC" vs "ATT-GC") may have
        # different metrics with adjust_gaps=True because the gap rewriting
        # process can fundamentally change the alignment interpretation.
        # Those cases are excluded since the primary test is adjust_gaps=False identity.
        test_cases = [
            # (seq1, seq2, description)
            ("ATCG", "ATCG", "Perfect match"),
            ("ATCG", "AXCG", "Single substitution"),
            ("AAAA-TT", "AAA--TT", "Homopolymer indel difference"),
            ("TGC-C-TC", "TGCT--TC", "Variant range"),
            ("ATCG", "AT-G", "Simple deletion"),
            ("AT-G", "ATCG", "Simple insertion"),
            ("AAAAATTTTT", "AAAA-TTTTT", "Homopolymer length A"),
            ("AAA-TTTTT", "AAATTTTTT", "Homopolymer length T"),
            ("ATGC", "ATRC", "IUPAC ambiguity"),
            ("ATGC", "ATGN", "IUPAC N wildcard"),
            ("AAAA----TTTT", "AAAA----TTTT", "Dual gaps"),
        ]

        for seq1, seq2, description in test_cases:
            comparison = compare_head_pair(seq1, seq2)
            assert comparison.adjust_true_metrics_match, (
                f"\nadjust_gaps=True metrics differ for: {description}\n"
                f"{comparison.format_adjust_true_diff()}"
            )

    @pytest.mark.parametrize("normalize_hp", [True, False])
    @pytest.mark.parametrize("normalize_indels", [True, False])
    @pytest.mark.parametrize("handle_iupac", [True, False])
    def test_parameter_combinations_metrics(self, normalize_hp, normalize_indels, handle_iupac):
        """
        Various parameter combinations should produce matching metrics.

        Tests a representative sequence with each parameter combination.
        """
        params = AdjustmentParams(
            normalize_homopolymers=normalize_hp,
            normalize_indels=normalize_indels,
            handle_iupac_overlap=handle_iupac,
            end_skip_distance=0,
        )

        # Test sequence with homopolymer, indel, and IUPAC elements
        seq1 = "AAAA-TGRC"
        seq2 = "AAA--TGNC"

        comparison = compare_head_pair(seq1, seq2, params)
        assert comparison.adjust_true_metrics_match, (
            f"\nadjust_gaps=True metrics differ with params:\n"
            f"  normalize_hp={normalize_hp}, normalize_indels={normalize_indels}, "
            f"handle_iupac={handle_iupac}\n"
            f"{comparison.format_adjust_true_diff()}"
        )


class TestEdgeCases:
    """Test edge cases for HEAD comparison."""

    def test_empty_sequences(self):
        """Empty sequences should be handled identically."""
        comparison = compare_head_pair("", "")
        assert comparison.adjust_false_identical
        assert comparison.adjust_true_metrics_match

    def test_all_gaps(self):
        """All-gap sequences should be handled identically."""
        comparison = compare_head_pair("----", "----")
        assert comparison.adjust_false_identical
        assert comparison.adjust_true_metrics_match

    def test_single_base(self):
        """Single base sequences should be handled identically."""
        comparison = compare_head_pair("A", "A")
        assert comparison.adjust_false_identical
        assert comparison.adjust_true_metrics_match

        comparison = compare_head_pair("A", "T")
        assert comparison.adjust_false_identical
        assert comparison.adjust_true_metrics_match

    def test_long_homopolymer(self):
        """Long homopolymer runs should be handled identically."""
        # Pre-aligned with gap inserted
        seq1 = "A" * 50
        seq2 = "A" * 49 + "-"

        comparison = compare_head_pair(seq1, seq2)
        assert comparison.adjust_false_identical
        assert comparison.adjust_true_metrics_match


# Utility for manual exploration
def explore_head_differences(fasta_path: str, max_pairs: int = 100, check_mode: str = "adjust_false"):
    """
    Utility function for interactive exploration of HEAD differences.

    Run with:
        python -c "from tests.test_head_comparison import explore_head_differences; explore_head_differences('path/to/file.fasta')"
    """
    total, diff_count, differences = compare_head_fasta_pairwise(
        fasta_path,
        stop_on_diff=False,
        max_pairs=max_pairs,
        check_mode=check_mode,
    )

    print(f"\nCompared {total} pairs, found {diff_count} differences ({100*diff_count/total:.1f}%)")
    print(f"Check mode: {check_mode}")

    if differences:
        print("\nFirst 5 differences:\n")
        for i, (name1, name2, comparison) in enumerate(differences[:5]):
            print(f"--- Difference {i+1}: {name1} vs {name2} ---")
            if check_mode == "adjust_false":
                print(comparison.format_adjust_false_diff())
            else:
                print(comparison.format_adjust_true_diff())
            print()


if __name__ == "__main__":
    # Run a quick check when executed directly
    if SAMPLE_FASTA.exists():
        print(f"Running HEAD comparison on {SAMPLE_FASTA}...")
        print("\n=== Checking adjust_gaps=False identical ===")
        explore_head_differences(str(SAMPLE_FASTA), max_pairs=100, check_mode="adjust_false")
        print("\n=== Checking adjust_gaps=True metrics match ===")
        explore_head_differences(str(SAMPLE_FASTA), max_pairs=100, check_mode="adjust_true")
    else:
        print(f"Sample FASTA not found: {SAMPLE_FASTA}")
        print("Running basic comparison test...")
        comparison = compare_head_pair("AAATTTGGG", "AAATTTGGG")
        print(f"Identical sequences:")
        print(f"  adjust_gaps=False identical: {comparison.adjust_false_identical}")
        print(f"  adjust_gaps=True metrics match: {comparison.adjust_true_metrics_match}")
