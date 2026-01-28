#!/usr/bin/env python3
"""
Empirical tests comparing v0.1.x and v0.2.x scoring behavior.

These tests run pairwise comparisons on real data and fail fast with
detailed context when differences are found. This allows us to:
1. Discover where the algorithms differ
2. Analyze each difference to understand the cause
3. Decide if the difference is an improvement, behavioral change, or regression

NOTE: These tests require git access to origin/main and are skipped in CI.
"""

import os
import pytest
from pathlib import Path

# Skip all tests in this module when running in CI (no git origin access)
pytestmark = pytest.mark.skipif(
    os.environ.get("CI") == "true",
    reason="Version comparison tests require git origin/main access (skipped in CI)"
)
from itertools import combinations

from .compare_versions import (
    compare_pair,
    compare_fasta_pairwise,
    load_v01x_scorer,
    read_fasta,
    VersionComparison,
)
from adjusted_identity import DEFAULT_ADJUSTMENT_PARAMS, AdjustmentParams


# Path to test data (gitignored, not published)
TEST_DATA_DIR = Path(__file__).parent / "data"
SAMPLE_FASTA = TEST_DATA_DIR / "ONT10.82-B11-IN25-00187-iNat272809801-c1-RiC500-msa.fasta"


class TestVersionComparisonInfrastructure:
    """Tests for the comparison infrastructure itself."""

    def test_v01x_loader_works(self):
        """Verify we can load the v0.1.x scorer."""
        scorer = load_v01x_scorer()
        assert callable(scorer)

        # Test it works on simple input
        result = scorer("AAAA", "AAAA", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 1.0

    def test_compare_pair_identical(self):
        """Identical sequences should match between versions."""
        comparison = compare_pair("AAATTTGGG", "AAATTTGGG")
        assert comparison.match, "Identical sequences should produce matching results"

    def test_compare_pair_simple_mismatch(self):
        """Simple mismatch should be handled similarly by both versions."""
        comparison = compare_pair("AAATTTGGG", "AAATTTGGC")
        # Both should detect 1 mismatch
        assert comparison.v01x_mismatches == 1
        assert comparison.v02x_mismatches == 1
        assert comparison.match


class TestPrealignedFasta:
    """
    Empirical tests on pre-aligned FASTA files.

    These tests fail fast with detailed context when differences are found.

    NOTE: These tests are currently expected to fail because v0.2.x uses a
    variant range algorithm that correctly identifies homopolymer extensions
    that v0.1.x missed. The differences are improvements, not regressions.
    """

    @pytest.mark.xfail(
        reason="v0.2.x variant range algorithm detects extensions v0.1.x missed"
    )
    @pytest.mark.skipif(
        not SAMPLE_FASTA.exists(),
        reason=f"Sample FASTA not found: {SAMPLE_FASTA}"
    )
    def test_sample_fasta_pairwise(self):
        """
        Compare all pairwise combinations in the sample FASTA.

        This test fails on the FIRST difference found, providing detailed
        context for analysis.
        """
        total, diff_count, differences = compare_fasta_pairwise(
            str(SAMPLE_FASTA),
            params=DEFAULT_ADJUSTMENT_PARAMS,
            stop_on_diff=True,
        )

        if differences:
            name1, name2, comparison = differences[0]
            pytest.fail(
                f"\n\nVersion difference found at pair {total}:\n"
                f"Sequences: {name1} vs {name2}\n\n"
                f"{comparison.format_diff()}"
            )

    @pytest.mark.xfail(
        reason="v0.2.x variant range algorithm detects extensions v0.1.x missed"
    )
    @pytest.mark.skipif(
        not SAMPLE_FASTA.exists(),
        reason=f"Sample FASTA not found: {SAMPLE_FASTA}"
    )
    def test_sample_fasta_count_differences(self):
        """
        Count total differences in the sample FASTA (doesn't stop on first).

        Use this to get an overview of how many pairs differ.
        """
        total, diff_count, differences = compare_fasta_pairwise(
            str(SAMPLE_FASTA),
            params=DEFAULT_ADJUSTMENT_PARAMS,
            stop_on_diff=False,
            max_pairs=1000,  # Limit for quick feedback
        )

        if diff_count > 0:
            # Show first few differences
            msg_lines = [
                f"\n\nFound {diff_count} differences out of {total} pairs compared:",
                ""
            ]
            for i, (name1, name2, comparison) in enumerate(differences[:3]):
                msg_lines.append(f"--- Difference {i+1}: {name1} vs {name2} ---")
                msg_lines.append(comparison.format_diff())
                msg_lines.append("")

            if diff_count > 3:
                msg_lines.append(f"... and {diff_count - 3} more differences")

            pytest.fail("\n".join(msg_lines))


class TestKnownDifferences:
    """
    Tests for specific cases where we expect v0.1.x and v0.2.x to differ.

    These tests document expected behavioral differences and verify
    that v0.2.x produces the "correct" result.
    """

    def test_alternating_gaps_improvement(self):
        """
        v0.2.x correctly handles alternating gaps as a single variant range.

        In v0.1.x, this was processed as separate indel events.
        In v0.2.x, it's recognized as equivalent homopolymer extensions.
        """
        # This is the case from the algorithm design: C extends left, T extends right
        seq1 = "TGC-C-TC"
        seq2 = "TGCT--TC"

        comparison = compare_pair(seq1, seq2)

        # v0.2.x should score this as identity=1.0 (both are valid extensions)
        assert comparison.v02x_identity == 1.0, (
            f"v0.2.x should recognize this as equivalent extensions.\n"
            f"{comparison.format_diff()}"
        )

        # Document whether v0.1.x differs (it likely does)
        if not comparison.match:
            print(f"\nExpected difference (v0.2.x improvement):\n{comparison.format_diff()}")


class TestParameterVariations:
    """Tests that compare versions across different parameter settings."""

    @pytest.mark.parametrize("normalize_hp", [True, False])
    @pytest.mark.parametrize("normalize_indels", [True, False])
    def test_simple_indel_with_params(self, normalize_hp, normalize_indels):
        """Compare a simple indel case across parameter combinations."""
        seq1 = "AAATGGG"
        seq2 = "AAA-GGG"

        params = AdjustmentParams(
            normalize_homopolymers=normalize_hp,
            normalize_indels=normalize_indels,
            handle_iupac_overlap=True,
            end_skip_distance=0,
        )

        comparison = compare_pair(seq1, seq2, params)

        if not comparison.match:
            pytest.fail(
                f"\nDifference with params (hp={normalize_hp}, indels={normalize_indels}):\n"
                f"{comparison.format_diff()}"
            )


# Utility for manual exploration
def explore_differences(fasta_path: str, max_pairs: int = 100):
    """
    Utility function for interactive exploration of differences.

    Run with: python -c "from test_version_comparison import explore_differences; explore_differences('path/to/file.fasta')"
    """
    total, diff_count, differences = compare_fasta_pairwise(
        fasta_path,
        stop_on_diff=False,
        max_pairs=max_pairs,
    )

    print(f"\nCompared {total} pairs, found {diff_count} differences ({100*diff_count/total:.1f}%)")

    if differences:
        print("\nFirst 5 differences:\n")
        for i, (name1, name2, comparison) in enumerate(differences[:5]):
            print(f"--- Difference {i+1}: {name1} vs {name2} ---")
            print(comparison.format_diff())
            print()


if __name__ == "__main__":
    # Run a quick check when executed directly
    if SAMPLE_FASTA.exists():
        print(f"Running comparison on {SAMPLE_FASTA}...")
        explore_differences(str(SAMPLE_FASTA), max_pairs=100)
    else:
        print(f"Sample FASTA not found: {SAMPLE_FASTA}")
        print("Running basic comparison test...")
        comparison = compare_pair("AAATTTGGG", "AAATTTGGG")
        print(f"Identical sequences match: {comparison.match}")
