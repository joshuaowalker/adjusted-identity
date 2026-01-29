"""
Test file documenting metric values with unified architecture.

This file documents the behavior of the unified single-pass architecture where
both adjust_gaps=True and adjust_gaps=False produce identical metrics.

Created as part of: Plan to Unify adjust_gaps Code Paths into Single-Pass Architecture
"""

import pytest
from adjusted_identity import (
    score_alignment,
    align_and_score,
    DEFAULT_ADJUSTMENT_PARAMS,
    RAW_ADJUSTMENT_PARAMS,
    AdjustmentParams,
)


class TestPureExtensionCases:
    """
    Document behavior when both alleles are pure extensions.

    These are cases where the variant range content can be fully explained
    as repeat extensions of the surrounding context.
    """

    def test_simple_homopolymer_extension(self):
        """Simple case: one extra A in homopolymer run."""
        # AAA vs AAAA - the extra A is a pure extension
        seq1 = "AAAA-TT"
        seq2 = "AAA--TT"

        result_false = score_alignment(seq1, seq2, adjust_gaps=False)
        result_true = score_alignment(seq1, seq2, adjust_gaps=True)

        # Document unified behavior:
        # Both adjust_gaps modes use the same _analyze_alignment pass
        # Only the output strings differ

        # Both should have identity=1.0 (homopolymer difference ignored)
        assert result_false.identity == 1.0
        assert result_true.identity == 1.0

        # Mismatches should be 0 in both cases (pure extension = 0 edits)
        assert result_false.mismatches == 0
        assert result_true.mismatches == 0

        # With unified architecture, scored_positions are always identical
        assert result_false.scored_positions == result_true.scored_positions
        print(f"scored_positions={result_false.scored_positions} (same for both modes)")

    def test_both_alleles_pure_extension(self):
        """Both alleles are pure extensions of context."""
        # Both sequences have content that extends the surrounding context
        # Left context: T, Right context: G
        # seq1 has "CC" which doesn't extend either context
        # seq2 has "T" which extends left context
        # This is NOT a both-pure-extension case
        seq1 = "TGC-C-TC"
        seq2 = "TGCT--TC"

        result_false = score_alignment(seq1, seq2, adjust_gaps=False)
        result_true = score_alignment(seq1, seq2, adjust_gaps=True)

        # Document current behavior
        print(f"TGC-C-TC vs TGCT--TC")
        print(f"adjust_gaps=False: identity={result_false.identity}, mismatches={result_false.mismatches}, scored={result_false.scored_positions}")
        print(f"adjust_gaps=True:  identity={result_true.identity}, mismatches={result_true.mismatches}, scored={result_true.scored_positions}")

        # Identity should be the same between both modes
        assert result_false.identity == result_true.identity

    def test_dinucleotide_pure_extension(self):
        """Dinucleotide repeat pure extension."""
        # ATAT context, one has extra AT
        seq1 = "GATAT--G"
        seq2 = "GATATATG"

        result_false = score_alignment(seq1, seq2, adjust_gaps=False)
        result_true = score_alignment(seq1, seq2, adjust_gaps=True)

        print(f"GATAT--G vs GATATATG")
        print(f"adjust_gaps=False: identity={result_false.identity}, mismatches={result_false.mismatches}, scored={result_false.scored_positions}")
        print(f"adjust_gaps=True:  identity={result_true.identity}, mismatches={result_true.mismatches}, scored={result_true.scored_positions}")

        # Both are pure extensions, should have identity=1.0
        assert result_false.identity == 1.0
        assert result_true.identity == 1.0


class TestMixedExtensionCoreCases:
    """
    Document behavior when alleles have mixed extension and core content.
    """

    def test_extension_plus_core(self):
        """One allele has extension + core, other is pure extension."""
        # seq1: G (extends right G context) + X (core mismatch)
        # seq2: empty (gap) - extends both contexts trivially
        seq1 = "AGX-G"
        seq2 = "AG--G"

        result_false = score_alignment(seq1, seq2, adjust_gaps=False)
        result_true = score_alignment(seq1, seq2, adjust_gaps=True)

        print(f"AGX-G vs AG--G")
        print(f"adjust_gaps=False: identity={result_false.identity}, mismatches={result_false.mismatches}, scored={result_false.scored_positions}")
        print(f"adjust_gaps=True:  identity={result_true.identity}, mismatches={result_true.mismatches}, scored={result_true.scored_positions}")

    def test_both_have_core_different_content(self):
        """Both alleles have core content that differs."""
        # Both have content that doesn't match context, and they differ
        seq1 = "AAA-X-GGG"
        seq2 = "AAAYY-GGG"

        result_false = score_alignment(seq1, seq2, adjust_gaps=False)
        result_true = score_alignment(seq1, seq2, adjust_gaps=True)

        print(f"AAA-X-GGG vs AAAYY-GGG")
        print(f"adjust_gaps=False: identity={result_false.identity}, mismatches={result_false.mismatches}, scored={result_false.scored_positions}")
        print(f"adjust_gaps=True:  identity={result_true.identity}, mismatches={result_true.mismatches}, scored={result_true.scored_positions}")

        # Identity should match
        assert result_false.identity == result_true.identity


class TestUnifiedMetricsCases:
    """
    Document cases that previously had different metrics for adjust_gaps=True
    vs adjust_gaps=False, but now produce identical metrics with unified architecture.
    """

    def test_opposite_direction_extensions(self):
        """
        Classic divergent case: opposite direction extensions.

        seq1: TGC-C-TC  -> variant range: "CC"
        seq2: TGCT--TC  -> variant range: "T"

        Context: left=C (from TGC), right=T (from TC)
        seq1 "CC": one C extends left, other C is core
        seq2 "T": extends right T context

        This case can produce different scored_positions because gap rewriting
        changes how the positions are interpreted.
        """
        seq1 = "TGC-C-TC"
        seq2 = "TGCT--TC"

        result_false = score_alignment(seq1, seq2, adjust_gaps=False)
        result_true = score_alignment(seq1, seq2, adjust_gaps=True)

        print(f"\nOpposite direction extensions:")
        print(f"seq1: {seq1}")
        print(f"seq2: {seq2}")
        print(f"adjust_gaps=False:")
        print(f"  identity={result_false.identity}")
        print(f"  mismatches={result_false.mismatches}")
        print(f"  scored_positions={result_false.scored_positions}")
        print(f"  score_aligned: {result_false.score_aligned}")
        print(f"adjust_gaps=True:")
        print(f"  identity={result_true.identity}")
        print(f"  mismatches={result_true.mismatches}")
        print(f"  scored_positions={result_true.scored_positions}")
        print(f"  score_aligned: {result_true.score_aligned}")

        # Document: these are equal due to unified architecture
        assert result_false.identity == result_true.identity
        assert result_false.mismatches == result_true.mismatches

    def test_complex_variant_region(self):
        """
        Complex variant region with multiple differences.

        These sequences have:
        - Homopolymer variation
        - Substitution
        """
        seq1 = "CACAGGCTGGTAAT"
        seq2 = "CACAAATTGGTAAT"

        result_false = align_and_score(seq1, seq2, adjust_gaps=False)
        result_true = align_and_score(seq1, seq2, adjust_gaps=True)

        print(f"\nComplex variant region:")
        print(f"seq1: {seq1}")
        print(f"seq2: {seq2}")
        print(f"adjust_gaps=False:")
        print(f"  aligned: {result_false.seq1_aligned}")
        print(f"           {result_false.seq2_aligned}")
        print(f"  identity={result_false.identity}")
        print(f"  mismatches={result_false.mismatches}")
        print(f"  scored_positions={result_false.scored_positions}")
        print(f"adjust_gaps=True:")
        print(f"  aligned: {result_true.seq1_aligned}")
        print(f"           {result_true.seq2_aligned}")
        print(f"  identity={result_true.identity}")
        print(f"  mismatches={result_true.mismatches}")
        print(f"  scored_positions={result_true.scored_positions}")

        # Document current values
        assert result_false.identity == result_true.identity
        assert result_false.mismatches == result_true.mismatches


class TestMetricValuesSnapshot:
    """
    Snapshot of exact metric values for key test cases with unified architecture.

    Both adjust_gaps=True and adjust_gaps=False produce identical metrics
    using the single _analyze_alignment pass.
    """

    def test_snapshot_homopolymer_extension(self):
        """Snapshot: AAAA-TT vs AAA--TT"""
        seq1 = "AAAA-TT"
        seq2 = "AAA--TT"

        result = score_alignment(seq1, seq2, adjust_gaps=False)

        # Document exact current values
        assert result.identity == 1.0
        assert result.mismatches == 0
        assert result.scored_positions == 5  # 3 A's + 2 T's (extension not counted)

    def test_snapshot_simple_mismatch(self):
        """Snapshot: ATCG vs AXCG"""
        seq1 = "ATCG"
        seq2 = "AXCG"

        result = score_alignment(seq1, seq2, adjust_gaps=False)

        assert result.identity == 0.75  # 3/4
        assert result.mismatches == 1
        assert result.scored_positions == 4

    def test_snapshot_dual_gaps(self):
        """Snapshot: sequences with dual gaps."""
        seq1 = "AA--TT"
        seq2 = "AA--TT"

        result = score_alignment(seq1, seq2, adjust_gaps=False)

        # Dual gaps are not scored
        assert result.identity == 1.0
        assert result.mismatches == 0
        assert result.scored_positions == 4  # Only the non-gap positions

    def test_snapshot_opposite_extensions(self):
        """Snapshot: TGC-C-TC vs TGCT--TC (opposite direction extensions)."""
        seq1 = "TGC-C-TC"
        seq2 = "TGCT--TC"

        result = score_alignment(seq1, seq2, adjust_gaps=False)

        # Document current behavior
        # This case has:
        # - Left context: C
        # - Right context: T
        # - seq1 variant: "CC" (one extends left, one is core)
        # - seq2 variant: "T" (extends right)
        # With HP normalization: 1 mismatch (the non-extension C in seq1 vs nothing in seq2)
        print(f"\nTGC-C-TC vs TGCT--TC snapshot:")
        print(f"  identity={result.identity}")
        print(f"  mismatches={result.mismatches}")
        print(f"  scored_positions={result.scored_positions}")
        print(f"  score_aligned: {result.score_aligned}")


class TestAlignAndScoreMetrics:
    """
    Test metrics from align_and_score (which performs alignment first).
    """

    def test_user_sequences_1(self):
        """Test from user-provided sequences issue."""
        seq1 = "TTTTCACAGGCTGGTAATGGCT"
        seq2 = "TTTTCACAAGTTGGTAATGGCT"

        result_false = align_and_score(seq1, seq2, adjust_gaps=False)
        result_true = align_and_score(seq1, seq2, adjust_gaps=True)

        print(f"\nUser sequences 1:")
        print(f"seq1: {seq1}")
        print(f"seq2: {seq2}")
        print(f"adjust_gaps=False: identity={result_false.identity}, mismatches={result_false.mismatches}, scored={result_false.scored_positions}")
        print(f"adjust_gaps=True:  identity={result_true.identity}, mismatches={result_true.mismatches}, scored={result_true.scored_positions}")

        # Both should have 1 mismatch
        assert result_false.mismatches == 1
        assert result_true.mismatches == 1

        # Identity should be ~95%
        assert result_false.identity == pytest.approx(20/21, rel=1e-4)
        assert result_true.identity == pytest.approx(20/21, rel=1e-4)

    def test_user_sequences_2(self):
        """Test from user-provided sequences issue 2."""
        seq1 = "CACAGGCTGGTAAT"
        seq2 = "CACAAATTGGTAAT"

        result_false = align_and_score(seq1, seq2, adjust_gaps=False)
        result_true = align_and_score(seq1, seq2, adjust_gaps=True)

        print(f"\nUser sequences 2:")
        print(f"seq1: {seq1}")
        print(f"seq2: {seq2}")
        print(f"adjust_gaps=False: identity={result_false.identity}, mismatches={result_false.mismatches}, scored={result_false.scored_positions}")
        print(f"adjust_gaps=True:  identity={result_true.identity}, mismatches={result_true.mismatches}, scored={result_true.scored_positions}")

        # Both should have 1 mismatch
        assert result_false.mismatches == 1
        assert result_true.mismatches == 1

        # Identity should be ~91.67%
        assert result_false.identity == pytest.approx(11/12, rel=1e-4)
        assert result_true.identity == pytest.approx(11/12, rel=1e-4)


if __name__ == "__main__":
    # Run with verbose output to see printed diagnostics
    pytest.main([__file__, "-v", "-s"])
