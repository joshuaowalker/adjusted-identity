"""
Tests for the unified single-pass architecture.

These tests verify that:
1. Both adjust_gaps=True and adjust_gaps=False produce identical metrics
2. Only the output alignment strings differ between modes
3. The single analysis produces correct metrics
"""

import pytest
from adjusted_identity import (
    score_alignment,
    align_and_score,
    DEFAULT_ADJUSTMENT_PARAMS,
    RAW_ADJUSTMENT_PARAMS,
    AdjustmentParams,
    _analyze_alignment,
)


class TestUnifiedMetrics:
    """Verify that both adjust_gaps modes produce identical metrics."""

    @pytest.mark.parametrize("seq1,seq2", [
        # Perfect match
        ("ATCG", "ATCG"),
        # Simple substitution
        ("ATCG", "AXCG"),
        # Simple insertion/deletion
        ("AT-CG", "ATXCG"),
        ("ATXCG", "AT-CG"),
        # Homopolymer length difference
        ("AAAA-TT", "AAA--TT"),
        ("AAA--TT", "AAAA-TT"),
        # Multiple homopolymer differences
        ("AAAA-TT-GGG", "AAA--TTGGGG"),
        # Opposite direction extensions
        ("TGC-C-TC", "TGCT--TC"),
        # Both sequences have different content in gap region
        ("AAA-GGG", "AAATGGG"),
        # Core mismatch scenario
        ("AAA-X-GGG", "AAAYY-GGG"),
        # Simple indel in middle of sequence
        ("ATCG-ATCG", "ATCGAATCG"),
        # Dual gaps
        ("AA--TT", "AA--TT"),
        # Complex variant with trailing gap (known limitation)
        # Biological interpretation: This is a pure homopolymer difference (7 T's vs 4 T's)
        # with identical GCA suffix. Ideally: identity=1.0, mismatches=0.
        # Actual: identity=10/11, mismatches=1. The trailing gap (A/-) at position 12 is
        # excluded from the scoring region (which requires both sequences to have content),
        # so the variant range is truncated to positions 9-11 ("TGC" vs "GCA") instead of
        # 9-12 ("TGCA" vs "GCA"). This causes core comparison of "GC" vs "GCA" = 1 edit.
        # We accept this compromise because scoring region truncation is essential for
        # handling the more common case of terminal differences from primer variation.
        ("GCATTTTTTTGCA", "GCAT-TTTTGCA-"),
    ])
    def test_metrics_identical_for_both_modes(self, seq1, seq2):
        """All metrics should be identical regardless of adjust_gaps mode."""
        result_false = score_alignment(seq1, seq2, adjust_gaps=False)
        result_true = score_alignment(seq1, seq2, adjust_gaps=True)

        assert result_false.identity == result_true.identity, \
            f"Identity differs: False={result_false.identity}, True={result_true.identity}"
        assert result_false.mismatches == result_true.mismatches, \
            f"Mismatches differ: False={result_false.mismatches}, True={result_true.mismatches}"
        assert result_false.scored_positions == result_true.scored_positions, \
            f"Scored positions differ: False={result_false.scored_positions}, True={result_true.scored_positions}"
        assert result_false.seq1_coverage == result_true.seq1_coverage, \
            f"seq1_coverage differs: False={result_false.seq1_coverage}, True={result_true.seq1_coverage}"
        assert result_false.seq2_coverage == result_true.seq2_coverage, \
            f"seq2_coverage differs: False={result_false.seq2_coverage}, True={result_true.seq2_coverage}"

    @pytest.mark.parametrize("normalize_hp", [True, False])
    @pytest.mark.parametrize("normalize_indels", [True, False])
    @pytest.mark.parametrize("handle_iupac", [True, False])
    def test_metrics_identical_with_various_params(self, normalize_hp, normalize_indels, handle_iupac):
        """Metrics should be identical with any parameter combination."""
        params = AdjustmentParams(
            normalize_homopolymers=normalize_hp,
            normalize_indels=normalize_indels,
            handle_iupac_overlap=handle_iupac,
            end_skip_distance=0,
        )

        seq1 = "AAAA-TGRC"
        seq2 = "AAA--TGNC"

        result_false = score_alignment(seq1, seq2, params, adjust_gaps=False)
        result_true = score_alignment(seq1, seq2, params, adjust_gaps=True)

        assert result_false.identity == result_true.identity
        assert result_false.mismatches == result_true.mismatches
        assert result_false.scored_positions == result_true.scored_positions


class TestOutputDifferences:
    """Verify that output strings differ as expected between modes."""

    def test_annotated_output_preserves_alignment(self):
        """adjust_gaps=False should return original alignment strings."""
        seq1 = "AAAA-TT"
        seq2 = "AAA--TT"

        result = score_alignment(seq1, seq2, adjust_gaps=False)

        assert result.seq1_aligned == seq1, "seq1_aligned should match input"
        assert result.seq2_aligned == seq2, "seq2_aligned should match input"

    def test_adjusted_output_may_differ(self):
        """adjust_gaps=True may produce different alignment strings."""
        seq1 = "AAAA-TT"
        seq2 = "AAA--TT"

        result_false = score_alignment(seq1, seq2, adjust_gaps=False)
        result_true = score_alignment(seq1, seq2, adjust_gaps=True)

        # Alignment strings may differ
        # (though in this simple case they might be the same)
        # The key point is that metrics are identical
        assert result_false.identity == result_true.identity

    def test_biological_content_preserved_in_both_modes(self):
        """Both modes should preserve biological (non-gap) content."""
        seq1 = "AAAA-TT"
        seq2 = "AAA--TT"

        result_false = score_alignment(seq1, seq2, adjust_gaps=False)
        result_true = score_alignment(seq1, seq2, adjust_gaps=True)

        # Extract non-gap content
        content_false_1 = result_false.seq1_aligned.replace('-', '')
        content_false_2 = result_false.seq2_aligned.replace('-', '')
        content_true_1 = result_true.seq1_aligned.replace('-', '')
        content_true_2 = result_true.seq2_aligned.replace('-', '')

        # Original content
        orig_seq1 = seq1.replace('-', '')
        orig_seq2 = seq2.replace('-', '')

        assert content_false_1 == orig_seq1
        assert content_false_2 == orig_seq2
        assert content_true_1 == orig_seq1
        assert content_true_2 == orig_seq2


class TestUnifiedAnalysisFunction:
    """Test the _analyze_alignment function directly."""

    def test_analysis_returns_correct_structure(self):
        """_analyze_alignment should return AlignmentAnalysis with all fields."""
        seq1 = "AAAA-TT"
        seq2 = "AAA--TT"

        analysis = _analyze_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)

        # Check required fields exist
        assert hasattr(analysis, 'identity')
        assert hasattr(analysis, 'mismatches')
        assert hasattr(analysis, 'scored_positions')
        assert hasattr(analysis, 'seq1_coverage')
        assert hasattr(analysis, 'seq2_coverage')
        assert hasattr(analysis, 'scoring_start')
        assert hasattr(analysis, 'scoring_end')
        assert hasattr(analysis, 'variant_ranges')

    def test_analysis_metrics_match_score_alignment(self):
        """Metrics from _analyze_alignment should match score_alignment output."""
        seq1 = "AAAA-TT"
        seq2 = "AAA--TT"

        analysis = _analyze_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)
        result = score_alignment(seq1, seq2, adjust_gaps=False)

        assert analysis.identity == result.identity
        assert analysis.mismatches == result.mismatches
        assert analysis.scored_positions == result.scored_positions
        assert analysis.seq1_coverage == result.seq1_coverage
        assert analysis.seq2_coverage == result.seq2_coverage

    def test_variant_ranges_captured(self):
        """Variant ranges should be captured in analysis."""
        seq1 = "ATCG-ATCG-ATCG"
        seq2 = "ATCGXATCGYATCG"

        analysis = _analyze_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)

        # Should have 2 variant ranges (at positions 4 and 9)
        assert len(analysis.variant_ranges) == 2

        # Check first variant range
        vr1 = analysis.variant_ranges[0]
        assert vr1.allele1 == ""  # gap in seq1
        assert vr1.allele2 == "X"

        # Check second variant range
        vr2 = analysis.variant_ranges[1]
        assert vr2.allele1 == ""  # gap in seq1
        assert vr2.allele2 == "Y"


class TestAlignAndScoreIntegration:
    """Test unified architecture with align_and_score function."""

    def test_align_and_score_unified_metrics(self):
        """align_and_score should produce identical metrics for both modes."""
        seq1 = "AAAATTTGGG"
        seq2 = "AAATTTGGG"

        result_false = align_and_score(seq1, seq2, adjust_gaps=False)
        result_true = align_and_score(seq1, seq2, adjust_gaps=True)

        assert result_false.identity == result_true.identity
        assert result_false.mismatches == result_true.mismatches
        assert result_false.scored_positions == result_true.scored_positions

    def test_real_sequence_comparison(self):
        """Test with realistic sequence comparison."""
        seq1 = "TTTTCACAGGCTGGTAATGGCT"
        seq2 = "TTTTCACAAGTTGGTAATGGCT"

        result_false = align_and_score(seq1, seq2, adjust_gaps=False)
        result_true = align_and_score(seq1, seq2, adjust_gaps=True)

        # Metrics should be identical
        assert result_false.identity == result_true.identity
        assert result_false.mismatches == result_true.mismatches
        assert result_false.scored_positions == result_true.scored_positions


class TestScoreStringGeneration:
    """Test score string generation for both modes."""

    def test_score_string_length_matches_alignment(self):
        """Score string length should match alignment length."""
        seq1 = "AAAA-TT"
        seq2 = "AAA--TT"

        result_false = score_alignment(seq1, seq2, adjust_gaps=False)
        result_true = score_alignment(seq1, seq2, adjust_gaps=True)

        # For annotated output
        assert len(result_false.score_aligned) == len(result_false.seq1_aligned)
        assert len(result_false.score_aligned) == len(result_false.seq2_aligned)

        # For adjusted output
        assert len(result_true.score_aligned) == len(result_true.seq1_aligned)
        assert len(result_true.score_aligned) == len(result_true.seq2_aligned)

    def test_match_markers_at_match_positions(self):
        """Match positions should have match markers."""
        seq1 = "ATCG"
        seq2 = "ATCG"

        result = score_alignment(seq1, seq2, adjust_gaps=False)

        # All positions should be matches
        assert result.score_aligned == "||||"

    def test_mismatch_markers_at_mismatch_positions(self):
        """Mismatch positions should have mismatch markers."""
        seq1 = "ATCG"
        seq2 = "AXCG"

        result = score_alignment(seq1, seq2, adjust_gaps=False)

        # Position 1 should be a mismatch (space)
        assert result.score_aligned[0] == '|'
        assert result.score_aligned[1] == ' '  # mismatch
        assert result.score_aligned[2] == '|'
        assert result.score_aligned[3] == '|'


class TestPerformanceConsiderations:
    """Verify performance characteristics of unified architecture."""

    def test_single_analysis_for_both_modes(self):
        """Both modes should use the same analysis pass internally."""
        # This test verifies the architecture by checking that both modes
        # produce identical metrics, which is only possible if they use
        # the same analysis logic.

        seq1 = "ATCG" * 50 + "AAA--" + "GCTA" * 50
        seq2 = "ATCG" * 50 + "AAAA-" + "GCTA" * 50

        result_false = score_alignment(seq1, seq2, adjust_gaps=False)
        result_true = score_alignment(seq1, seq2, adjust_gaps=True)

        # If both used the same analysis, metrics must be identical
        assert result_false.identity == result_true.identity
        assert result_false.mismatches == result_true.mismatches
        assert result_false.scored_positions == result_true.scored_positions

    def test_no_double_analysis_for_adjust_gaps_true(self):
        """adjust_gaps=True should not re-analyze the adjusted alignment."""
        # This is tested implicitly by verifying metrics are identical.
        # If adjust_gaps=True re-analyzed the adjusted alignment, metrics
        # would differ in some cases (like the old implementation).

        # This test case specifically triggered different metrics in the old code
        seq1 = "GCATTTTTTTGCA"
        seq2 = "GCAT-TTTTGCA-"

        result_false = score_alignment(seq1, seq2, adjust_gaps=False)
        result_true = score_alignment(seq1, seq2, adjust_gaps=True)

        # Metrics must be identical (proves single analysis pass)
        assert result_false.identity == result_true.identity
        assert result_false.mismatches == result_true.mismatches
        assert result_false.scored_positions == result_true.scored_positions


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
