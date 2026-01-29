"""
Tests for the adjust_gaps parameter in score_alignment and align_and_score.

These tests verify that:
1. adjust_gaps=True produces identical metrics to adjust_gaps=False
2. Adjusted scoring strings are position-by-position interpretable
3. Specific scenarios like opposite direction extensions work correctly
4. Backward compatibility is maintained (adjust_gaps=False is default)
"""

import pytest
from adjusted_identity import (
    score_alignment,
    align_and_score,
    DEFAULT_ADJUSTMENT_PARAMS,
    RAW_ADJUSTMENT_PARAMS,
    AdjustmentParams,
    AlignmentResult,
)


class TestMetricsPreservation:
    """Verify adjust_gaps=True produces identical metrics to adjust_gaps=False."""

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
        # Longer sequences with proper alignment
        ("ATCGATCG----ATCG", "ATCGATCGATCG----"),
        # Opposite direction extensions (key test case from plan)
        ("TGC-C-TC", "TGCT--TC"),
        # Both sequences have different content in gap region
        ("AAA-GGG", "AAATGGG"),
        # Core mismatch scenario
        ("AAA-X-GGG", "AAAYY-GGG"),
        # Simple indel in middle of sequence
        ("ATCG-ATCG", "ATCGAATCG"),
    ])
    def test_identity_metrics_identical(self, seq1, seq2):
        """Identity value should be identical with or without gap adjustment.

        Note: Some pathological alignments where gap adjustment relocates gaps
        may have different scored_positions due to changes in overlap region.
        We test realistic alignment patterns that a real aligner would produce.
        """
        result_normal = score_alignment(seq1, seq2, adjust_gaps=False)
        result_adjusted = score_alignment(seq1, seq2, adjust_gaps=True)

        assert result_normal.identity == result_adjusted.identity, \
            f"Identity differs: normal={result_normal.identity}, adjusted={result_adjusted.identity}"

    @pytest.mark.parametrize("seq1,seq2", [
        ("ATCG", "ATCG"),
        ("ATCG", "AXCG"),
        ("AAAA-TT", "AAA--TT"),
        ("TGC-C-TC", "TGCT--TC"),
        ("AAA-X-GGG", "AAAYY-GGG"),
    ])
    def test_mismatches_identical(self, seq1, seq2):
        """Mismatch count should be identical with or without gap adjustment."""
        result_normal = score_alignment(seq1, seq2, adjust_gaps=False)
        result_adjusted = score_alignment(seq1, seq2, adjust_gaps=True)

        assert result_normal.mismatches == result_adjusted.mismatches, \
            f"Mismatches differ: normal={result_normal.mismatches}, adjusted={result_adjusted.mismatches}"

    @pytest.mark.parametrize("seq1,seq2", [
        ("ATCG", "ATCG"),
        ("AAAA-TT", "AAA--TT"),
        ("TGC-C-TC", "TGCT--TC"),
    ])
    def test_coverage_identical(self, seq1, seq2):
        """Coverage values should be identical with or without gap adjustment."""
        result_normal = score_alignment(seq1, seq2, adjust_gaps=False)
        result_adjusted = score_alignment(seq1, seq2, adjust_gaps=True)

        assert result_normal.seq1_coverage == result_adjusted.seq1_coverage
        assert result_normal.seq2_coverage == result_adjusted.seq2_coverage


class TestScoringStringOptimality:
    """Verify adjusted scoring strings are position-by-position interpretable."""

    def test_extension_markers_have_one_gap(self):
        """Extension positions should have exactly one gap and one nucleotide."""
        # Homopolymer difference: one A is extension
        result = score_alignment("AAAA-TT", "AAA--TT", adjust_gaps=True)

        for i, marker in enumerate(result.score_aligned):
            if marker == '=':  # Extension marker
                c1, c2 = result.seq1_aligned[i], result.seq2_aligned[i]
                gaps = (c1 == '-') + (c2 == '-')
                assert gaps == 1, \
                    f"Position {i}: extension marker '=' should have 1 gap, got {gaps} ({c1}/{c2})"

    def test_match_markers_have_no_gaps(self):
        """Match positions should have content in both sequences."""
        result = score_alignment("ATCGATCG", "ATCGATCG", adjust_gaps=True)

        for i, marker in enumerate(result.score_aligned):
            if marker == '|':
                assert result.seq1_aligned[i] != '-', \
                    f"Position {i}: match marker '|' but seq1 has gap"
                assert result.seq2_aligned[i] != '-', \
                    f"Position {i}: match marker '|' but seq2 has gap"

    def test_no_dual_gaps_in_adjusted_scoring_region(self):
        """Adjusted alignment should not contain dual-gaps in the scoring region."""
        # Test with a case that might have dual-gaps in original
        result = score_alignment("AA--TT", "AA--TT", adjust_gaps=True)

        # Find scoring region (skip end-trimmed positions)
        for i in range(len(result.seq1_aligned)):
            if result.score_aligned[i] not in ('.', ):  # Not end-trimmed or dual-gap
                if result.seq1_aligned[i] == '-' and result.seq2_aligned[i] == '-':
                    # Dual-gaps should be marked with '.' in score_aligned
                    assert result.score_aligned[i] == '.', \
                        f"Dual-gap at position {i} not properly marked"

    def test_adjusted_alignment_same_length_components(self):
        """Adjusted seq1_aligned, seq2_aligned, and score_aligned should all be same length."""
        test_cases = [
            ("AAAA-TT", "AAA--TT"),
            ("TGC-C-TC", "TGCT--TC"),
            ("ATCG", "ATCG"),
        ]

        for seq1, seq2 in test_cases:
            result = score_alignment(seq1, seq2, adjust_gaps=True)
            assert len(result.seq1_aligned) == len(result.seq2_aligned), \
                f"seq1_aligned and seq2_aligned lengths differ for {seq1} vs {seq2}"
            assert len(result.seq1_aligned) == len(result.score_aligned), \
                f"seq1_aligned and score_aligned lengths differ for {seq1} vs {seq2}"


class TestGapAdjustmentScenarios:
    """Test specific gap adjustment scenarios."""

    def test_homopolymer_length_difference(self):
        """Different homopolymer lengths should result in extensions."""
        result = score_alignment("AAAA-TT", "AAA--TT", adjust_gaps=True)

        # Should have perfect identity (homopolymer differences ignored)
        assert result.identity == 1.0

        # Extra A should be shown as extension
        assert '=' in result.score_aligned

    def test_simple_substitution_unchanged(self):
        """Simple substitutions should be unaffected by gap adjustment."""
        seq1 = "ATCG"
        seq2 = "AXCG"

        result_normal = score_alignment(seq1, seq2, adjust_gaps=False)
        result_adjusted = score_alignment(seq1, seq2, adjust_gaps=True)

        # Aligned sequences should be identical (no gaps to adjust)
        assert result_adjusted.seq1_aligned == result_normal.seq1_aligned
        assert result_adjusted.seq2_aligned == result_normal.seq2_aligned
        assert result_adjusted.score_aligned == result_normal.score_aligned

    def test_perfect_match_unchanged(self):
        """Perfect matches should be unaffected by gap adjustment."""
        seq1 = "ATCGATCG"
        seq2 = "ATCGATCG"

        result_normal = score_alignment(seq1, seq2, adjust_gaps=False)
        result_adjusted = score_alignment(seq1, seq2, adjust_gaps=True)

        assert result_adjusted.seq1_aligned == result_normal.seq1_aligned
        assert result_adjusted.seq2_aligned == result_normal.seq2_aligned

    def test_extension_visualization(self):
        """Extensions should be visualized with extension marker."""
        # Simple homopolymer case
        result = score_alignment("AAA-TT", "AA--TT", adjust_gaps=True)

        # The extra A should be shown as extension
        extension_count = result.score_aligned.count('=')
        assert extension_count >= 1, \
            f"Expected at least one extension marker, got {extension_count}"

    def test_gap_adjustment_preserves_biological_content(self):
        """Gap adjustment should preserve the biological (non-gap) content."""
        seq1 = "AAAA-TT"
        seq2 = "AAA--TT"

        result = score_alignment(seq1, seq2, adjust_gaps=True)

        # Extract non-gap content from adjusted alignment
        adj_seq1_content = result.seq1_aligned.replace('-', '')
        adj_seq2_content = result.seq2_aligned.replace('-', '')

        # Extract non-gap content from original
        orig_seq1_content = seq1.replace('-', '')
        orig_seq2_content = seq2.replace('-', '')

        assert adj_seq1_content == orig_seq1_content, \
            f"Biological content changed: {orig_seq1_content} -> {adj_seq1_content}"
        assert adj_seq2_content == orig_seq2_content, \
            f"Biological content changed: {orig_seq2_content} -> {adj_seq2_content}"


class TestBackwardCompatibility:
    """Verify adjust_gaps=False (default) preserves exact existing behavior."""

    def test_default_is_false(self):
        """Default value for adjust_gaps should be False."""
        seq1 = "AAAA-TT"
        seq2 = "AAA--TT"

        # Call without adjust_gaps parameter
        result_default = score_alignment(seq1, seq2)
        result_explicit_false = score_alignment(seq1, seq2, adjust_gaps=False)

        # Should be identical
        assert result_default.seq1_aligned == result_explicit_false.seq1_aligned
        assert result_default.seq2_aligned == result_explicit_false.seq2_aligned
        assert result_default.score_aligned == result_explicit_false.score_aligned
        assert result_default.identity == result_explicit_false.identity

    def test_original_alignment_preserved_when_false(self):
        """When adjust_gaps=False, input alignment should be preserved."""
        seq1 = "AAAA-TT"
        seq2 = "AAA--TT"

        result = score_alignment(seq1, seq2, adjust_gaps=False)

        # Output should match input exactly
        assert result.seq1_aligned == seq1
        assert result.seq2_aligned == seq2

    def test_explicit_false_same_as_omitted(self):
        """Explicitly passing adjust_gaps=False should match default."""
        seq1 = "TGC-C-TC"
        seq2 = "TGCT--TC"

        result1 = score_alignment(seq1, seq2)
        result2 = score_alignment(seq1, seq2, adjust_gaps=False)

        assert result1 == result2


class TestAlignAndScoreIntegration:
    """Test adjust_gaps with align_and_score function."""

    def test_align_and_score_accepts_adjust_gaps(self):
        """align_and_score should accept adjust_gaps parameter."""
        seq1 = "AAAATTTGGG"
        seq2 = "AAATTTGGG"

        # Should not raise
        result_false = align_and_score(seq1, seq2, adjust_gaps=False)
        result_true = align_and_score(seq1, seq2, adjust_gaps=True)

        # Metrics should be identical
        assert result_false.identity == result_true.identity

    def test_align_and_score_default_is_false(self):
        """align_and_score default for adjust_gaps should be False."""
        seq1 = "AAAATTTGGG"
        seq2 = "AAATTTGGG"

        result_default = align_and_score(seq1, seq2)
        result_false = align_and_score(seq1, seq2, adjust_gaps=False)

        assert result_default.seq1_aligned == result_false.seq1_aligned
        assert result_default.seq2_aligned == result_false.seq2_aligned


class TestEdgeCases:
    """Test edge cases for gap adjustment."""

    def test_empty_variant_range(self):
        """Handle alignments with no variant ranges."""
        seq1 = "ATCGATCG"
        seq2 = "ATCGATCG"

        # Should work without error
        result = score_alignment(seq1, seq2, adjust_gaps=True)
        assert result.identity == 1.0

    def test_all_gaps_one_sequence(self):
        """Handle case where one part of alignment is all gaps."""
        seq1 = "ATCG----"
        seq2 = "----ATCG"

        result = score_alignment(seq1, seq2, adjust_gaps=True)
        # Should complete without error
        assert result is not None

    def test_single_position_variant(self):
        """Handle single-position variants."""
        seq1 = "ATG"
        seq2 = "ACG"

        result_normal = score_alignment(seq1, seq2, adjust_gaps=False)
        result_adjusted = score_alignment(seq1, seq2, adjust_gaps=True)

        assert result_normal.identity == result_adjusted.identity

    def test_with_raw_params(self):
        """Gap adjustment should work with RAW_ADJUSTMENT_PARAMS."""
        seq1 = "AAAA-TT"
        seq2 = "AAA--TT"

        # Should work without error
        result = score_alignment(seq1, seq2, RAW_ADJUSTMENT_PARAMS, adjust_gaps=True)
        assert result is not None

    def test_with_custom_params(self):
        """Gap adjustment should work with custom AdjustmentParams."""
        params = AdjustmentParams(
            normalize_homopolymers=True,
            handle_iupac_overlap=False,
            normalize_indels=False,
            end_skip_distance=0
        )

        seq1 = "AAAA-TT"
        seq2 = "AAA--TT"

        result = score_alignment(seq1, seq2, params, adjust_gaps=True)
        assert result is not None


class TestPerformanceRegression:
    """Ensure adjust_gaps doesn't significantly regress performance."""

    def test_default_path_minimal_overhead(self):
        """The default (adjust_gaps=False) path should have minimal overhead."""
        import time

        # Generate test sequences with some variant regions
        seq1 = "ATCG" * 50 + "AAA--" + "GCTA" * 50
        seq2 = "ATCG" * 50 + "AAAA-" + "GCTA" * 50

        # Time the default path
        start = time.perf_counter()
        for _ in range(100):
            score_alignment(seq1, seq2, adjust_gaps=False)
        default_time = time.perf_counter() - start

        # The default path should complete in reasonable time
        assert default_time < 5.0, f"Default path took {default_time:.2f}s for 100 iterations"

    def test_adjust_gaps_reasonable_overhead(self):
        """The adjust_gaps=True path should complete in reasonable time."""
        import time

        seq1 = "ATCG" * 50 + "AAA--" + "GCTA" * 50
        seq2 = "ATCG" * 50 + "AAAA-" + "GCTA" * 50

        start = time.perf_counter()
        for _ in range(100):
            score_alignment(seq1, seq2, adjust_gaps=True)
        adjusted_time = time.perf_counter() - start

        # Adjusted path can be slower, but should complete reasonably
        assert adjusted_time < 10.0, f"Adjusted path took {adjusted_time:.2f}s for 100 iterations"

    def test_overhead_ratio(self):
        """Compare overhead between adjust_gaps=True and adjust_gaps=False."""
        import time

        seq1 = "ATCG" * 50 + "AAA--" + "GCTA" * 50
        seq2 = "ATCG" * 50 + "AAAA-" + "GCTA" * 50

        # Time default path
        start = time.perf_counter()
        for _ in range(100):
            score_alignment(seq1, seq2, adjust_gaps=False)
        default_time = time.perf_counter() - start

        # Time adjusted path
        start = time.perf_counter()
        for _ in range(100):
            score_alignment(seq1, seq2, adjust_gaps=True)
        adjusted_time = time.perf_counter() - start

        # Report the overhead ratio (informational, not a strict assertion)
        if default_time > 0:
            overhead_ratio = adjusted_time / default_time
            print(f"\nPerformance: default={default_time:.3f}s, adjusted={adjusted_time:.3f}s, ratio={overhead_ratio:.2f}x")

        # Adjusted path should be no more than 5x slower (generous allowance)
        assert adjusted_time < default_time * 5, \
            f"Adjusted path is {adjusted_time/default_time:.1f}x slower than default"


class TestRealWorldScenarios:
    """Test real-world-like scenarios from mycological sequencing."""

    def test_its_like_homopolymer_variation(self):
        """Simulate ITS-like homopolymer variation."""
        # Common pattern: poly-T tract varies in length
        # seq1: GCA + 7T's + GCA (13 chars)
        # seq2: GCA + T + gap + 4T's + GCA + gap (alignment representation)
        seq1 = "GCATTTTTTTGCA"
        seq2 = "GCAT-TTTTGCA-"

        result = score_alignment(seq1, seq2, adjust_gaps=True)

        # The variant range at position 4-4 (extra T in seq1) is a pure extension
        # But the variant range at 9-11 has cores "GC" vs "GCA" that differ
        # So there is 1 mismatch due to the length difference in cores
        # With unified analysis, both adjust_gaps modes produce the same metrics
        assert result.mismatches == 1
        assert result.scored_positions == 11
        # Identity = 10/11 = 0.909...
        assert result.identity == pytest.approx(10/11, rel=1e-6)

    def test_multiple_variant_ranges(self):
        """Test alignment with multiple separate variant ranges."""
        seq1 = "ATCG-ATCG-ATCG"
        seq2 = "ATCGXATCGYATCG"

        result_normal = score_alignment(seq1, seq2, adjust_gaps=False)
        result_adjusted = score_alignment(seq1, seq2, adjust_gaps=True)

        # Metrics should match
        assert result_normal.identity == result_adjusted.identity
        assert result_normal.mismatches == result_adjusted.mismatches

    def test_user_provided_sequences_issue(self):
        """Test user-provided sequences with complex variant region.

        These sequences have multiple differences in a short region:
        - Homopolymer variation (A vs AA)
        - Homopolymer variation (GG vs G)
        - Substitution (C vs T)

        The gap adjustment should preserve metrics while providing
        cleaner visualization.
        """
        seq1 = "TTTTCACAGGCTGGTAATGGCT"
        seq2 = "TTTTCACAAGTTGGTAATGGCT"

        # Use align_and_score since these are unaligned sequences
        result_normal = align_and_score(seq1, seq2, adjust_gaps=False)
        result_adjusted = align_and_score(seq1, seq2, adjust_gaps=True)

        # Metrics should be identical
        assert result_normal.identity == result_adjusted.identity, \
            f"Identity differs: normal={result_normal.identity}, adjusted={result_adjusted.identity}"
        assert result_normal.mismatches == result_adjusted.mismatches, \
            f"Mismatches differ: normal={result_normal.mismatches}, adjusted={result_adjusted.mismatches}"
        assert result_normal.scored_positions == result_adjusted.scored_positions, \
            f"Scored positions differ: normal={result_normal.scored_positions}, adjusted={result_adjusted.scored_positions}"

        # Verify biological content is preserved
        adj_seq1_content = result_adjusted.seq1_aligned.replace('-', '')
        adj_seq2_content = result_adjusted.seq2_aligned.replace('-', '')
        assert adj_seq1_content == seq1, "seq1 biological content changed"
        assert adj_seq2_content == seq2, "seq2 biological content changed"

        # The sequences have one true mismatch (C vs T at position ~12)
        assert result_normal.mismatches == 1
        # Identity should be ~95% (20/21 positions match after adjustments)
        assert result_normal.identity == pytest.approx(20/21, rel=1e-4)

    def test_user_provided_sequences_issue_2(self):
        """Test user-provided sequences with extensions and core mismatch.

        seq1: CACAGGCTGGTAAT
        seq2: CACAAATTGGTAAT

        The variant region contains:
        - AA extension in seq2 (extends left A context)
        - GGC vs empty core in seq1 (after extensions consumed)
        - T extension in seq2 (extends right T context)
        """
        seq1 = "CACAGGCTGGTAAT"
        seq2 = "CACAAATTGGTAAT"

        result_normal = align_and_score(seq1, seq2, adjust_gaps=False)
        result_adjusted = align_and_score(seq1, seq2, adjust_gaps=True)

        # Metrics should be identical
        assert result_normal.identity == result_adjusted.identity, \
            f"Identity differs: normal={result_normal.identity}, adjusted={result_adjusted.identity}"
        assert result_normal.mismatches == result_adjusted.mismatches, \
            f"Mismatches differ: normal={result_normal.mismatches}, adjusted={result_adjusted.mismatches}"
        assert result_normal.scored_positions == result_adjusted.scored_positions, \
            f"Scored positions differ: normal={result_normal.scored_positions}, adjusted={result_adjusted.scored_positions}"

        # Verify biological content is preserved
        adj_seq1_content = result_adjusted.seq1_aligned.replace('-', '')
        adj_seq2_content = result_adjusted.seq2_aligned.replace('-', '')
        assert adj_seq1_content == seq1, "seq1 biological content changed"
        assert adj_seq2_content == seq2, "seq2 biological content changed"

        # The sequences have one mismatch (core content differs)
        assert result_normal.mismatches == 1
        # Identity should be ~91.67% (11/12 scored positions)
        assert result_normal.identity == pytest.approx(11/12, rel=1e-4)
