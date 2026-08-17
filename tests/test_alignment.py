#!/usr/bin/env python3
"""
Tests for alignment functions (align_edlib_bidirectional, align_and_score).

These tests focus on the alignment algorithms and end-to-end functionality.
"""

import pytest
from adjusted_identity import (
    align_edlib_bidirectional,
    _align_nw,
    align_and_score,
    AdjustmentParams,
    ScoringMode,
    DEFAULT_ADJUSTMENT_PARAMS,
    RAW_ADJUSTMENT_PARAMS,
)


class TestAlignEdlibBidirectional:
    """Test the bidirectional alignment algorithm."""
    
    def test_simple_alignment(self):
        """Basic alignment should work correctly."""
        result = align_edlib_bidirectional("ATCG", "ATCG")
        assert result is not None
        assert result['aligned_seq1'] == "ATCG"
        assert result['aligned_seq2'] == "ATCG"
    
    def test_alignment_with_gaps(self):
        """Alignment with insertions/deletions."""
        result = align_edlib_bidirectional("ATCG", "ATCCG")
        assert result is not None
        # Should have alignment with gap or mismatch
        assert len(result['aligned_seq1']) == len(result['aligned_seq2'])
    
    def test_empty_sequence_handling(self):
        """Empty sequences should return None."""
        assert align_edlib_bidirectional("", "ATCG") is None
        assert align_edlib_bidirectional("ATCG", "") is None
        assert align_edlib_bidirectional("", "") is None
    
    def test_very_different_sequences(self):
        """Very different sequences should still align."""
        result = align_edlib_bidirectional("ATCG", "CGAT")
        assert result is not None
        assert len(result['aligned_seq1']) == len(result['aligned_seq2'])


class TestAlignAndScore:
    """Test the main align_and_score function."""
    
    def test_identical_sequences(self):
        """Identical sequences should have perfect identity."""
        result = align_and_score("ATCGATCG", "ATCGATCG")
        assert result.identity == 1.0
        assert result.mismatches == 0
        assert result.seq1_coverage == 1.0
        assert result.seq2_coverage == 1.0
    
    def test_sequences_with_mismatches(self):
        """Sequences with clear mismatches should have lower identity."""
        # Use sequences that align fully but have mismatches in the middle
        result = align_and_score("ATCGXXCG", "ATCGAACG")
        assert result.identity < 1.0  # Should not be perfect due to XX vs AA
        assert result.mismatches > 0
    
    def test_coverage_calculation(self):
        """Coverage should be calculated correctly."""
        # Sequences where not all of each sequence participates in alignment
        result = align_and_score("ATCGATCGATCG", "GATCG")
        assert 0.0 <= result.seq1_coverage <= 1.0
        assert 0.0 <= result.seq2_coverage <= 1.0
        # Shorter sequence should have higher coverage
        assert result.seq2_coverage >= result.seq1_coverage
    
    def test_empty_sequence_handling(self):
        """Empty sequences should be handled gracefully."""
        result = align_and_score("", "ATCG")
        assert result.identity == 0.0
        assert result.seq1_coverage == 0.0
        assert result.seq2_coverage == 0.0
        
        result = align_and_score("ATCG", "")
        assert result.identity == 0.0
        assert result.seq1_coverage == 0.0
        assert result.seq2_coverage == 0.0
    
    def test_parameter_passing(self):
        """Custom parameters should be passed through correctly."""
        seq1, seq2 = "AAATTT", "AAA-TTT"
        
        # With homopolymer adjustment
        adjusted = align_and_score(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)
        
        # Without homopolymer adjustment
        raw = align_and_score(seq1, seq2, RAW_ADJUSTMENT_PARAMS)
        
        # Results should be different
        assert adjusted.identity != raw.identity


class TestRealWorldScenarios:
    """Test scenarios based on real-world mycological sequence data."""
    
    def test_its_sequence_comparison(self):
        """Simulated ITS sequence comparison with common issues."""
        # Simulated ITS sequences with homopolymer differences and ambiguities
        its1 = "TCCGTAGGTGAACCTGCGGAAGGATCATTACCGAGTTTA"
        its2 = "TCCGTAGGTGAACCTGCGGAAGGATCATTACCGAGTTTTA"  # Extra T
        
        result = align_and_score(its1, its2, DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity > 0.95  # Should be very similar with adjustment
        assert result.seq1_coverage > 0.9
        assert result.seq2_coverage > 0.9
    
    def test_barcode_gap_sequence(self):
        """Test with barcode gap (sequence quality issues)."""
        # Common scenario: good sequence with a region of poor quality
        good_seq = "ATCGATCGATCGATCGATCG"
        gap_seq = "ATCGATCG--------GATCG"  # Poor quality region replaced with gaps
        
        result = align_and_score(good_seq, gap_seq)
        # Should still have reasonable identity due to adjustments
        assert result.identity > 0.5
        assert result.seq1_coverage > 0.5
        assert result.seq2_coverage > 0.5
    
    def test_primer_region_differences(self):
        """Test sequences with different primer regions (ends)."""
        # Common in amplicon sequencing: same gene, different primers
        seq1 = "PRIMER1" + "ATCGATCGATCGATCG" + "PRIMER2"
        seq2 = "XXXXXX" + "ATCGATCGATCGATCG" + "YYYYYY"  # Different primers
        
        # With end trimming, primer differences should be ignored
        result = align_and_score(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)
        # Core region identity should be preserved
        assert result.identity > 0.5


class TestPerformanceAndRobustness:
    """Test performance with various sequence characteristics."""
    
    def test_long_sequences(self):
        """Test with longer sequences (simulating gene sequences)."""
        # Create longer sequences
        seq1 = "ATCGATCGATCG" * 50  # 600bp
        seq2 = "ATCGATCGATCG" * 49 + "ATCGATCGATCC"  # One difference
        
        result = align_and_score(seq1, seq2)
        assert result.identity > 0.99  # Should be very high
        assert result.seq1_coverage > 0.95
        assert result.seq2_coverage > 0.95
    
    def test_highly_repetitive_sequences(self):
        """Test with highly repetitive sequences."""
        seq1 = "AT" * 100  # 200bp of AT repeats
        seq2 = "AT" * 99 + "AG"  # One difference at end
        
        result = align_and_score(seq1, seq2)
        assert result.identity > 0.95
        assert result is not None
    
    def test_low_complexity_sequences(self):
        """Test with low complexity sequences (homopolymers)."""
        seq1 = "A" * 50
        seq2 = "A" * 48  # Shorter homopolymer
        
        result = align_and_score(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)
        # With homopolymer adjustment, should be identical
        assert result.identity == 1.0
        
        result_raw = align_and_score(seq1, seq2, RAW_ADJUSTMENT_PARAMS)
        # With RAW_ADJUSTMENT_PARAMS (end_skip_distance=0), only overlap region is scored.
        # The 2 missing A's are in overhang region and ignored, so identity should be 1.0
        assert result_raw.identity == 1.0

    def test_alignment_length_mismatch_regression(self):
        """Regression test for alignment length mismatch bug.
        
        This test uses sequences extracted from real FASTA data that previously
        caused 'Aligned sequences must have same length' errors due to incorrect
        sequence trimming in align_edlib_bidirectional().
        """
        # Problematic pair 1: Long vs short sequence with complex indels
        seq1 = 'GTGGTTGTAGCTGGCCCCTATTAAGGGGATGTGCACACTGTCTCTTTCTCTTGCTTGTTTTTTTCATTCTTTCCACTTGTGCACTGCTTGTAGGCAGCCTGGCATTGTTCAGGTTGTCTATGATTTTCTTTACATACATGAATGATTGTTGTACAGAATGTGATGAAAAAAAAAGTAATACAACTTTCAACAACGGATCTCTTGGGTCTCGCATCGATGAAAAACGCACCGAAATGCGATAAGTAATGTGAATTGCGGAATTCTCTGAATCATCNAATCTTTGAACGCATCTTGCGCTCCTTGGGATTCCGAGGAGCATGCCTGTTTGAGTGTCATTAAATTCTGTCAAGACATGCACTTGAGTGTGTTTTGGATTGTGGGAGTGTCTGCTGGATTCTTTATATATATGAGCCAGCTCTCCTGAAAGACATTAGCTTTGGAGGGATGTGCCAAGTCACTTCTGCCTTTCCATTGGTGTGATAGATGAATAAACTTATCTACGCCAGGAAAGCAGGTTGCAGGTGATGCACTATGATCTCTCTGCTCTCTAATTGACATTTGTCTGATAACTTGACCTCAAATCAGGTAGGACTACCCGCTGAACTTAAGCATATCAATAAGCGGAGGAAAAGAAACTAACAAGGATTCCCCTAGTAACTGCGAGTGAAGAGGGAAGAGCTCAAATTTGAAATCTGGCAGTGTTTTGCTGTCCGAGTTGTAATCTAGAGAAGTGCTGCCCGTGCTGGACCATGTACAAGTCTCCTGGAATGGAGCGTCACAGANGGTGAGAATCCCGTCTTTGACATGGACTGCCAGTGCATTTGTGGTGTGCTCTCAAAGAGTCGAGTTGTTTGGGAATGCAGCTCTAAATGGGTGGTAAACTCCATCTAAAGCTAAATA'
        seq2 = 'GTGACCTGCGGAAGGANCATTATTGAATTTTATAAAGACAACTAGTAGGGAGTCTGTTGCTGGCTCCTCTTTGGAGGCGCATGTGCACGTCTTTTTT'
        
        # This should now work without error
        result = align_and_score(seq1, seq2)
        assert result is not None
        assert 0.0 <= result.identity <= 1.0
        
        # Problematic pair 2: Sequences with IUPAC codes and complex alignment
        seq1 = 'TMASYMTTMTTCGTMAGAWGAACTGCGGAAGGATCATTATTGAAGGAGAATGGGTGGCAAGGGCTGTTGCTGGCTTGAATGAGCATGTGCACGTCTGTTGCTRCTTATTTCATTCATWTTTCCTCCTGTGCAYGTTTTGTAGACACTTGGGAATGAGAGGTTGGTTGTAATGTAATGAATTGACCTCTTGAGGTSARTCTGGGWGTCTATGACMTTTTTATWAACMCSGCTGSWTGTGTATGGAATGAGAYTGKAGGTTTTTAATKAAAAAMCCTKTWAARRKAAAAARYAACACTTTMCAACACGGATCTTGTGGCTCTCSYCTCTRAAAAAAACSCCGCAAAATGSTAAAAAATGTGKAAATTGGAAAATTCWGWGAAATCTCAAATTTTTAAACSCCTTGTGCTCCCCGTGGTTTTCAGAAGAAYGYGCCTGTTTGAGTGWTTTTCRAATCTCTCAAAATGTTWTGTGCWATTMTTGGWGCWTGGGAATTTTGGAAGTTGGGGGTTGCTGGTCAAKTGKAAAGSTGTTCGGTTCTATAAAAAASSATGAGCTKGGGGCTCTCTRCWCKCGTGATGGGKGATCTACGCTCTGAGACGWGWGATGAGGTGTCTGCTGTCTRCTGATCCTCAGTGCMCAAKATGATAAACTTGACATCTGATCTCGTATGACTACRCGCTGCMCTYWAGCATATCATATTMYGARGAARRRRGRAAAAAAA'
        seq2 = 'ACCTGCGGAAGGATCATTATTGAAGGAGAATGGGTGGCAAGGGCTGTTGCTGGCTTGAATGAGCATGTGCACGTCTGTTGCTRCTTATTTCATTCATWTTTCCTCCTGTGCAYGTTTTGTAGACACTTGGGAATGAGAGGTTGGTTGTAATG'

        # This should now work without error
        result = align_and_score(seq1, seq2)
        assert result is not None
        assert 0.0 <= result.identity <= 1.0


class TestNWAlignment:
    """Test the Needleman-Wunsch alignment wrapper."""

    def test_identical_sequences(self):
        """NW alignment of identical sequences."""
        result = _align_nw("ATCG", "ATCG")
        assert result is not None
        assert result['aligned_seq1'] == "ATCG"
        assert result['aligned_seq2'] == "ATCG"

    def test_different_lengths(self):
        """NW alignment with different lengths produces gaps."""
        result = _align_nw("ATCG", "ATCGATCG")
        assert result is not None
        assert len(result['aligned_seq1']) == len(result['aligned_seq2'])
        # Shorter sequence should have gaps
        assert '-' in result['aligned_seq1']

    def test_empty_sequences_return_none(self):
        """Empty sequences should return None."""
        assert _align_nw("", "ATCG") is None
        assert _align_nw("ATCG", "") is None
        assert _align_nw("", "") is None

    def test_with_mismatches(self):
        """NW alignment with mismatches."""
        result = _align_nw("ATCG", "AXCG")
        assert result is not None
        assert len(result['aligned_seq1']) == len(result['aligned_seq2'])


class TestGlobalModeEndToEnd:
    """End-to-end tests for GLOBAL scoring mode via align_and_score."""

    def test_identical_sequences_global(self):
        """Identical sequences should have perfect identity in GLOBAL mode."""
        params = AdjustmentParams(scoring_mode=ScoringMode.GLOBAL)
        result = align_and_score("ATCGATCG", "ATCGATCG", params)
        assert result.identity == 1.0
        assert result.mismatches == 0

    def test_asymmetric_length_global_vs_local(self):
        """Key scenario: short vs long sequence — GLOBAL penalizes, LOCAL doesn't."""
        short = "ATCG" * 55   # 220bp
        long = "ATCG" * 165   # 660bp

        local_result = align_and_score(short, long)
        global_result = align_and_score(
            short, long,
            AdjustmentParams(scoring_mode=ScoringMode.GLOBAL, normalize_indels=False)
        )

        # LOCAL should report ~100% identity (overlap region matches)
        assert local_result.identity > 0.95

        # GLOBAL should report ~33% identity (terminal gaps penalized)
        assert global_result.identity < 0.5
        assert global_result.identity > 0.2

        # GLOBAL mode reports coverage based on overlap within the full NW alignment
        assert global_result.seq2_coverage < 0.5

    def test_global_empty_sequence(self):
        """Empty sequence handling in GLOBAL mode."""
        params = AdjustmentParams(scoring_mode=ScoringMode.GLOBAL)
        result = align_and_score("", "ATCG", params)
        assert result.identity == 0.0

class TestFlankReattachment:
    """Sequence regions beyond the overlap must be preserved in the alignment.

    Regression tests: edlib HW alignment omits target (seq2) regions outside
    the aligned location range from getNiceAlignment output. Previously these
    flanks were silently dropped from seq2_aligned, so seq2_coverage read 1.0
    no matter how much of seq2 was unaligned, and coverage was asymmetric
    between align_and_score(a, b) and align_and_score(b, a).
    """

    def test_seq2_tail_overhang(self):
        """seq2 extending beyond seq1's 3' end must appear in the alignment."""
        result = align_and_score("AAATTTGGG", "AAATTTGGGCCCCC")
        assert result.seq1_aligned == "AAATTTGGG-----"
        assert result.seq2_aligned == "AAATTTGGGCCCCC"
        assert result.identity == 1.0
        assert result.seq1_coverage == 1.0
        assert result.seq2_coverage == pytest.approx(9 / 14)

    def test_seq2_head_overhang(self):
        """seq2 extending beyond seq1's 5' end must appear in the alignment."""
        result = align_and_score("AAATTTGGG", "CCCCCAAATTTGGG")
        assert result.seq1_aligned == "-----AAATTTGGG"
        assert result.seq2_aligned == "CCCCCAAATTTGGG"
        assert result.identity == 1.0
        assert result.seq1_coverage == 1.0
        assert result.seq2_coverage == pytest.approx(9 / 14)

    def test_staggered_overhangs(self):
        """Both sequences overhanging opposite ends must both be preserved."""
        result = align_and_score("CCCCCAAATTTGGG", "AAATTTGGGAAAAA")
        assert result.seq1_aligned == "CCCCCAAATTTGGG-----"
        assert result.seq2_aligned == "-----AAATTTGGGAAAAA"
        assert result.identity == 1.0
        assert result.seq1_coverage == pytest.approx(9 / 14)
        assert result.seq2_coverage == pytest.approx(9 / 14)

    def test_coverage_symmetric_under_argument_swap(self):
        """Coverage must be consistent whichever argument has the overhang."""
        pairs = [
            ("AAATTTGGG", "AAATTTGGGCCCCC"),
            ("AAATTTGGG", "CCCCCAAATTTGGG"),
            ("CCCCCAAATTTGGG", "AAATTTGGGAAAAA"),
        ]
        for seq_a, seq_b in pairs:
            forward = align_and_score(seq_a, seq_b)
            reverse = align_and_score(seq_b, seq_a)
            assert forward.seq1_coverage == pytest.approx(reverse.seq2_coverage)
            assert forward.seq2_coverage == pytest.approx(reverse.seq1_coverage)
            assert forward.identity == pytest.approx(reverse.identity)

    def test_biological_content_preserved_with_overhangs(self):
        """Aligned strings must contain every input nucleotide."""
        seq1, seq2 = "CCCCCAAATTTGGG", "AAATTTGGGAAAAA"
        result = align_and_score(seq1, seq2)
        assert result.seq1_aligned.replace('-', '') == seq1
        assert result.seq2_aligned.replace('-', '') == seq2


class TestDissimilarSequenceFailure:
    """Completely dissimilar sequences must fail gracefully, not raise.

    Regression test: when the RC-stage CIGAR reported the entire query as a
    terminal insertion, trimming emptied one sequence and the forward
    alignment crashed inside edlib.getNiceAlignment.
    """

    def test_align_and_score_returns_failure_sentinel(self):
        result = align_and_score("AAAAAAAAAA", "TTTTTTTTTT")
        assert result.identity == 0.0
        assert result.mismatches == -1  # alignment-failure sentinel
        assert result.scored_positions == 0
        assert result.seq1_aligned == ""
        assert result.seq2_aligned == ""

    def test_bidirectional_returns_none(self):
        assert align_edlib_bidirectional("AAAAAAAAAA", "TTTTTTTTTT") is None
