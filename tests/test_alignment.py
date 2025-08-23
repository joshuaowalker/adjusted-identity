#!/usr/bin/env python3
"""
Tests for alignment functions (align_edlib_bidirectional, align_and_score).

These tests focus on the alignment algorithms and end-to-end functionality.
"""

import pytest
from adjusted_identity import (
    align_edlib_bidirectional,
    align_and_score,
    AdjustmentParams,
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
    
    def test_completely_different_sequences(self):
        """Completely different sequences should have low identity."""
        result = align_and_score("ATCG", "CGAT")
        assert result.identity < 1.0  # Should not be perfect
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
        # Without adjustment, should show difference
        assert result_raw.identity < 1.0