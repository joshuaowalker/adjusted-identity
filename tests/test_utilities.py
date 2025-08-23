#!/usr/bin/env python3
"""
Tests for utility functions and helper classes.

Tests the various helper functions, data classes, and internal utilities.
"""

import pytest
from adjusted_identity import (
    AdjustmentParams,
    ScoringFormat,
    AlignmentResult,
    IUPAC_CODES,
    _are_nucleotides_equivalent,
    _parse_suffix_gap_from_cigar,
    _find_scoring_region,
)


class TestDataClasses:
    """Test the dataclasses used in the package."""
    
    def test_adjustment_params_defaults(self):
        """Test default AdjustmentParams values."""
        params = AdjustmentParams()
        assert params.normalize_homopolymers is True
        assert params.handle_iupac_overlap is True
        assert params.normalize_indels is True
        assert params.end_skip_distance == 20
    
    def test_adjustment_params_custom(self):
        """Test custom AdjustmentParams values."""
        params = AdjustmentParams(
            normalize_homopolymers=False,
            handle_iupac_overlap=False,
            normalize_indels=False,
            end_skip_distance=0
        )
        assert params.normalize_homopolymers is False
        assert params.handle_iupac_overlap is False
        assert params.normalize_indels is False
        assert params.end_skip_distance == 0
    
    def test_scoring_format_defaults(self):
        """Test default ScoringFormat values."""
        fmt = ScoringFormat()
        assert fmt.match == '|'
        assert fmt.substitution == ' '
        assert fmt.indel_start == ' '
        assert fmt.indel_extension == '-'
        assert fmt.homopolymer_extension == '='
        assert fmt.end_trimmed == '.'
    
    def test_scoring_format_validation(self):
        """Test ScoringFormat validation."""
        # Valid single characters should work
        fmt = ScoringFormat(match='*', substitution='X')
        assert fmt.match == '*'
        assert fmt.substitution == 'X'
        
        # Invalid (non-single character) should raise error
        with pytest.raises(ValueError, match="single character"):
            ScoringFormat(match="too_long")
        
        with pytest.raises(ValueError, match="single character"):
            ScoringFormat(substitution="")
    
    def test_alignment_result_immutable(self):
        """Test that AlignmentResult is immutable (frozen)."""
        result = AlignmentResult(
            identity=0.5,
            mismatches=2,
            scored_positions=4,
            seq1_coverage=0.8,
            seq2_coverage=0.9,
            seq1_aligned="ATCG",
            seq2_aligned="ATCG",
            score_aligned="||||"
        )
        
        # Should not be able to modify fields
        with pytest.raises(Exception):  # FrozenInstanceError or similar
            result.identity = 0.6


class TestNucleotideEquivalence:
    """Test IUPAC nucleotide equivalence function."""
    
    def test_exact_matches(self):
        """Test exact nucleotide matches."""
        assert _are_nucleotides_equivalent('A', 'A')
        assert _are_nucleotides_equivalent('T', 'T')
        assert _are_nucleotides_equivalent('C', 'C')
        assert _are_nucleotides_equivalent('G', 'G')
        assert _are_nucleotides_equivalent('N', 'N')
    
    def test_case_insensitive(self):
        """Test case insensitive matching."""
        assert _are_nucleotides_equivalent('a', 'A')
        assert _are_nucleotides_equivalent('A', 'a')
        assert _are_nucleotides_equivalent('r', 'R')
    
    def test_iupac_intersections_enabled(self):
        """Test IUPAC code intersections when enabled."""
        # R (AG) and K (GT) both contain G
        assert _are_nucleotides_equivalent('R', 'K', enable_iupac_intersection=True)
        
        # R (AG) contains A
        assert _are_nucleotides_equivalent('R', 'A', enable_iupac_intersection=True)
        
        # Y (CT) and S (GC) both contain C
        assert _are_nucleotides_equivalent('Y', 'S', enable_iupac_intersection=True)
        
        # R (AG) and Y (CT) have no overlap
        assert not _are_nucleotides_equivalent('R', 'Y', enable_iupac_intersection=True)
    
    def test_iupac_intersections_disabled(self):
        """Test IUPAC code handling when intersections disabled."""
        # Different ambiguity codes should not match
        assert not _are_nucleotides_equivalent('R', 'K', enable_iupac_intersection=False)
        
        # But standard nucleotide vs ambiguity should still work
        assert _are_nucleotides_equivalent('R', 'A', enable_iupac_intersection=False)
        assert _are_nucleotides_equivalent('A', 'R', enable_iupac_intersection=False)
    
    def test_gap_handling(self):
        """Test gap character handling."""
        assert _are_nucleotides_equivalent('-', '-')
        assert not _are_nucleotides_equivalent('-', 'A')
        assert not _are_nucleotides_equivalent('A', '-')
    
    def test_unknown_codes(self):
        """Test handling of unknown nucleotide codes."""
        # Unknown codes should only match themselves exactly
        assert _are_nucleotides_equivalent('X', 'X')
        assert not _are_nucleotides_equivalent('X', 'A')
        assert not _are_nucleotides_equivalent('A', 'X')


class TestCigarParsing:
    """Test CIGAR string parsing functions."""
    
    def test_parse_suffix_gap_deletion(self):
        """Test parsing suffix gaps (deletions in query)."""
        # Query has gaps at end (D = deletion from query)
        result = _parse_suffix_gap_from_cigar("6=3D")
        assert result == (3, True)  # 3 bp gap in query
    
    def test_parse_suffix_gap_insertion(self):
        """Test parsing suffix gaps (insertions in query)."""
        # Target has gaps at end (I = insertion in query)
        result = _parse_suffix_gap_from_cigar("6=3I")
        assert result == (3, False)  # 3 bp gap in target
    
    def test_parse_no_suffix_gap(self):
        """Test CIGAR with no suffix gap."""
        result = _parse_suffix_gap_from_cigar("6=")
        assert result is None
        
        result = _parse_suffix_gap_from_cigar("3=2D3=")  # Gap in middle
        assert result is None
    
    def test_parse_mixed_suffix_gaps(self):
        """Test CIGAR with mixed gap types at end (should reset)."""
        result = _parse_suffix_gap_from_cigar("6=2D1I")
        assert result is None  # Mixed types, not clean suffix
    
    def test_parse_empty_cigar(self):
        """Test empty or invalid CIGAR strings."""
        assert _parse_suffix_gap_from_cigar("") is None
        assert _parse_suffix_gap_from_cigar("invalid") is None


class TestScoringRegion:
    """Test scoring region identification for end trimming."""
    
    def test_no_trimming_short_sequence(self):
        """Short sequences should not be trimmed."""
        seq1 = "ATCGATCG"
        seq2 = "ATCGATCG"
        start, end = _find_scoring_region(seq1, seq2, end_skip_distance=20)
        assert start == 0
        assert end == 7  # Full sequence
    
    def test_trimming_long_sequence(self):
        """Long sequences should have ends trimmed."""
        # Create sequences with 25bp on each side
        seq1 = "A" * 25 + "ATCG" + "T" * 25
        seq2 = "A" * 25 + "ATCG" + "T" * 25
        start, end = _find_scoring_region(seq1, seq2, end_skip_distance=20)
        
        # Should skip first 20 and last 20 from each sequence (0-indexed, so >= 19)
        assert start >= 19
        assert end <= len(seq1) - 20  # end can be at the boundary
    
    def test_trimming_with_gaps(self):
        """Trimming should account for gaps in sequences."""
        # Sequence with gaps at start
        seq1 = "----" + "A" * 25 + "ATCG" + "T" * 25
        seq2 = "AAAA" + "A" * 25 + "ATCG" + "T" * 25
        start, end = _find_scoring_region(seq1, seq2, end_skip_distance=20)
        
        # Should account for gap positions
        assert start > 0
        assert end < len(seq1) - 1


class TestIUPACCodes:
    """Test IUPAC code definitions."""
    
    def test_iupac_code_definitions(self):
        """Test that IUPAC codes are defined correctly."""
        assert 'A' in IUPAC_CODES['R']  # R contains A
        assert 'G' in IUPAC_CODES['R']  # R contains G
        assert 'C' in IUPAC_CODES['Y']  # Y contains C
        assert 'T' in IUPAC_CODES['Y']  # Y contains T
        assert len(IUPAC_CODES['N']) == 4  # N contains all nucleotides
    
    def test_standard_nucleotides(self):
        """Test standard nucleotide definitions."""
        assert IUPAC_CODES['A'] == {'A'}
        assert IUPAC_CODES['T'] == {'T'}
        assert IUPAC_CODES['C'] == {'C'}
        assert IUPAC_CODES['G'] == {'G'}
    
    def test_gap_definition(self):
        """Test gap character definition."""
        assert IUPAC_CODES['-'] == {'-'}