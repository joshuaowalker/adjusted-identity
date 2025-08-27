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
    _reverse_complement,
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
        # Standard nucleotides should match exactly and not be ambiguous
        assert _are_nucleotides_equivalent('A', 'A') == (True, False)
        assert _are_nucleotides_equivalent('T', 'T') == (True, False)
        assert _are_nucleotides_equivalent('C', 'C') == (True, False)
        assert _are_nucleotides_equivalent('G', 'G') == (True, False)
        # N=N is ambiguous since N is an ambiguity code
        assert _are_nucleotides_equivalent('N', 'N') == (True, True)
    
    def test_case_insensitive(self):
        """Test case insensitive matching."""
        assert _are_nucleotides_equivalent('a', 'A') == (True, False)
        assert _are_nucleotides_equivalent('A', 'a') == (True, False)
        assert _are_nucleotides_equivalent('r', 'R') == (True, True)  # R=R is ambiguous
    
    def test_iupac_intersections_enabled(self):
        """Test IUPAC code intersections when enabled."""
        # R (AG) and K (GT) both contain G - ambiguous match
        assert _are_nucleotides_equivalent('R', 'K', enable_iupac_intersection=True) == (True, True)
        
        # R (AG) contains A - ambiguous match
        assert _are_nucleotides_equivalent('R', 'A', enable_iupac_intersection=True) == (True, True)
        
        # Y (CT) and S (GC) both contain C - ambiguous match
        assert _are_nucleotides_equivalent('Y', 'S', enable_iupac_intersection=True) == (True, True)
        
        # R (AG) and Y (CT) have no overlap
        assert _are_nucleotides_equivalent('R', 'Y', enable_iupac_intersection=True) == (False, False)
    
    def test_iupac_intersections_disabled(self):
        """Test IUPAC code handling when intersections disabled."""
        # Different ambiguity codes should not match
        assert _are_nucleotides_equivalent('R', 'K', enable_iupac_intersection=False) == (False, False)
        
        # But standard nucleotide vs ambiguity should still work - ambiguous match
        assert _are_nucleotides_equivalent('R', 'A', enable_iupac_intersection=False) == (True, True)
        assert _are_nucleotides_equivalent('A', 'R', enable_iupac_intersection=False) == (True, True)
    
    def test_gap_handling(self):
        """Test gap character handling."""
        # Gap matching gap is not ambiguous (it's a special case)
        assert _are_nucleotides_equivalent('-', '-') == (True, True)  # Gap is not a standard nucleotide
        assert _are_nucleotides_equivalent('-', 'A') == (False, False)
        assert _are_nucleotides_equivalent('A', '-') == (False, False)
    
    def test_unknown_codes(self):
        """Test handling of unknown nucleotide codes."""
        # Unknown codes should only match themselves exactly - ambiguous since not standard
        assert _are_nucleotides_equivalent('X', 'X') == (True, True)
        assert _are_nucleotides_equivalent('X', 'A') == (False, False)
        assert _are_nucleotides_equivalent('A', 'X') == (False, False)


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


class TestReverseComplement:
    """Test custom reverse complement implementation."""
    
    def test_standard_nucleotides(self):
        """Test reverse complement of standard nucleotides."""
        assert _reverse_complement('A') == 'T'
        assert _reverse_complement('T') == 'A'
        assert _reverse_complement('G') == 'C' 
        assert _reverse_complement('C') == 'G'
        
        # Case variants
        assert _reverse_complement('a') == 't'
        assert _reverse_complement('t') == 'a'
        assert _reverse_complement('g') == 'c'
        assert _reverse_complement('c') == 'g'
    
    def test_standard_sequences(self):
        """Test reverse complement of standard DNA sequences."""
        assert _reverse_complement('ATCG') == 'CGAT'
        assert _reverse_complement('AAATTT') == 'AAATTT'  # Palindrome
        assert _reverse_complement('GCATGC') == 'GCATGC'  # Palindrome
        assert _reverse_complement('ATCGATCG') == 'CGATCGAT'
    
    def test_iupac_ambiguity_codes(self):
        """Test reverse complement of IUPAC ambiguity codes."""
        # Purine/Pyrimidine exchange
        assert _reverse_complement('R') == 'Y'  # AG -> CT
        assert _reverse_complement('Y') == 'R'  # CT -> AG
        
        # Keto/Amino exchange  
        assert _reverse_complement('K') == 'M'  # GT -> AC
        assert _reverse_complement('M') == 'K'  # AC -> GT
        
        # Complement sets exchange
        assert _reverse_complement('B') == 'V'  # CGT -> ACG
        assert _reverse_complement('V') == 'B'  # ACG -> CGT
        assert _reverse_complement('D') == 'H'  # AGT -> ACT
        assert _reverse_complement('H') == 'D'  # ACT -> AGT
        
        # Palindromic codes (self-complementary)
        assert _reverse_complement('S') == 'S'  # GC -> GC
        assert _reverse_complement('W') == 'W'  # AT -> AT
        assert _reverse_complement('N') == 'N'  # Any -> Any
    
    def test_case_variants(self):
        """Test case variants of IUPAC codes."""
        assert _reverse_complement('r') == 'y'
        assert _reverse_complement('y') == 'r'
        assert _reverse_complement('k') == 'm'
        assert _reverse_complement('m') == 'k'
        assert _reverse_complement('b') == 'v'
        assert _reverse_complement('v') == 'b'
        assert _reverse_complement('d') == 'h'
        assert _reverse_complement('h') == 'd'
        assert _reverse_complement('s') == 's'
        assert _reverse_complement('w') == 'w'
        assert _reverse_complement('n') == 'n'
    
    def test_mixed_case_sequences(self):
        """Test sequences with mixed upper/lower case."""
        assert _reverse_complement('AtCg') == 'cGaT'
        assert _reverse_complement('ATCGRGTC') == 'GACYCGAT'
        assert _reverse_complement('atcgrgtc') == 'gacycgat'
    
    def test_gaps_and_unknowns(self):
        """Test gap characters and unknown characters."""
        assert _reverse_complement('ATC-G') == 'C-GAT'
        assert _reverse_complement('---') == '---'
        
        # Unknown characters pass through unchanged
        assert _reverse_complement('ATCXG') == 'CXGAT'
        assert _reverse_complement('AT@CG') == 'CG@AT'
    
    def test_roundtrip_property(self):
        """Test that reverse complement of reverse complement gives original."""
        test_sequences = [
            'ATCG',
            'ATCGATCGATCG',
            'AAATTTGGGCCC',
            'ATCGRGTC',    # With IUPAC codes
            'RYSWKMBVDHN', # All IUPAC codes
            'ryswkmbvdhn', # Lowercase
            'ATC-G-TC',    # With gaps
            'ATCXRYTC',    # With unknown character
        ]
        
        for seq in test_sequences:
            roundtrip = _reverse_complement(_reverse_complement(seq))
            assert roundtrip == seq, f"Roundtrip failed for {seq}: got {roundtrip}"
    
    def test_empty_sequence(self):
        """Test empty sequence."""
        assert _reverse_complement('') == ''
    
    def test_single_characters(self):
        """Test all single characters with known mappings."""
        # Test known mappings directly
        mappings = {
            # Standard nucleotides
            'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
            'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
            # IUPAC codes
            'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
            'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
            'D': 'H', 'H': 'D', 'N': 'N', '-': '-',
            'r': 'y', 'y': 'r', 's': 's', 'w': 'w',
            'k': 'm', 'm': 'k', 'b': 'v', 'v': 'b', 
            'd': 'h', 'h': 'd', 'n': 'n'
        }
        
        for char, expected_complement in mappings.items():
            assert _reverse_complement(char) == expected_complement
            # Test roundtrip for all mapped characters
            assert _reverse_complement(_reverse_complement(char)) == char
    
    def test_palindromic_sequences(self):
        """Test sequences that are their own reverse complement."""
        palindromes = [
            'ATAT',
            'GCGC', 
            'AAATTT',
            'CCCGGG',
            'SWWS',   # IUPAC palindrome
            'NWWN',   # Mixed IUPAC palindrome
        ]
        
        for seq in palindromes:
            assert _reverse_complement(seq) == seq, f"Palindrome test failed for {seq}"
    
    def test_biological_realism(self):
        """Test with realistic biological sequences."""
        # ITS-like sequence with common IUPAC codes
        its_seq = "TCCGTAGGTGAACCTGCGGAAGGATCATTACCGAGTTTAAR"
        rc_its = _reverse_complement(its_seq)
        
        # Should be same length
        assert len(rc_its) == len(its_seq)
        
        # Should roundtrip correctly
        assert _reverse_complement(rc_its) == its_seq
        
        # Manual check of a few positions
        assert rc_its[0] == 'Y'  # Last R -> Y (at start of RC)
        assert rc_its[-1] == 'A'  # First T -> A (at end of RC)


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