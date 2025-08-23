#!/usr/bin/env python3
"""
Comprehensive test suite for score_alignment function.

This test suite serves as both documentation and validation for the various
adjustment parameters and scoring behaviors. Each test case demonstrates
expected behavior for specific sequence patterns and parameter combinations.
"""

import pytest
from adjusted_identity import (
    score_alignment,
    AdjustmentParams,
    ScoringFormat,
    DEFAULT_ADJUSTMENT_PARAMS,
    RAW_ADJUSTMENT_PARAMS,
)


class TestBasicMatching:
    """Test basic sequence matching without adjustments."""
    
    def test_perfect_match(self):
        """Perfect match should have 100% identity."""
        result = score_alignment("ATCG", "ATCG", RAW_ADJUSTMENT_PARAMS)
        assert result.identity == 1.0
        assert result.mismatches == 0
        assert result.scored_positions == 4
        assert result.score_aligned == "||||"
    
    def test_single_substitution(self):
        """Single nucleotide substitution."""
        result = score_alignment("ATCG", "ATCX", RAW_ADJUSTMENT_PARAMS)
        assert result.identity == 0.75  # 3/4
        assert result.mismatches == 1
        assert result.scored_positions == 4
        assert result.score_aligned == "||| "
    
    def test_multiple_substitutions(self):
        """Multiple substitutions."""
        result = score_alignment("ATCG", "XXXX", RAW_ADJUSTMENT_PARAMS)
        assert result.identity == 0.0
        assert result.mismatches == 4
        assert result.scored_positions == 4
        assert result.score_aligned == "    "


class TestIndelScoring:
    """Test indel (insertion/deletion) scoring with and without normalization."""
    
    def test_single_insertion_raw(self):
        """Single insertion without normalization - each gap position counted."""
        result = score_alignment("ATC-G", "ATCXG", RAW_ADJUSTMENT_PARAMS)
        assert result.identity == 0.8  # 4/5
        assert result.mismatches == 1
        assert result.scored_positions == 5
        assert result.score_aligned == "||| |"
    
    def test_single_insertion_normalized(self):
        """Single insertion with normalization - single position indel behavior same as raw."""
        result = score_alignment("ATC-G", "ATCXG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 0.8  # 4/5 (single position indel)
        assert result.mismatches == 1
        assert result.scored_positions == 5
        assert result.score_aligned == "||| |"
    
    def test_multi_position_indel_raw(self):
        """Multi-position indel without normalization."""
        result = score_alignment("AT---G", "ATXXXG", RAW_ADJUSTMENT_PARAMS)
        assert result.identity == 0.5  # 3/6
        assert result.mismatches == 3
        assert result.scored_positions == 6
        assert result.score_aligned == "||   |"
    
    def test_multi_position_indel_normalized(self):
        """Multi-position indel with normalization - counts as single event."""
        result = score_alignment("AT---G", "ATXXXG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 0.75  # 3/4
        assert result.mismatches == 1
        assert result.scored_positions == 4
        assert result.score_aligned == "|| --|"
    
    def test_multiple_separate_indels_normalized(self):
        """Separate indels should each count as one event when normalized."""
        result = score_alignment("A-TC-G", "AXTCXG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == pytest.approx(0.667, abs=0.001)  # 4/6
        assert result.mismatches == 2
        assert result.scored_positions == 6
        assert result.score_aligned == "| || |"


class TestHomopolymerAdjustment:
    """Test homopolymer length normalization."""
    
    def test_homopolymer_extension_adjustment_enabled(self):
        """Homopolymer extension should be ignored when adjustment enabled."""
        # Extra A extends the A homopolymer
        result = score_alignment("AAA-TTT", "AAAATTT", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 1.0  # Homopolymer extension ignored
        assert result.mismatches == 0
        assert result.scored_positions == 6  # Extension not counted in denominator
        assert result.score_aligned == "|||=|||"
    
    def test_homopolymer_extension_adjustment_disabled(self):
        """Homopolymer extension should count as mismatch when adjustment disabled."""
        result = score_alignment("AAA-TTT", "AAAATTT", RAW_ADJUSTMENT_PARAMS)
        assert result.identity == pytest.approx(6/7, abs=0.001)  # 6 matches, 1 indel
        assert result.mismatches == 1
        assert result.scored_positions == 7
        assert result.score_aligned == "||| |||"
    
    def test_homopolymer_deletion(self):
        """Homopolymer shortening should be ignored when adjustment enabled."""
        result = score_alignment("AAAATTT", "AAA-TTT", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 1.0
        assert result.mismatches == 0
        assert result.scored_positions == 6
        assert result.score_aligned == "|||=|||"
    
    def test_non_homopolymer_indel_in_homopolymer_region(self):
        """Non-homopolymer indels in homopolymer regions should still count."""
        # G insertion in A homopolymer region is not a homopolymer extension
        result = score_alignment("AAA-AAA", "AAAGAAA", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == pytest.approx(6/7, abs=0.001)  # Should count as regular indel
        assert result.mismatches == 1
        assert result.scored_positions == 7
        assert result.score_aligned == "||| |||"
    
    def test_complex_homopolymer_scenario(self):
        """Multiple homopolymer extensions in same sequence."""
        result = score_alignment("AA--TTTT--GG", "AAAATTTTTTGG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 1.0
        assert result.mismatches == 0
        assert result.scored_positions == 8  # 2A + 4T + 2G, extensions excluded
        assert result.score_aligned == "||==||||==||"


class TestIUPACAdjustment:
    """Test IUPAC ambiguity code handling."""
    
    def test_exact_iupac_match(self):
        """Same IUPAC codes should always match."""
        result = score_alignment("ANRG", "ANRG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 1.0
        assert result.mismatches == 0
        assert result.scored_positions == 4
        assert result.score_aligned == "||||"
    
    def test_iupac_intersection_match(self):
        """Different IUPAC codes with overlapping nucleotides should match."""
        # R (AG) and K (GT) both contain G
        result = score_alignment("ATRG", "ATKG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 1.0
        assert result.mismatches == 0
        assert result.scored_positions == 4
        assert result.score_aligned == "||||"
    
    def test_iupac_intersection_disabled(self):
        """IUPAC intersection disabled should not match different ambiguity codes."""
        params = AdjustmentParams(handle_iupac_overlap=False)
        result = score_alignment("ATRG", "ATKG", params)
        assert result.identity == 0.75  # 3/4
        assert result.mismatches == 1
        assert result.scored_positions == 4
        assert result.score_aligned == "|| |"
    
    def test_iupac_no_intersection(self):
        """IUPAC codes with no overlap should not match even with adjustment enabled."""
        # R (AG) and Y (CT) have no overlap
        result = score_alignment("ATRG", "ATYG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 0.75
        assert result.mismatches == 1
        assert result.scored_positions == 4
        assert result.score_aligned == "|| |"
    
    def test_standard_vs_iupac(self):
        """Standard nucleotide vs IUPAC should match if overlap exists."""
        # R (AG) contains A
        result = score_alignment("ATRG", "ATAG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 1.0
        assert result.mismatches == 0
        assert result.scored_positions == 4
        assert result.score_aligned == "||||"


class TestEndTrimming:
    """Test end region mismatch skipping."""
    
    def test_end_trimming_enabled_default(self):
        """Mismatches near ends should be skipped with default 20bp trimming."""
        # Create sequences longer than 40bp to enable trimming
        seq1 = "A" * 21 + "XXXX" + "T" * 21  # Mismatches in middle
        seq2 = "A" * 21 + "TTTT" + "T" * 21
        result = score_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)
        
        # After trimming 20bp from each end, 8 positions remain in scoring region
        assert result.mismatches == 4
        assert result.scored_positions == 8
    
    def test_end_trimming_disabled(self):
        """All positions should be scored when end trimming disabled."""
        params = AdjustmentParams(end_skip_distance=0)
        seq1 = "XXXXTTTTXXXX"
        seq2 = "AAAAAAAAAAAA"
        result = score_alignment(seq1, seq2, params)
        
        assert result.mismatches == 12
        assert result.scored_positions == 12
        assert result.score_aligned == "            "
    
    def test_short_sequence_no_trimming(self):
        """Short sequences should not be affected by end trimming."""
        # Sequence too short for 20bp trimming on each end
        result = score_alignment("ATCGXXXX", "ATCGTTTT", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.mismatches == 4
        assert result.scored_positions == 8


class TestCombinedAdjustments:
    """Test combinations of different adjustments."""
    
    def test_all_adjustments_enabled(self):
        """Complex case with all adjustments enabled."""
        # Homopolymer extension + IUPAC + indel normalization
        result = score_alignment("AAA-TTRG", "AAAATTKG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 1.0  # All differences adjusted away
        assert result.mismatches == 0
        assert result.scored_positions == 7
        assert result.score_aligned == "|||=||||"
    
    def test_all_adjustments_disabled(self):
        """Same case with all adjustments disabled."""
        result = score_alignment("AAA-TTRG", "AAAATTKG", RAW_ADJUSTMENT_PARAMS)
        assert result.identity < 1.0  # Raw scoring shows differences
        assert result.mismatches > 0
    
    def test_selective_adjustments(self):
        """Enable only some adjustments."""
        params = AdjustmentParams(
            normalize_homopolymers=True,
            handle_iupac_overlap=False,  # Disabled
            normalize_indels=True,
            end_skip_distance=0
        )
        result = score_alignment("AAA-TTRG", "AAAATTKG", params)
        # Homopolymer and indel adjustments work, but IUPAC mismatch counts
        assert result.identity == pytest.approx(6/7, abs=0.001)  # One IUPAC mismatch remains
        assert result.mismatches == 1


class TestEdgeCases:
    """Test edge cases and error conditions."""
    
    def test_empty_sequences(self):
        """Empty sequences should be handled gracefully."""
        result = score_alignment("", "", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 1.0
        assert result.mismatches == 0
        assert result.scored_positions == 0
    
    def test_one_empty_sequence(self):
        """One empty sequence should result in all gaps."""
        result = score_alignment("----", "ATCG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 0.0
        assert result.mismatches == 1  # Normalized as single indel
        assert result.scored_positions == 1
    
    def test_all_gaps(self):
        """All-gap alignment should be handled."""
        result = score_alignment("----", "----", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 1.0  # No mismatches in all-gap alignment
        assert result.mismatches == 0
        assert result.scored_positions == 4  # Gap characters are still scored as matches
    
    def test_unequal_length_sequences(self):
        """Aligned sequences of different lengths should raise error."""
        with pytest.raises(ValueError, match="same length"):
            score_alignment("ATCG", "ATCGX", DEFAULT_ADJUSTMENT_PARAMS)


class TestScoringFormatCustomization:
    """Test custom scoring format codes."""
    
    def test_custom_scoring_format(self):
        """Custom scoring format should be used in output."""
        custom_format = ScoringFormat(
            match='*',
            substitution='X',
            indel_start='I',
            homopolymer_extension='H'
        )
        result = score_alignment("AAA-TTT", "AAAATTT", DEFAULT_ADJUSTMENT_PARAMS, custom_format)
        assert result.score_aligned == "***H***"
    
    def test_invalid_scoring_format(self):
        """Invalid scoring format should raise error."""
        with pytest.raises(ValueError):
            ScoringFormat(match="too_long")


class TestDocumentationExamples:
    """Test cases that serve as clear documentation examples."""
    
    def test_mycology_example_homopolymer(self):
        """Example: Sequencing artifact in fungal ITS region."""
        # Common scenario: different homopolymer lengths in ITS sequences
        its_seq1 = "ATCGAAAAATGTC"  # 5 A's
        its_seq2 = "ATCGAA-AATGTC"  # 4 A's (gap represents shorter homopolymer)
        
        # With adjustment: sequences are considered identical
        adjusted = score_alignment(its_seq1, its_seq2, DEFAULT_ADJUSTMENT_PARAMS)
        assert adjusted.identity == 1.0
        
        # Without adjustment: homopolymer difference counts as mismatch
        raw = score_alignment(its_seq1, its_seq2, RAW_ADJUSTMENT_PARAMS)
        assert raw.identity < 1.0
    
    def test_ambiguous_barcoding_example(self):
        """Example: IUPAC codes in DNA barcoding."""
        # Sequencing may produce ambiguous calls at same position
        barcode1 = "ATCGRGTC"  # R = A or G
        barcode2 = "ATCGKGTC"  # K = G or T
        
        # Both R and K contain G, so with IUPAC adjustment they match
        result = score_alignment(barcode1, barcode2, DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 1.0
    
    def test_end_trimming_example(self):
        """Example: Poor quality sequence ends."""
        # Poor quality at sequence ends (common with Sanger sequencing)
        seq1 = "N" * 25 + "ATCGATCGATCG" + "N" * 25  # Good sequence in middle
        seq2 = "X" * 25 + "ATCGATCGATCG" + "Y" * 25  # Same middle, bad ends
        
        # With end trimming: some end mismatches ignored but not all
        result = score_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 0.75  # Some end mismatches still counted
        
        # Without end trimming: end mismatches count
        no_trim = AdjustmentParams(end_skip_distance=0)
        result_no_trim = score_alignment(seq1, seq2, no_trim)
        assert result_no_trim.identity < result.identity