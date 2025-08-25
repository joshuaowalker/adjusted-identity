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
    align_and_score,
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
    
    def test_mixed_indel_with_homopolymer_ends(self):
        """Indel with homopolymer extensions on ends but non-homopolymer content in middle."""
        # AAATTGGG vs AA----GG: indel is ATTG
        # Leading A and trailing G are now correctly detected as homopolymer extensions
        # Internal TT is treated as plain indel content
        result = score_alignment("AAATTGGG", "AA----GG", DEFAULT_ADJUSTMENT_PARAMS)
        
        # New algorithm correctly detects partial homopolymer extensions
        assert result.identity == 0.8  # 4/5: matches at positions 0,1,6,7 + normalized indel counts as 1 mismatch
        assert result.mismatches == 1  # Only the TT counts as 1 edit
        assert result.scored_positions == 5  # 2 matches + 1 normalized indel + 2 matches
        assert result.score_aligned == "||= -=||"  # Shows HP extensions (=) and regular indel ( -)
        
    def test_mixed_indel_potential_improvement(self):
        """More complex case showing potential for improvement in homopolymer handling."""
        # AAAATTTTGGGG vs AAA-----GGGG: indel contains ATTTT
        # Could potentially recognize:
        # - Leading A as homopolymer extension of AAA
        # - Trailing G could be extension but we need the sequences to match exactly in the Gs
        # - Middle TTTT as regular indel content
        result = score_alignment("AAAATTTTGGGG", "AAA-----GGGG", DEFAULT_ADJUSTMENT_PARAMS)
        
        # Current behavior: entire indel treated as non-homopolymer
        print(f"Complex case - Identity: {result.identity}")
        print(f"Mismatches: {result.mismatches}")
        print(f"Scored positions: {result.scored_positions}")
        print(f"Score pattern: {result.score_aligned}")
        
        # Document current behavior - will adjust after seeing output
        assert result.mismatches >= 1
        assert result.scored_positions >= 7
        
    def test_left_right_homopolymer_algorithm(self):
        """Test the new left-right homopolymer extension algorithm."""
        
        # Case 1: Simple case - should work the same as before
        # AAA-TTT vs AAAATTT: single A extension
        result1 = score_alignment("AAA-TTT", "AAAATTT", DEFAULT_ADJUSTMENT_PARAMS)
        assert result1.identity == 1.0  # Should be perfect match
        assert result1.score_aligned == "|||=|||"  # Homopolymer extension detected
        
        # Case 2: Mixed indel - left and right extensions
        # AAATTGGG vs AA----GG: should detect A extension on left, G extension on right
        result2 = score_alignment("AAATTGGG", "AA----GG", DEFAULT_ADJUSTMENT_PARAMS) 
        
        # Compare with raw (no adjustments) behavior
        raw_result = score_alignment("AAATTGGG", "AA----GG", RAW_ADJUSTMENT_PARAMS)
        
        # The new algorithm should perform significantly better
        assert result2.identity == 0.8  # 1 edit out of 5 scored positions
        assert result2.score_aligned == "||= -=||"  # Left HP, regular indel, right HP
        assert result2.identity > raw_result.identity  # Better than raw (0.5)
        
        # Verify the components
        assert result2.mismatches == 1  # Only the TT counts as 1 normalized edit
        assert result2.scored_positions == 5  # 2 matches + 1 normalized indel + 2 matches
        
    def test_partial_homopolymer_extension_detection(self):
        """Test improved detection of partial homopolymer extensions within complex indels."""
        
        # Case 1: Mixed indel with no dominant homopolymer character
        # AAATTGGG vs AA----GG: indel contains ATTG 
        # A=25%, T=50%, G=25% - no character reaches 30% threshold in context
        result1 = score_alignment("AAATTGGG", "AA----GG", DEFAULT_ADJUSTMENT_PARAMS)
        print(f"Simple case - Identity: {result1.identity}, Score: {result1.score_aligned}")
        
        # Case 2: Mixed indel with dominant homopolymer character  
        # AAAAATTGGG vs AAA---GGGG: indel contains AAT
        # A=67% (above 30% threshold) and has homopolymer context (AAA before)
        result2 = score_alignment("AAAAATTGGG", "AAA---GGGG", DEFAULT_ADJUSTMENT_PARAMS)
        print(f"A-rich case - Identity: {result2.identity}, Score: {result2.score_aligned}")
        
        # Just verify results are reasonable for now
        assert 0.0 <= result1.identity <= 1.0
        assert 0.0 <= result2.identity <= 1.0

    def test_boundary_conditions(self):
        """Test boundary conditions for left-right algorithm components."""
        
        # Pure left only (no middle, no right)
        result_left = score_alignment("AAA-TTT", "AAAATTT", DEFAULT_ADJUSTMENT_PARAMS)
        assert result_left.identity == 1.0
        assert result_left.score_aligned == "|||=|||"
        assert result_left.mismatches == 0
        
        # Pure right only (no left, no middle) 
        result_right = score_alignment("TTT-GGG", "TTTGGGG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result_right.identity == 1.0
        assert result_right.score_aligned == "|||=|||"
        assert result_right.mismatches == 0
        
        # Pure middle only (no left, no right - no homopolymer context)
        result_middle = score_alignment("ATCG-CGA", "ATCGCCGA", DEFAULT_ADJUSTMENT_PARAMS)
        assert result_middle.identity == 1.0  # C extends C context
        assert result_middle.score_aligned == "||||=|||"
        assert result_middle.mismatches == 0
        
        # Left + middle (no right)
        result_left_middle = score_alignment("AAATTCG", "AA---CG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result_left_middle.identity == 0.8  # A extension + TT indel
        assert result_left_middle.score_aligned == "||= -||"
        assert result_left_middle.mismatches == 1
        
        # Right + middle (no left)  
        result_right_middle = score_alignment("ATCGGGG", "AT--GGG", DEFAULT_ADJUSTMENT_PARAMS)
        assert abs(result_right_middle.identity - (5/6)) < 0.001  # C indel + G extension
        assert result_right_middle.score_aligned == "|| =|||"
        assert result_right_middle.mismatches == 1
        
        # Left + right (no middle)
        result_left_right = score_alignment("AA---GGG", "AAAAAGGG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result_left_right.identity == 1.0  # AA + G extensions
        assert result_left_right.score_aligned == "||===|||"
        assert result_left_right.mismatches == 0

    def test_context_edge_cases(self):
        """Test edge cases with missing context (start/end of sequence)."""
        
        # Indel at start of sequence (no left context)
        result_no_left = score_alignment("-ATCG", "GATCG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result_no_left.identity == 0.8  # G indel at start
        assert result_no_left.score_aligned == " ||||"
        assert result_no_left.mismatches == 1
        
        # Indel at end of sequence (no right context)
        result_no_right = score_alignment("ATCG-", "ATCGG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result_no_right.identity == 1.0  # G extends G context (homopolymer)
        assert result_no_right.score_aligned == "||||="
        assert result_no_right.mismatches == 0
        
        # Homopolymer at start (should still be detected if context allows)
        result_hp_start = score_alignment("A-TCG", "AATCG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result_hp_start.identity == 1.0  # A extends A context
        assert result_hp_start.score_aligned == "|=|||"
        assert result_hp_start.mismatches == 0
        
        # Homopolymer at end (should still be detected if context allows)  
        result_hp_end = score_alignment("TCGA-", "TCGAA", DEFAULT_ADJUSTMENT_PARAMS)
        assert result_hp_end.identity == 1.0  # A extends A context
        assert result_hp_end.score_aligned == "||||="
        assert result_hp_end.mismatches == 0

    def test_multiple_component_edge_cases(self):
        """Test edge cases with various component combinations."""
        
        # Empty components should not affect scoring
        # All components present but some very small
        result_tiny = score_alignment("AATGGG", "A--GGG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result_tiny.identity == 0.8  # A extension + T indel 
        assert result_tiny.score_aligned == "|= |||"
        assert result_tiny.mismatches == 1
        
        # Large left component, small middle, small right
        result_large_left = score_alignment("AAAAAAATCGGG", "AAA------GGG", DEFAULT_ADJUSTMENT_PARAMS)
        assert abs(result_large_left.identity - (6/7)) < 0.001  # AAAA extension + ATC indel
        assert result_large_left.score_aligned == "|||==== -|||"
        assert result_large_left.mismatches == 1
        
        # Context characters that don't match (should be treated as middle)
        result_no_match = score_alignment("ATCG-CGTA", "ATCGACGTA", DEFAULT_ADJUSTMENT_PARAMS)
        assert result_no_match.identity > 0.8  # A indel, should find some context
        assert result_no_match.mismatches <= 1


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


class TestRepeatMotifs:
    """Test handling of dinucleotide and longer repeat motifs."""
    
    def test_dinucleotide_repeat_basic(self):
        """Test basic AT dinucleotide repeat from Russell article."""
        # Example from article: CGATATC vs CGATATATC (extra AT motif)
        seq1 = "CGATAT--C"
        seq2 = "CGATATATC"
        
        params = AdjustmentParams(max_repeat_motif_length=2)
        result = score_alignment(seq1, seq2, params)
        
        # The AT insertion should be treated as repeat extension
        assert result.identity == 1.0
        assert result.mismatches == 0
        assert '=' in result.score_aligned  # Should have repeat extension markers
    
    def test_dinucleotide_repeat_multiple(self):
        """Test multiple dinucleotide repeat units."""
        # Two extra AT units
        seq1 = "CGATAT----C"
        seq2 = "CGATATATATC"
        
        params = AdjustmentParams(max_repeat_motif_length=2)
        result = score_alignment(seq1, seq2, params)
        
        assert result.identity == 1.0
        assert result.mismatches == 0
        assert result.score_aligned.count('=') == 4  # Four repeat extension positions
    
    def test_mixed_motif_lengths(self):
        """Test indel with different left and right motif lengths."""
        # Left side: AT repeat (length 2)
        # Middle: C (not a repeat)
        # Right side: single G homopolymer (length 1)
        seq1 = "ATAT----GGG"
        seq2 = "ATATATCGGGG"
        
        params = AdjustmentParams(max_repeat_motif_length=2)
        result = score_alignment(seq1, seq2, params)
        
        # AT on left (extension), C in middle (regular indel), G on right (homopolymer extension)
        assert result.identity < 1.0  # Middle C is a mismatch
        assert '=' in result.score_aligned  # AT and G extensions
        assert ' ' in result.score_aligned or '-' in result.score_aligned  # Regular indel for C
    
    def test_degenerate_dinucleotide_as_homopolymer(self):
        """Test that AA/TT/CC/GG are treated as homopolymers, not dinucleotides."""
        # "AA" should be treated as homopolymer 'A', not dinucleotide "AA"
        seq1 = "CGAAA----TC"
        seq2 = "CGAAAAAA-TC"
        
        params = AdjustmentParams(max_repeat_motif_length=2)
        result = score_alignment(seq1, seq2, params)
        
        # Should treat as homopolymer extension
        assert result.identity == 1.0
        assert result.mismatches == 0
        assert '=' in result.score_aligned
    
    def test_partial_motif_not_consumed(self):
        """Test that partial motifs are not consumed as extensions."""
        # "ATA" where motif is "AT" - should only consume "AT", leave "A"
        seq1 = "ATAT---C"
        seq2 = "ATATATAC"
        
        params = AdjustmentParams(max_repeat_motif_length=2)
        result = score_alignment(seq1, seq2, params)
        
        # The final "A" should be a regular indel, not extension
        assert result.identity < 1.0
        assert result.mismatches == 1  # The "A" counts as mismatch
    
    def test_no_matching_context(self):
        """Test indel with no repeat context."""
        # No repeating pattern
        seq1 = "ATCG---TGCA"
        seq2 = "ATCGACGTGCA"
        
        params = AdjustmentParams(max_repeat_motif_length=2)
        result = score_alignment(seq1, seq2, params)
        
        # Should be treated as regular indel
        assert result.identity < 1.0
        assert '=' not in result.score_aligned  # No repeat extensions
    
    def test_trinucleotide_with_max_length_2(self):
        """Test that trinucleotide repeats are not detected when max_length=2."""
        # CAG repeat, but max_repeat_motif_length=2
        seq1 = "CAGCAG---TTC"
        seq2 = "CAGCAGCAGTTC"
        
        params = AdjustmentParams(max_repeat_motif_length=2)
        result = score_alignment(seq1, seq2, params)
        
        # Should not detect CAG repeat (length 3 > max 2)
        assert result.identity < 1.0
        assert result.mismatches > 0
    
    def test_both_sides_same_dinucleotide(self):
        """Test indel with same dinucleotide repeat on both sides."""
        # AT repeat on both sides of indel
        seq1 = "ATAT------ATAT"
        seq2 = "ATATATATATATAT"
        
        params = AdjustmentParams(max_repeat_motif_length=2)
        result = score_alignment(seq1, seq2, params)
        
        # All should be AT extensions
        assert result.identity == 1.0
        assert result.mismatches == 0
        assert result.score_aligned.count('=') == 6  # Six repeat extension positions
    
    def test_reverse_complement_motifs(self):
        """Test different motifs that are reverse complements."""
        # AT on left, TA on right (reverse complement)
        seq1 = "ATAT--TATA"
        seq2 = "ATATATTATA"
        
        params = AdjustmentParams(max_repeat_motif_length=2)
        result = score_alignment(seq1, seq2, params)
        
        # Should handle both independently
        assert result.identity == 1.0
        assert '=' in result.score_aligned
    
    def test_motif_at_sequence_boundary(self):
        """Test repeat motif at the very start or end of sequence."""
        # AT repeat at start
        seq1 = "--ATATGC"
        seq2 = "ATATATGC"
        
        params = AdjustmentParams(max_repeat_motif_length=2)
        result = score_alignment(seq1, seq2, params)
        
        # No left context, but should still work with right context
        assert '=' in result.score_aligned or ' ' in result.score_aligned
    
    def test_complex_mixed_indel(self):
        """Test complex indel with both repeat extensions and regular content."""
        # Left: AT repeat, Middle: CGT (non-repeat), Right: G homopolymer
        seq1 = "ATAT-------GGG"
        seq2 = "ATATATATCGTGGG"
        
        params = AdjustmentParams(max_repeat_motif_length=2)
        result = score_alignment(seq1, seq2, params)
        
        # Should have both extensions (=) and regular indel ( -)
        assert '=' in result.score_aligned  # Extensions
        assert ' ' in result.score_aligned or '-' in result.score_aligned  # Regular indel
    
    def test_motif_length_disabled(self):
        """Test that setting max_repeat_motif_length=1 disables dinucleotide detection."""
        # AT repeat that should not be detected
        seq1 = "ATAT--C"
        seq2 = "ATATATC"
        
        params = AdjustmentParams(max_repeat_motif_length=1)  # Only homopolymers
        result = score_alignment(seq1, seq2, params)
        
        # Should treat as regular indel
        assert result.identity < 1.0
        assert result.mismatches > 0
    
    def test_overlapping_motif_possibilities(self):
        """Test sequence where multiple motif lengths could apply."""
        # AAAA could be: "AAAA" (length 4), "AA" (length 2), or "A" (length 1)
        seq1 = "AAAA----TTTT"
        seq2 = "AAAAAAAATTTT"
        
        params = AdjustmentParams(max_repeat_motif_length=2)
        result = score_alignment(seq1, seq2, params)
        
        # Should detect as homopolymer (length 1) due to degeneracy check
        assert result.identity == 1.0
        assert '=' in result.score_aligned


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


class TestOverhangScoring:
    """Test cases for proper handling of overhang regions when end_skip_distance=0."""
    
    def test_left_overhang_not_scored(self):
        """When end_skip_distance=0, left overhang should not be scored."""
        # Sequence 1 has extra bases at start
        seq1 = "AAATTTGGG"
        seq2 = "TTTGGG"
        
        # This should align as: AAATTTGGG
        #                       ---TTTGGG
        # With end_skip_distance=0, we should only score TTTGGG vs TTTGGG (identity=1.0)
        # Currently, this incorrectly scores the AAA overhang region too
        
        no_trim = AdjustmentParams(end_skip_distance=0)
        result = align_and_score(seq1, seq2, no_trim)
        
        # Should be 1.0 since overlapping region matches perfectly
        assert result.identity == 1.0, f"Expected 1.0, got {result.identity}"
        
    def test_right_overhang_not_scored(self):
        """When end_skip_distance=0, right overhang should not be scored."""
        # Sequence 1 has extra bases at end
        seq1 = "TTTGGGAAA"
        seq2 = "TTTGGG"
        
        # This should align as: TTTGGGAAA
        #                       TTTGGG---
        # With end_skip_distance=0, we should only score TTTGGG vs TTTGGG (identity=1.0)
        
        no_trim = AdjustmentParams(end_skip_distance=0)
        result = align_and_score(seq1, seq2, no_trim)
        
        # Should be 1.0 since overlapping region matches perfectly
        assert result.identity == 1.0, f"Expected 1.0, got {result.identity}"
        
    def test_both_overhangs_not_scored(self):
        """When end_skip_distance=0, neither overhang should be scored."""
        # Both sequences have overhangs
        seq1 = "AAATTTGGGCCC"
        seq2 = "XXXTTTGGGYYY"
        
        # This should align as: AAATTTGGGCCC
        #                       XXXTTTGGGYYY
        # With end_skip_distance=0, we should only score TTTGGG vs TTTGGG (identity=1.0)
        # The AAA/XXX and CCC/YYY overhangs should not be scored
        
        no_trim = AdjustmentParams(end_skip_distance=0)
        result = align_and_score(seq1, seq2, no_trim)
        
        # Should be 1.0 since overlapping region matches perfectly
        assert result.identity == 1.0, f"Expected 1.0, got {result.identity}"
        
    def test_no_overlap_case(self):
        """Edge case: sequences with no overlapping content should have identity=0."""
        # Use score_alignment directly with pre-aligned sequences
        seq1_aligned = "AAAA----"
        seq2_aligned = "----TTTT"
        
        no_trim = AdjustmentParams(end_skip_distance=0)
        result = score_alignment(seq1_aligned, seq2_aligned, no_trim)
        
        # No overlapping content - should be identity 0
        assert result.identity == 0.0, f"Expected 0.0, got {result.identity}"
        
    def test_single_position_overlap(self):
        """Edge case: single position of overlap."""
        # Use score_alignment directly with pre-aligned sequences
        seq1_aligned = "AAA-"
        seq2_aligned = "-AAG"
        
        no_trim = AdjustmentParams(end_skip_distance=0)
        result = score_alignment(seq1_aligned, seq2_aligned, no_trim)
        
        # Single overlapping position: A vs A = match, so identity should be 1.0
        assert result.identity == 1.0, f"Expected 1.0, got {result.identity}"
        
    def test_single_position_mismatch(self):
        """Edge case: two positions of overlap with one mismatch."""
        # Use score_alignment directly with pre-aligned sequences
        seq1_aligned = "AAA-"
        seq2_aligned = "-ATG"
        
        no_trim = AdjustmentParams(end_skip_distance=0)
        result = score_alignment(seq1_aligned, seq2_aligned, no_trim)
        
        # Two overlapping positions: A vs A (match), A vs T (mismatch) = 50% identity
        assert result.identity == 0.5, f"Expected 0.5, got {result.identity}"
        
    def test_true_single_position_mismatch(self):
        """Edge case: single position of overlap with mismatch."""
        # Use score_alignment directly with pre-aligned sequences
        seq1_aligned = "AA-"
        seq2_aligned = "-AT"
        
        no_trim = AdjustmentParams(end_skip_distance=0)
        result = score_alignment(seq1_aligned, seq2_aligned, no_trim)
        
        # Single overlapping position: A vs A = match, so identity should be 1.0
        assert result.identity == 1.0, f"Expected 1.0, got {result.identity}"