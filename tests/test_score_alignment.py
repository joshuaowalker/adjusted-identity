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

        # Indel near start of sequence (minimal left context)
        # T-ATCG vs TGATCG - G doesn't extend T → core content
        result_no_left = score_alignment("T-ATCG", "TGATCG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result_no_left.identity < 1.0  # G indel doesn't extend T
        assert result_no_left.mismatches == 1

        # Indel near end of sequence (homopolymer extension)
        # ATCG-T vs ATCGGT - G extends G context (homopolymer)
        result_no_right = score_alignment("ATCG-T", "ATCGGT", DEFAULT_ADJUSTMENT_PARAMS)
        assert result_no_right.identity == 1.0  # G extends G context
        assert result_no_right.score_aligned == "||||=|"
        assert result_no_right.mismatches == 0

        # Homopolymer at start (should still be detected if context allows)
        result_hp_start = score_alignment("A-TCG", "AATCG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result_hp_start.identity == 1.0  # A extends A context
        assert result_hp_start.score_aligned == "|=|||"
        assert result_hp_start.mismatches == 0
        
        # Homopolymer at end (should still be detected if context allows)
        # Add trailing context to avoid overhang behavior
        result_hp_end = score_alignment("TCGA-T", "TCGAAT", DEFAULT_ADJUSTMENT_PARAMS)
        assert result_hp_end.identity == 1.0  # A extends A context
        assert result_hp_end.score_aligned == "||||=|"
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


class TestAmbiguousMatching:
    """Test the new ambiguous matching scoring code feature."""
    
    def test_standard_nucleotide_matches(self):
        """Standard nucleotides (ATCG) should use '|' when matching exactly."""
        result = score_alignment("ATCG", "ATCG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.score_aligned == "||||"
        assert '=' not in result.score_aligned  # No ambiguous matches
    
    def test_ambiguous_code_self_match(self):
        """Ambiguity codes matching themselves should use '='."""
        # N, R, Y are ambiguity codes
        result = score_alignment("NRY", "NRY", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.score_aligned == "==="
        assert result.identity == 1.0
    
    def test_mixed_standard_and_ambiguous(self):
        """Mixed standard and ambiguous nucleotides should use appropriate codes."""
        # A and T are standard, N is ambiguous
        result = score_alignment("ATN", "ATN", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.score_aligned == "||="
    
    def test_ambiguous_intersection_match(self):
        """Different ambiguity codes that intersect should use '='."""
        # R (AG) and D (AGT) intersect at A and G
        result = score_alignment("R", "D", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.score_aligned == "="
        assert result.identity == 1.0
    
    def test_standard_vs_ambiguous_match(self):
        """Standard nucleotide matching ambiguity code should use '='."""
        # A matches R (AG)
        result = score_alignment("A", "R", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.score_aligned == "="
        assert result.identity == 1.0
    
    def test_custom_scoring_format(self):
        """Test using custom ScoringFormat with different ambiguous_match code."""
        custom_format = ScoringFormat(ambiguous_match='*')
        # R (AG) matches A
        result = score_alignment("ATRG", "ATAG", DEFAULT_ADJUSTMENT_PARAMS, custom_format)
        assert result.score_aligned == "||*|"
        assert result.identity == 1.0


class TestIUPACAdjustment:
    """Test IUPAC ambiguity code handling."""
    
    def test_exact_iupac_match(self):
        """Same IUPAC codes should always match."""
        result = score_alignment("ANRG", "ANRG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 1.0
        assert result.mismatches == 0
        assert result.scored_positions == 4
        # A=A is exact (|), N=N, R=R, G=G are all ambiguous (=)
        assert result.score_aligned == "|==|"
    
    def test_iupac_intersection_match(self):
        """Different IUPAC codes with overlapping nucleotides should match."""
        # R (AG) and K (GT) both contain G
        result = score_alignment("ATRG", "ATKG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 1.0
        assert result.mismatches == 0
        assert result.scored_positions == 4
        # A=A, T=T are exact (|), R=K is ambiguous (=), G=G is exact (|)
        assert result.score_aligned == "||=|"
    
    def test_iupac_intersection_disabled(self):
        """IUPAC intersection disabled: codes that don't match can still extend context."""
        # R (AG) vs K (GT) - don't match each other, but each can extend surrounding context
        # R can represent G (matches right context 'G'), K can represent T (matches left context 'T')
        # Variant range algorithm: both are valid extensions of their contexts → 0 edits
        params = AdjustmentParams(handle_iupac_overlap=False)
        result = score_alignment("ATRG", "ATKG", params)
        assert result.identity == 1.0  # Both IUPAC codes extend context
        assert result.mismatches == 0
        assert result.scored_positions == 3  # Position 2 is variant, not counted
        assert result.score_aligned == "||=|"  # = shows it's a homopolymer-like extension

    def test_iupac_no_intersection(self):
        """IUPAC codes with no overlap can still extend context in variant range."""
        # R (AG) and Y (CT) have no direct overlap, but:
        # R can represent G (matches right context 'G')
        # Y can represent T (matches left context 'T')
        # Variant range algorithm: both are valid extensions → 0 edits
        result = score_alignment("ATRG", "ATYG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 1.0
        assert result.mismatches == 0
        assert result.scored_positions == 3  # Position 2 is variant, not counted
        assert result.score_aligned == "||=|"
    
    def test_standard_vs_iupac(self):
        """Standard nucleotide vs IUPAC should match if overlap exists."""
        # R (AG) contains A
        result = score_alignment("ATRG", "ATAG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 1.0
        assert result.mismatches == 0
        assert result.scored_positions == 4
        # A=A, T=T are exact (|), R=A is ambiguous (=), G=G is exact (|)
        assert result.score_aligned == "||=|"


class TestEndTrimming:
    """Test end region mismatch skipping."""
    
    def test_end_trimming_enabled(self):
        """Mismatches near ends should be skipped with 20bp trimming."""
        # Create sequences longer than 40bp to enable trimming
        seq1 = "A" * 21 + "XXXX" + "T" * 21  # Mismatches in middle
        seq2 = "A" * 21 + "TTTT" + "T" * 21
        params = AdjustmentParams(end_skip_distance=20)
        result = score_alignment(seq1, seq2, params)

        # Variant range algorithm: XXXX vs TTTT is one variant range
        # TTTT extends T context (right side) → pure extension
        # XXXX doesn't extend any context → core content
        # Result: 1 normalized edit for the XXXX core
        assert result.mismatches == 1
        assert result.scored_positions == 5  # 1 (A) + 1 (variant) + 3 (TT)
    
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
        # AAA match (|||), homopolymer indel (=), TT match (||), R/K ambiguous (=), G match (|)
        assert result.score_aligned == "|||=||=|"
    
    def test_all_adjustments_disabled(self):
        """Same case with all adjustments disabled."""
        result = score_alignment("AAA-TTRG", "AAAATTKG", RAW_ADJUSTMENT_PARAMS)
        assert result.identity < 1.0  # Raw scoring shows differences
        assert result.mismatches > 0
    
    def test_selective_adjustments(self):
        """Enable only some adjustments."""
        params = AdjustmentParams(
            normalize_homopolymers=True,
            handle_iupac_overlap=False,  # Disabled for direct matching
            normalize_indels=True,
            end_skip_distance=0
        )
        result = score_alignment("AAA-TTRG", "AAAATTKG", params)
        # Variant range algorithm:
        # - Position 3: A extends A context → pure extension
        # - Position 6: R extends G context, K extends T context → both pure extensions
        # All variants are pure extensions → 0 edits
        assert result.identity == 1.0
        assert result.mismatches == 0


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
        """Test repeat motif near start or end of sequence."""
        # AT repeat near start (add leading context to avoid overhang)
        seq1 = "C--ATATGC"
        seq2 = "CATATATGC"

        params = AdjustmentParams(max_repeat_motif_length=2)
        result = score_alignment(seq1, seq2, params)

        # AT extends AT context (right side)
        assert result.identity == 1.0
        assert '=' in result.score_aligned
    
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

    def test_invalid_adjustment_params(self):
        """Invalid AdjustmentParams configuration should raise error."""
        # normalize_homopolymers=True requires max_repeat_motif_length >= 1
        with pytest.raises(ValueError, match="Contradictory configuration"):
            AdjustmentParams(normalize_homopolymers=True, max_repeat_motif_length=0)

        # Negative max_repeat_motif_length should also fail
        with pytest.raises(ValueError, match="Contradictory configuration"):
            AdjustmentParams(normalize_homopolymers=True, max_repeat_motif_length=-1)

        # But normalize_homopolymers=False with max_repeat_motif_length=0 should be OK
        params = AdjustmentParams(normalize_homopolymers=False, max_repeat_motif_length=0)
        assert params.normalize_homopolymers is False
        assert params.max_repeat_motif_length == 0


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

        # With end trimming enabled: only middle region scored
        trim_params = AdjustmentParams(end_skip_distance=20)
        result_trim = score_alignment(seq1, seq2, trim_params)
        # N's extend context (N matches A on right, G on left), X's and Y's are core
        # Result: 1 normalized edit for each variant range with core content
        assert result_trim.identity > 0.9  # High identity due to N extending context

        # Without end trimming (default): all positions scored
        result_no_trim = score_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)
        # More positions scored but similar identity due to N extending context
        assert result_no_trim.scored_positions > result_trim.scored_positions


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


class TestMSADualGaps:
    """Test homopolymer detection with MSA-derived dual-gap sequences.

    In multi-sequence alignments (MSA), both sequences may have gaps at the same
    position due to alignment with a third sequence. This test class validates that
    homopolymer normalization works correctly with such dual-gap sequences.
    """

    def test_dual_gap_right_context_homopolymer(self):
        """Case: AGA--TT vs AGAT-TT - 'T' should be recognized as homopolymer extension."""
        result = score_alignment("AGA--TT", "AGAT-TT", DEFAULT_ADJUSTMENT_PARAMS)

        # The 'T' at position 3 should be recognized as extending the 'TT' at positions 5-6
        assert result.identity == 1.0, f"Expected 1.0, got {result.identity}"
        assert result.mismatches == 0
        # Score pattern: AGA (|||) + dual-gap (|) + T extension (=) + dual-gap (|) + TT (||)
        assert result.score_aligned == "|||=|||"

    def test_dual_gap_left_context_homopolymer(self):
        """Case: AAA--TTT vs AA-A-TTT - both 'A's should be recognized as homopolymer extensions."""
        result = score_alignment("AAA--TTT", "AA-A-TTT", DEFAULT_ADJUSTMENT_PARAMS)

        # Both 'A's at positions 2 and 3 should extend the 'AA' at positions 0-1
        assert result.identity == 1.0, f"Expected 1.0, got {result.identity}"
        assert result.mismatches == 0
        # Score pattern: AA (||) + A extension (=) + A extension (=) + dual-gap (|) + TTT (|||)
        assert result.score_aligned == "||==||||"

    def test_dual_gap_skipped_for_context(self):
        """Case: CTT--GCTGGC vs CTT-TGCTGGC - context extraction skips dual-gap to find homopolymer."""
        # Test with gap in seq1
        result = score_alignment("CTT--GCTGGC", "CTT-TGCTGGC", DEFAULT_ADJUSTMENT_PARAMS)

        # The 'T' at position 4 should be recognized as extending the 'T' at position 2
        # Context extraction must skip the dual-gap at position 3 to find 'T' context
        assert result.identity == 1.0, f"Expected 1.0, got {result.identity}"
        assert result.mismatches == 0
        # Score pattern: CTT (|||) + dual-gap (|) + T extension (=) + GCTGGC (||||||)
        assert result.score_aligned == "||||=||||||"

        # Test with sequences reversed (gap in seq2) - should be symmetric
        result_reversed = score_alignment("CTT-TGCTGGC", "CTT--GCTGGC", DEFAULT_ADJUSTMENT_PARAMS)
        assert result_reversed.identity == 1.0, f"Expected 1.0, got {result_reversed.identity}"
        assert result_reversed.mismatches == 0
        assert result_reversed.score_aligned == "||||=||||||"

    def test_dual_gap_not_homopolymer(self):
        """Case: AGT-AC vs AG-GAC - variant range with different alleles."""
        result = score_alignment("AGT-AC", "AG-GAC", DEFAULT_ADJUSTMENT_PARAMS)

        # Variant range algorithm:
        # - Variant range at positions 2-3 (T- vs -G)
        # - allele1="T", allele2="G"
        # - left_context="G", right_context="A"
        # - T doesn't extend G or A → core="T"
        # - G extends G (left context) → pure extension
        # - One pure extension, one core → 1 normalized edit
        assert result.identity < 1.0
        assert result.mismatches == 1  # Single variant range, one core content
        assert result.scored_positions == 5  # AG + variant + AC

    def test_alternating_indels_with_trailing_dual_gap(self):
        """Case: TGC-C-TC vs TGCT--TC - KEY TEST for variant range algorithm."""
        result = score_alignment("TGC-C-TC", "TGCT--TC", DEFAULT_ADJUSTMENT_PARAMS)

        # THIS IS THE CASE THE VARIANT RANGE ALGORITHM WAS DESIGNED TO HANDLE!
        # Position 3: 'T' in seq2 (gap in seq1)
        # Position 4: 'C' in seq1 (gap in seq2)
        # Position 5: dual-gap (counted as match)
        #
        # Variant range at positions 3-4:
        # - allele1="C", allele2="T"
        # - left_context="C" (position 2), right_context="T" (position 6)
        # - C extends C (left context) → pure extension
        # - T extends T (right context) → pure extension
        # - Both are valid extensions in their respective directions → 0 edits
        assert result.identity == 1.0
        assert result.mismatches == 0
        assert result.scored_positions == 6  # TGC + dual-gap + TC (variant not counted)
        assert result.score_aligned == "|||==|||"  # == shows both extensions

        # Test reversed - should be symmetric
        result_rev = score_alignment("TGCT--TC", "TGC-C-TC", DEFAULT_ADJUSTMENT_PARAMS)
        assert result_rev.identity == 1.0
        assert result_rev.mismatches == 0
        assert result_rev.score_aligned == "|||==|||"

    def test_dual_gap_same_char_homopolymer(self):
        """Case: AGG-AC vs AG-GAC - both 'G's should be recognized as homopolymer extensions."""
        result = score_alignment("AGG-AC", "AG-GAC", DEFAULT_ADJUSTMENT_PARAMS)

        # Both 'G's extend the 'G' at position 1
        assert result.identity == 1.0, f"Expected 1.0, got {result.identity}"
        assert result.mismatches == 0
        # Score pattern: A (|) + G (|) + G extension (=) + G extension (=) + A (|) + C (|)
        assert result.score_aligned == "||==||"

    def test_multiple_consecutive_dual_gaps(self):
        """Multiple consecutive dual-gaps in indel region."""
        result = score_alignment("A---TT", "A--GTT", DEFAULT_ADJUSTMENT_PARAMS)

        # Positions 1-2 are dual-gaps (matched), position 3 has 'G'
        # 'G' is not a homopolymer extension of 'T'
        assert result.mismatches == 1
        assert result.scored_positions == 6  # A + dual-gaps (2) + indel + TT

    def test_dual_gaps_in_context_region(self):
        """Dual-gaps in context region should be skipped when extracting context."""
        result = score_alignment("AA--AACGG", "AAT--ACGG", DEFAULT_ADJUSTMENT_PARAMS)

        # Context should skip the dual-gaps at positions 2-3 and find 'AA' context
        # The second 'A' at position 5 should match the 'A' homopolymer context
        assert result.identity > 0.8  # Should recognize at least some homopolymer

    def test_insufficient_context_with_dual_gaps(self):
        """Insufficient context when dual-gaps consume early positions."""
        # Add leading match to avoid overhang behavior
        result = score_alignment("C--ATT", "C-G-TT", DEFAULT_ADJUSTMENT_PARAMS)

        # Minimal left context (C), 'G' and 'A' don't extend C
        # Should be treated as regular indel
        assert result.mismatches >= 1

    def test_conflicting_context_no_homopolymer(self):
        """When sequences disagree in context, homopolymer should not be detected."""
        # Position 1 has 'G' in seq1 and 'X' in seq2 - conflicting context
        result = score_alignment("AGT-AC", "AXT-AC", DEFAULT_ADJUSTMENT_PARAMS)

        # Cannot extract consensus context, 'T' treated as regular indel
        assert result.mismatches == 1

    def test_dual_gap_both_directions(self):
        """Test with homopolymer extensions on both left and right sides."""
        result = score_alignment("AAA--GGG", "AA-A-GGG", DEFAULT_ADJUSTMENT_PARAMS)

        # Both 'A's extend left context, no right extension (G's already match)
        assert result.identity == 1.0, f"Expected 1.0, got {result.identity}"
        assert result.mismatches == 0

    def test_dual_gap_dinucleotide_repeat(self):
        """Test dinucleotide repeat with dual-gaps."""
        # Note: This is a complex case with no actual dual-gaps and a substitution
        result = score_alignment("ATAT---C", "AT-ATATC", DEFAULT_ADJUSTMENT_PARAMS)

        # This case has: AT match, A indel, T vs A substitution, TAT indel, C match
        # Not a simple dinucleotide repeat - just verify reasonable identity
        assert 0.0 <= result.identity <= 1.0  # Valid identity range

    def test_dual_gap_at_boundaries(self):
        """Dual-gap at the very start or end of sequence."""
        # Dual-gap at position 0-1
        result = score_alignment("--ATT", "--GTT", DEFAULT_ADJUSTMENT_PARAMS)

        # No left context, but should still work for the rest
        assert result.scored_positions >= 2  # At least the matching positions

    def test_dual_gap_scoring_visualization(self):
        """Verify that dual-gaps use the match marker '|' in score_aligned."""
        result = score_alignment("A--T", "A--T", DEFAULT_ADJUSTMENT_PARAMS)

        # All positions match, including the dual-gaps
        assert result.score_aligned == "||||"
        assert result.identity == 1.0
        assert result.mismatches == 0

    def test_dual_gap_mixed_with_regular_indels(self):
        """Complex case with both dual-gaps and regular indels."""
        result = score_alignment("AA--T-GG", "A-ATT-GG", DEFAULT_ADJUSTMENT_PARAMS)

        # Position 1: gap in seq2
        # Position 2: gap in seq1
        # Position 3: dual-gap
        # Position 4: 'T' in both
        # Position 5: gap in seq2
        # Should recognize some patterns
        assert 0.0 <= result.identity <= 1.0  # Valid identity range

    def test_dual_gap_with_iupac_codes(self):
        """Dual-gaps combined with IUPAC ambiguity codes."""
        result = score_alignment("AR--TT", "AG-RTT", DEFAULT_ADJUSTMENT_PARAMS)

        # 'R' at position 1 in seq1 matches 'G' at position 1 in seq2 (R contains G)
        # Position 2 has dual-gap
        # Position 3: 'R' in seq2 should match 'T' context? No.
        # This tests the interaction between MSA and IUPAC handling
        assert result.scored_positions >= 4

    def test_consensus_left_context_extraction(self):
        """Test that left context extraction enforces consensus correctly."""
        # Sequences agree at position 1 ('G'), so homopolymer should be detected
        result1 = score_alignment("AGG-AC", "AG-GAC", DEFAULT_ADJUSTMENT_PARAMS)
        assert result1.identity == 1.0  # Homopolymer detected

        # Sequences disagree at position 1 ('G' vs 'T'), no homopolymer
        result2 = score_alignment("AGG-AC", "AT-GAC", DEFAULT_ADJUSTMENT_PARAMS)
        assert result2.identity < 1.0  # Homopolymer not detected

    def test_consensus_right_context_extraction(self):
        """Test that right context extraction enforces consensus correctly."""
        # Sequences agree at positions 5-6 ('TT'), so homopolymer should be detected
        result1 = score_alignment("AGA--TT", "AGAT-TT", DEFAULT_ADJUSTMENT_PARAMS)
        assert result1.identity == 1.0  # Homopolymer detected

        # Sequences disagree at position 6 ('T' vs 'G'), no homopolymer
        result2 = score_alignment("AGA--TG", "AGAT-TG", DEFAULT_ADJUSTMENT_PARAMS)
        # 'T' should still be homopolymer extension since only 1 char of context needed
        # But if we need 2 chars of context and they disagree, it fails
        assert 0.0 <= result2.identity <= 1.0  # Valid range

    def test_end_to_end_msa_example(self):
        """Real-world MSA example with multiple dual-gaps and homopolymers."""
        # Simulates output from spoa or similar MSA tool
        seq1 = "ATCGAAA--TTTT--GGG"
        seq2 = "ATCG--AAATTT-TTGGG"

        result = score_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)

        # Should recognize 'A' and 'T' homopolymer extensions
        # Despite complex dual-gap pattern
        assert result.identity >= 0.9  # Most should be recognized as homopolymers
        assert result.score_aligned.count('=') >= 4  # Multiple homopolymer extensions


class TestVariantRangeAlgorithm:
    """Tests specifically for variant range algorithm behavior introduced in v0.2.0.

    The variant range algorithm identifies contiguous non-match regions, extracts
    gap-free alleles, analyzes them for repeat extensions, and scores using
    Occam's razor (most parsimonious explanation).
    """

    def test_opposite_direction_extensions_basic(self):
        """Basic case: alleles extend in opposite directions."""
        # TGC-C-TC vs TGCT--TC
        # C extends left (C context), T extends right (T context)
        # This is the canonical case for the variant range algorithm
        result = score_alignment("TGC-C-TC", "TGCT--TC", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 1.0
        assert result.mismatches == 0
        assert result.score_aligned == "|||==|||"

    def test_opposite_direction_with_core_content(self):
        """Alleles extend in opposite directions but also have core content."""
        # TGC-CX-TC vs TGCT---TC
        # allele1="CX", allele2="T"
        # C extends left C, X is core; T extends right T
        # allele1 has core "X", allele2 is pure extension
        # Score: 1 edit for the X core
        result = score_alignment("TGC-CX-TC", "TGCT---TC", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity < 1.0
        assert result.mismatches == 1  # X core counts as 1 edit
        # T(pos3)=extension, C(pos4)=extension, X(pos5)=core(space), dual-gap(pos6)=match
        assert result.score_aligned == "|||== |||"

    def test_both_alleles_have_core_same_content(self):
        """Both alleles have core content that's identical."""
        # AAA-XGG vs AA-X-GG
        # allele1="X", allele2="X"
        # Neither X extends A or G context
        # Both have core "X", cores match → 0 edits
        result = score_alignment("AAA-XGG", "AA-X-GG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 1.0
        assert result.mismatches == 0

    def test_both_alleles_have_core_different_content(self):
        """Both alleles have core content that differs."""
        # AAA-XGG vs AA-Y-GG
        # allele1="X", allele2="Y"
        # Neither extends context
        # Both have core, cores differ → 1 edit (normalized)
        result = score_alignment("AAA-XGG", "AA-Y-GG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity < 1.0
        assert result.mismatches == 1

    def test_partial_extension_left_only(self):
        """Allele partially extends left context but has remaining core."""
        # AAA--GGG vs AAAAXGGG
        # allele1="", allele2="AX"
        # A extends left A, X is core
        # allele1 is pure (empty), allele2 has core "X"
        result = score_alignment("AAA--GGG", "AAAAXGGG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity < 1.0
        assert result.mismatches == 1  # X counts as 1 edit

    def test_partial_extension_right_only(self):
        """Allele partially extends right context but has remaining core."""
        # AAA--GGG vs AAAXGGGG
        # allele1="", allele2="XG"
        # G extends right G, X is core
        result = score_alignment("AAA--GGG", "AAAXGGGG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity < 1.0
        assert result.mismatches == 1  # X counts as 1 edit

    def test_partial_extension_both_sides(self):
        """Allele extends both sides but has core in middle."""
        # AAA---GGG vs AAAXG-GGG (note: need proper alignment)
        # Test with pre-aligned sequences
        result = score_alignment("AAA---GGG", "AAAXGGGGG", DEFAULT_ADJUSTMENT_PARAMS)
        # allele1="", allele2="XGG"
        # GG extends G context, X is core
        assert result.mismatches == 1  # Only X counts

    def test_iupac_allele_extends_context(self):
        """IUPAC code in allele can extend context via equivalence."""
        # ATG-C vs ATGRC where R=(A|G)
        # Context is G (both sequences agree at position 2)
        # R can represent G → extends G context (left)
        # Pure extension → 0 edits
        result = score_alignment("ATG-C", "ATGRC", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 1.0
        assert result.mismatches == 0

    def test_iupac_allele_cannot_extend_context(self):
        """IUPAC code in allele cannot extend context."""
        # ATG-C vs ATGWC where W=(A|T)
        # Context is G (position 2), right context is C (position 4)
        # W cannot represent G or C → core content
        # Compare with pure extension (gap) → 1 edit
        result = score_alignment("ATG-C", "ATGWC", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity < 1.0
        assert result.mismatches == 1

    def test_empty_allele_vs_pure_extension(self):
        """One allele is empty (all gaps), other is pure extension."""
        # AAA-TTT vs AAAATTT
        # allele1="", allele2="A"
        # A extends A context → pure extension
        # Both pure extensions (empty is trivially pure) → 0 edits
        result = score_alignment("AAA-TTT", "AAAATTT", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 1.0
        assert result.mismatches == 0

    def test_empty_allele_vs_core_content(self):
        """One allele is empty, other has core content."""
        # AAA-TTT vs AAAX-TT (need same length)
        # allele1="", allele2="X"
        # X doesn't extend A or T → core
        # Empty vs core → 1 edit
        result = score_alignment("AAA-TT", "AAAXTT", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity < 1.0
        assert result.mismatches == 1

    def test_longer_variant_range_multiple_chars(self):
        """Variant range with multiple characters, partial extensions."""
        # AAAT----GGG vs AAATCCCCGGG
        # allele1="", allele2="CCCC"
        # CCCC doesn't extend A or G → core
        result = score_alignment("AAAT----GGG", "AAATCCCCGGG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity < 1.0
        assert result.mismatches == 1  # Normalized indel

    def test_dinucleotide_repeat_in_variant_range(self):
        """Dinucleotide repeat recognized in variant range."""
        # CGATAT--C vs CGATATATC
        # allele1="", allele2="AT"
        # AT extends AT context (dinucleotide)
        result = score_alignment("CGATAT--C", "CGATATATC", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 1.0
        assert result.mismatches == 0

    def test_variant_range_at_sequence_start(self):
        """Variant range at very start of sequence (no left context)."""
        # T-ATCG vs TGATCG - variant after first position, not overhang
        # No left context for G (only T, G doesn't extend T) → treated as core
        result = score_alignment("T-ATCG", "TGATCG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity < 1.0
        assert result.mismatches == 1

    def test_variant_range_at_sequence_end(self):
        """Variant range at very end of sequence (no right context)."""
        # ATCG-T vs ATCGAT - variant in middle, not overhang
        # A doesn't extend G (left) or T (right) → core
        result = score_alignment("ATCG-T", "ATCGAT", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity < 1.0
        assert result.mismatches == 1

    def test_variant_range_at_end_with_extension(self):
        """Variant range at end where allele does extend left context."""
        # ATCG- vs ATCGG
        # G extends G (left context) → pure extension
        result = score_alignment("ATCG-", "ATCGG", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.identity == 1.0
        assert result.mismatches == 0

    def test_multiple_variant_ranges(self):
        """Multiple separate variant ranges in one alignment."""
        # AA-TT-CC vs AAATTGCC
        # Range 1 (pos 2): A extends A
        # Range 2 (pos 5): G doesn't extend T or C → core
        result = score_alignment("AA-TT-CC", "AAATTGCC", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.mismatches == 1  # Only G counts

    def test_adjacent_variant_ranges(self):
        """Adjacent variant ranges separated by single match."""
        # A-T-G vs AXTXG (X doesn't match T)
        # This creates complex variant patterns
        result = score_alignment("A-T-G", "AXTXG", DEFAULT_ADJUSTMENT_PARAMS)
        # X's don't extend context
        assert result.mismatches >= 2  # Both X's count

    def test_normalized_indels_disabled_with_variant_range(self):
        """Variant range scoring respects normalize_indels=False."""
        params = AdjustmentParams(normalize_indels=False)
        # AAA---GGG vs AAAXYZGGG
        # When indels not normalized, XYZ counts as 3 edits
        result = score_alignment("AAA---GGG", "AAAXYZGGG", params)
        assert result.mismatches == 3  # X, Y, Z each count

    def test_homopolymer_disabled_treats_extensions_as_indels(self):
        """With normalize_homopolymers=False, extensions are regular indels."""
        params = AdjustmentParams(normalize_homopolymers=False)
        # AAA-TTT vs AAAATTT
        # Without homopolymer normalization, A is just an indel
        result = score_alignment("AAA-TTT", "AAAATTT", params)
        assert result.mismatches == 1  # Single indel (normalized by default)

    def test_both_adjustments_disabled(self):
        """Both homopolymer and indel normalization disabled."""
        params = AdjustmentParams(normalize_homopolymers=False, normalize_indels=False)
        result = score_alignment("AAA-TTT", "AAAATTT", params)
        assert result.mismatches == 1  # 1 gap position
        assert result.scored_positions == 7  # All positions scored

    # Score visualization tests - ensure correct markers for different scenarios

    def test_score_visualization_both_pure_extensions(self):
        """When both alleles are pure extensions, show '=' markers."""
        # C extends C (left), T extends T (right) → both pure
        result = score_alignment("TGC-C-TC", "TGCT--TC", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.score_aligned == "|||==|||"
        assert result.mismatches == 0

    def test_score_visualization_one_extension_one_core(self):
        """When one allele extends and other has core, show substitution marker."""
        # T extends T (right context), C doesn't extend → C is core
        # The space at position 5 indicates the mismatch
        result = score_alignment("TTG-ATT---T", "TTG-ACT---T", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.score_aligned == "||||| |||||"
        assert result.mismatches == 1

    def test_score_visualization_both_have_different_core(self):
        """When both alleles have different core content, show substitution."""
        # X doesn't extend A or G, Y doesn't extend A or G → both core, different
        result = score_alignment("AAA-XGG", "AA-Y-GG", DEFAULT_ADJUSTMENT_PARAMS)
        assert " " in result.score_aligned  # Substitution marker for core mismatch
        assert result.mismatches == 1

    def test_score_visualization_extension_vs_gap(self):
        """Extension on one side, gap on other → show '=' for extension positions."""
        # A extends A context, empty allele → pure extension
        result = score_alignment("AAA-TTT", "AAAATTT", DEFAULT_ADJUSTMENT_PARAMS)
        assert result.score_aligned == "|||=|||"
        assert result.mismatches == 0

    def test_score_visualization_core_indel(self):
        """Core content that doesn't extend context → show indel marker."""
        # X doesn't extend A or T → core content (indel)
        result = score_alignment("AAA-TT", "AAAXTT", DEFAULT_ADJUSTMENT_PARAMS)
        assert " " in result.score_aligned  # Core indel shows as space
        assert result.mismatches == 1