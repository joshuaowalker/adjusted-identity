#!/usr/bin/env python3
"""
Specification compliance tests for the v0.2.x variant range algorithm.

These tests are harvested from empirical analysis of real-world MSA data.
Each test represents a distinct pattern found in the data and validates
that the algorithm behaves according to the specification.

Test cases reference sections in docs/SCORING_SPEC.md.
"""

import pytest
from adjusted_identity import score_alignment, DEFAULT_ADJUSTMENT_PARAMS, AdjustmentParams


class TestDualGapHandling:
    """Tests for dual-gap (MSA artifact) handling.

    Spec reference: Section 2.2 (Dual-Gap Handling), Section 7.1.1

    Dual-gap positions (both sequences have '-') are:
    - Treated as matches for variant range boundary detection
    - Marked with '.' in visualization
    - NOT counted in scored_positions
    """

    def test_dual_gap_not_scored(self):
        """Spec 7.1.1: Dual-gaps are not counted in scored_positions.

        Harvested from: ONT10.82-B11 FASTA pairwise comparison
        """
        seq1 = "AA--TT"
        seq2 = "AA--TT"
        result = score_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)

        assert result.score_aligned == "||..||"
        assert result.scored_positions == 4  # Only A, A, T, T
        assert result.identity == 1.0
        assert result.mismatches == 0

    def test_dual_gap_visualization(self):
        """Spec 7.1.1: Dual-gaps are marked with '.' in score_aligned.

        Harvested from: ONT10.82-B11 FASTA pairwise comparison
        """
        seq1 = "A-T"
        seq2 = "A-T"
        result = score_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)

        assert result.score_aligned == "|.|"
        assert "." in result.score_aligned

    def test_dual_gap_does_not_split_variant_range(self):
        """Spec 2.2: Dual-gaps treated as matches for boundary detection.

        A dual-gap within a variant range should NOT split it into two ranges.

        Harvested from: ONT10.82-B11 FASTA pairwise comparison
        """
        seq1 = "AA--T-CC"
        seq2 = "AAT---CC"
        result = score_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)

        # The T positions should be in the same variant range
        # Both T's extend the same direction -> pure extensions
        assert result.mismatches == 0
        assert result.identity == 1.0


class TestExtensionRecognition:
    """Tests for homopolymer/repeat extension recognition.

    Spec reference: Section 1.3 (Extension Detection), Section 5.1

    Extensions are only detected from positions of agreement (match boundaries).
    """

    def test_opposite_direction_extensions(self):
        """Spec 1.3, 1.5: Both alleles extend in opposite directions.

        This is the canonical case the variant range algorithm handles:
        - C extends left (C context)
        - T extends right (T context)
        - Both are pure extensions -> 0 edits

        Harvested from: ONT10.82-B11 FASTA, alternating gap pattern
        """
        seq1 = "TGC-C-TC"
        seq2 = "TGCT--TC"
        result = score_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)

        assert result.identity == 1.0
        assert result.mismatches == 0
        assert result.scored_positions == 5  # TGC + TC, variant range not scored

    def test_homopolymer_extension_same_direction(self):
        """Spec 1.3: Both alleles extend in same direction (left).

        Harvested from: ONT10.82-B11 FASTA pairwise comparison
        """
        seq1 = "AAG-G-CC"
        seq2 = "AAGG--CC"
        result = score_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)

        # Both G's extend the G context to the left
        assert result.identity == 1.0
        assert result.mismatches == 0

    def test_extension_from_match_boundary_only(self):
        """Spec 1.3: Extensions only detected from positions of agreement.

        The context for extension detection must come from matched positions
        outside the variant range, not from within variant ranges.

        Harvested from: ONT10.82-B11 FASTA, position 588 analysis
        """
        seq1 = "GG-ACG-TAG"
        seq2 = "GG-CGC-CAG"
        result = score_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)

        # A/C and G/C are substitutions, not extensions
        # Context G doesn't allow C to extend
        assert result.mismatches >= 2

    def test_gap_free_substitution_as_extensions(self):
        """Spec 8.7: Gap-free substitution recognized as opposite-direction extensions.

        Standard aligners don't implement homopolymer-aware scoring, so they may
        produce gap-free alignments where a "substitution" actually represents
        equivalent homopolymer expansions. If the aligner had homopolymer scoring,
        it might produce ATT-CA vs AT-CCA instead.

        Harvested from: ONT10.82-B11 FASTA, larger dataset analysis
        """
        seq1 = "ATTCA"
        seq2 = "ATCCA"
        result = score_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)

        # T extends left T context (homopolymer)
        # C extends right C context (homopolymer)
        # Both are pure extensions -> 0 edits
        assert result.identity == 1.0
        assert result.mismatches == 0
        assert result.score_aligned == "||=||"


class TestCoreContentScoring:
    """Tests for core content (non-extension) scoring.

    Spec reference: Section 1.4 (Core Content), Section 5.2 (Core Comparison)
    """

    def test_both_alleles_have_core(self):
        """Spec 5.2: Both alleles have core content that differs.

        When neither allele is a pure extension, compare cores.

        Harvested from: ONT10.82-B11 FASTA pairwise comparison
        """
        seq1 = "AAA-X-TTT"
        seq2 = "AAA--YTTT"
        result = score_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)

        # X doesn't extend A or T -> core
        # Y doesn't extend A or T -> core
        # Different cores -> mismatch
        assert result.mismatches >= 1

    def test_one_pure_extension_one_core(self):
        """Spec 5.1: One allele is pure extension, other has core.

        Harvested from: ONT10.82-B11 FASTA pairwise comparison
        """
        seq1 = "AAA-A-TTT"  # A extends A context
        seq2 = "AAA--XTTT"  # X doesn't extend A or T
        result = score_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)

        # A is pure extension, X is core
        # Core content X counts as mismatch
        assert result.mismatches == 1

    def test_substitution_not_extension(self):
        """Spec 1.3: Substitutions are not treated as extensions.

        When nucleotide doesn't match context, it's core content not extension.

        Harvested from: ONT10.82-B11 FASTA, regression analysis position 588
        """
        seq1 = "GG-T-AA"
        seq2 = "GG-C-AA"
        result = score_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)

        # T doesn't extend G or A -> core
        # C doesn't extend G or A -> core
        # T != C -> substitution
        assert result.mismatches == 1

    def test_insertion_not_extending_context(self):
        """Spec 9.2: Insertion not extending context is correctly counted as mismatch.

        When an inserted nucleotide doesn't match the surrounding context,
        it should be counted as a mismatch, not treated as an extension.

        Harvested from: edlib alignment comparison, FASTQ pairwise analysis
        """
        # G insertion between C (left context) and A (right context)
        # G doesn't extend C or A, so it's a mismatch
        seq1 = "AAC-ATT"
        seq2 = "AACGATT"
        result = score_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)

        # G doesn't extend C (left) or A (right) → mismatch
        assert result.mismatches == 1
        assert " " in result.score_aligned  # Space indicates mismatch

        # Contrast: G insertion that DOES extend context
        seq1_ext = "AAG-ATT"
        seq2_ext = "AAGGATT"
        result_ext = score_alignment(seq1_ext, seq2_ext, DEFAULT_ADJUSTMENT_PARAMS)

        # G extends G (left context) → extension, 0 mismatches
        assert result_ext.mismatches == 0
        assert "=" in result_ext.score_aligned  # Extension marker


class TestVariantRangeBoundaries:
    """Tests for variant range boundary detection.

    Spec reference: Section 2 (Variant Range Detection)
    """

    def test_match_ends_variant_range(self):
        """Spec 2.1: Match position ends a variant range.

        Harvested from: ONT10.82-B11 FASTA pairwise comparison
        """
        seq1 = "AA-G-TT"
        seq2 = "AAG--TT"
        result = score_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)

        # Variant range is positions 2-4
        # G extends G context -> pure extension
        assert result.mismatches == 0

    def test_multiple_variant_ranges(self):
        """Spec 2: Multiple separate variant ranges are scored independently.

        Harvested from: ONT10.82-B11 FASTA pairwise comparison
        """
        seq1 = "AA-G-CC-T-GG"
        seq2 = "AAG--CCT--GG"
        result = score_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)

        # Two variant ranges, each with pure extensions
        assert result.mismatches == 0


class TestContextExtraction:
    """Tests for context extraction from match boundaries.

    Spec reference: Section 3 (Context Extraction)
    """

    def test_context_skips_dual_gaps(self):
        """Spec 3.1: Dual-gaps are skipped when extracting context.

        Harvested from: ONT10.82-B11 FASTA pairwise comparison
        """
        seq1 = "TT--T-GG"
        seq2 = "TT---TGG"
        result = score_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)

        # Context extraction should skip dual-gap to find T context
        # T extends T -> pure extension
        assert result.mismatches == 0

    def test_consensus_context_with_single_gap(self):
        """Spec 3.1: Use nucleotide when one sequence has gap.

        Harvested from: ONT10.82-B11 FASTA pairwise comparison
        """
        seq1 = "AT-T-GG"
        seq2 = "A-TT-GG"
        result = score_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)

        # Left context: T (consensus from positions with gaps)
        # T extends T -> pure extension
        assert result.mismatches == 0


class TestNormalizationParameters:
    """Tests for parameter interactions.

    Spec reference: Section 6 (Parameter Interactions)
    """

    def test_homopolymer_normalization_disabled(self):
        """Spec 6.1: With normalize_homopolymers=False, no extension detection.

        Harvested from: ONT10.82-B11 FASTA, parameter variation tests
        """
        params = AdjustmentParams(
            normalize_homopolymers=False,
            normalize_indels=True,
            handle_iupac_overlap=True,
            end_skip_distance=0,
        )
        seq1 = "AAA-TTT"
        seq2 = "AAAATTT"
        result = score_alignment(seq1, seq2, params)

        # Without HP normalization, the A is counted as indel not extension
        assert result.mismatches == 1  # One indel event

    def test_indel_normalization_disabled(self):
        """Spec 6.2: With normalize_indels=False, each position counts.

        Harvested from: ONT10.82-B11 FASTA, parameter variation tests
        """
        params = AdjustmentParams(
            normalize_homopolymers=False,
            normalize_indels=False,
            handle_iupac_overlap=True,
            end_skip_distance=0,
        )
        seq1 = "AAA---TTT"
        seq2 = "AAACCCTTT"
        result = score_alignment(seq1, seq2, params)

        # Each of 3 indel positions counts separately
        assert result.mismatches == 3


class TestVisualization:
    """Tests for score string visualization markers.

    Spec reference: Section 7 (Visualization)
    """

    def test_extension_marker(self):
        """Spec 7.1: Extension positions marked with '='.
        """
        seq1 = "AAA-TTT"
        seq2 = "AAAATTT"
        result = score_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)

        assert "=" in result.score_aligned

    def test_match_marker(self):
        """Spec 7.1: Match positions marked with '|'.
        """
        seq1 = "AAATTT"
        seq2 = "AAATTT"
        result = score_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)

        assert result.score_aligned == "||||||"

    def test_mismatch_marker(self):
        """Spec 7.1: Mismatch positions marked with ' ' (space).
        """
        seq1 = "AAATTT"
        seq2 = "AAACTT"
        result = score_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)

        assert " " in result.score_aligned
        assert result.mismatches == 1

    def test_iupac_ambiguous_match_marker(self):
        """Spec 8.8: IUPAC ambiguous matches marked with '='.

        When one sequence has an IUPAC ambiguity code (e.g., R={A,G}) and
        the other has a nucleotide within that set, it's an ambiguous match.
        """
        seq1 = "AAARTTT"  # R = A or G
        seq2 = "AAAGTTT"  # G is in {A, G}
        result = score_alignment(seq1, seq2, DEFAULT_ADJUSTMENT_PARAMS)

        assert result.score_aligned == "|||=|||"
        assert result.identity == 1.0
        assert result.mismatches == 0
        assert result.scored_positions == 7
