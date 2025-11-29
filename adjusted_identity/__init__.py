#!/usr/bin/env python3
"""
Copyright (c) 2025, Josh Walker

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Adjusted Identity Calculator for DNA Sequences

This module provides functions to calculate sequence identity metrics that
account for homopolymer length differences, commonly used in mycological
DNA barcoding applications.

Author: Josh Walker

Based on the MycoBLAST algorithm developed by Stephen Russell and Mycota Lab.
See: https://mycotalab.substack.com/p/why-ncbi-blast-identity-scores-can
"""

import edlib
import re
from dataclasses import dataclass
# Custom reverse complement implementation (replaces BioPython dependency)
# Uses optimized str.translate() for performance that exceeds BioPython
_RC_TRANSLATE_TABLE = str.maketrans(
    'ATGCRYSWKMBDHVNatgcryswkmbdhvn-',
    'TACGYRSWMKVHDBNtacgyrswmkvhdbn-'
)


def _reverse_complement(seq):
    """
    Generate reverse complement of DNA sequence with full IUPAC support.
    
    Handles all standard nucleotides (ATGC) and IUPAC ambiguity codes.
    Unknown characters are left unchanged. Uses optimized str.translate()
    for performance comparable to BioPython.
    
    Args:
        seq (str): DNA sequence to reverse complement
        
    Returns:
        str: Reverse complement sequence
        
    Examples:
        >>> _reverse_complement('ATCG')
        'CGAT'
        >>> _reverse_complement('ATCGRGTC')  # R = A/G, becomes Y = C/T
        'GACYCGAT'
    """
    return seq.translate(_RC_TRANSLATE_TABLE)[::-1]

@dataclass(frozen=True)
class AlignmentResult:
    """Result of sequence alignment with identity calculations.

    This dataclass contains alignment results with a single set of identity metrics
    based on the specified adjustment parameters. Use different AdjustmentParams
    to get raw vs adjusted results.

    Fields:
        identity: Identity score (0.0-1.0) based on adjustment parameters
        mismatches: Number of mismatches/edits counted
        scored_positions: Number of positions used for identity calculation (denominator)
        seq1_coverage: Fraction of seq1 covered by alignment (0.0-1.0)
        seq2_coverage: Fraction of seq2 covered by alignment (0.0-1.0)
        seq1_aligned: Aligned sequence 1 with gap characters
        seq2_aligned: Aligned sequence 2 with gap characters
        score_aligned: Scoring visualization for seq1 (kept for backwards compatibility)
        score_aligned_seq2: Scoring visualization for seq2
    """
    identity: float
    mismatches: int
    scored_positions: int
    seq1_coverage: float
    seq2_coverage: float
    seq1_aligned: str
    seq2_aligned: str
    score_aligned: str
    score_aligned_seq2: str = ""  # Default empty for backwards compatibility
    


@dataclass(frozen=True)
class ScoringFormat:
    """Format codes for alignment scoring visualization."""
    match: str = '|'                    # Exact match (A=A, C=C, G=G, T=T)
    ambiguous_match: str = '='          # Ambiguous nucleotide match (any IUPAC code match)
    substitution: str = ' '             # Nucleotide substitution
    indel_start: str = ' '              # First position of indel (scored)
    indel_extension: str = '-'          # Indel positions skipped due to normalization
    homopolymer_extension: str = '='    # Homopolymer length difference
    end_trimmed: str = '.'              # Position outside scoring region (end trimmed)
    
    def __post_init__(self):
        """Validate that all scoring codes are single characters."""
        for field_name, value in self.__dict__.items():
            if not isinstance(value, str) or len(value) != 1:
                raise ValueError(f"Scoring code '{field_name}' must be a single character, got: {value!r}")


@dataclass(frozen=True)
class AdjustmentParams:
    """
    Parameters for MycoBLAST-style sequence adjustments.
    
    Attributes:
        normalize_homopolymers: Ignore homopolymer length differences (e.g., "AAA" vs "AAAA")
        handle_iupac_overlap: Allow different ambiguity codes to match via nucleotide intersection
        normalize_indels: Count contiguous indels as single evolutionary events
        end_skip_distance: Number of nucleotides (not positions) to skip from each sequence end.
                          Only activates when sequences have ≥ 2×end_skip_distance nucleotides.
                          Set to 0 to disable end trimming completely.
        max_repeat_motif_length: Maximum length of repeat motifs to detect (e.g., 2 for dinucleotides).
                                Set to 1 to only detect homopolymers, 2 for dinucleotides, etc.
    """
    normalize_homopolymers: bool = True      # Ignore homopolymer length differences
    handle_iupac_overlap: bool = True        # Allow different ambiguity codes to match via intersection
    normalize_indels: bool = True            # Count contiguous indels as single events
    end_skip_distance: int = 0               # Nucleotides to skip from each end (0 = disabled by default)
    max_repeat_motif_length: int = 2         # Maximum repeat motif length to detect (default: dinucleotides)

    def __post_init__(self):
        """Validate parameter combinations."""
        if self.normalize_homopolymers and self.max_repeat_motif_length < 1:
            raise ValueError(
                f"Contradictory configuration: normalize_homopolymers=True requires "
                f"max_repeat_motif_length >= 1 to detect homopolymers, got {self.max_repeat_motif_length}"
            )


@dataclass(frozen=True)
class AlleleAnalysis:
    """Analysis of an allele's composition within a variant range.

    Used by the variant range algorithm to track how much of an allele
    can be explained as repeat extensions of surrounding context.

    Attributes:
        left_extension_count: Number of characters consumed as left repeat extension
        right_extension_count: Number of characters consumed as right repeat extension
        core_content: Remaining characters not explained by repeat extensions
        is_pure_extension: True if entire allele is repeat extensions (core is empty)
    """
    left_extension_count: int
    right_extension_count: int
    core_content: str
    is_pure_extension: bool


# Default adjustment parameters (all adjustments enabled)
DEFAULT_ADJUSTMENT_PARAMS = AdjustmentParams()

# Raw parameters (no adjustments - equivalent to traditional sequence identity)
RAW_ADJUSTMENT_PARAMS = AdjustmentParams(
    normalize_homopolymers=False,
    handle_iupac_overlap=False, 
    normalize_indels=False,
    end_skip_distance=0
)

# Default scoring format
DEFAULT_SCORING_FORMAT = ScoringFormat()

# IUPAC nucleotide ambiguity codes
IUPAC_CODES = {
    '-': {'-'},
    'A': {'A'},
    'T': {'T'},
    'C': {'C'},
    'G': {'G'},
    'R': {'A', 'G'},      # puRine
    'Y': {'C', 'T'},      # pYrimidine
    'S': {'G', 'C'},      # Strong (3 H bonds)
    'W': {'A', 'T'},      # Weak (2 H bonds)
    'K': {'G', 'T'},      # Keto
    'M': {'A', 'C'},      # aMino
    'B': {'C', 'G', 'T'}, # not A
    'D': {'A', 'G', 'T'}, # not C
    'H': {'A', 'C', 'T'}, # not G
    'V': {'A', 'C', 'G'}, # not T
    'N': {'A', 'C', 'G', 'T'}, # aNy
}


def _are_nucleotides_equivalent(nuc1, nuc2, enable_iupac_intersection=True):
    """
    Check if two nucleotides are equivalent according to IUPAC ambiguity codes.
    
    Args:
        nuc1 (str): First nucleotide (single character)
        nuc2 (str): Second nucleotide (single character)
        enable_iupac_intersection (bool): Allow different ambiguity codes to match via intersection
        
    Returns:
        tuple: (is_match, is_ambiguous) where:
            - is_match: True if nucleotides are equivalent
            - is_ambiguous: True if match involves ambiguity codes (not exact A=A, C=C, G=G, T=T)
    """
    # Convert to uppercase
    nuc1 = nuc1.upper()
    nuc2 = nuc2.upper()
    
    # Check for exact match first (including gap characters)
    if nuc1 == nuc2:
        # Exact match - check if it's one of the standard nucleotides or a gap
        # Gaps matching are treated as non-ambiguous matches (MSA support)
        is_standard_or_gap = nuc1 in {'A', 'T', 'C', 'G', '-'}
        return (True, not is_standard_or_gap)
    
    # Get possible nucleotides for each code
    possible1 = IUPAC_CODES.get(nuc1, {nuc1})
    possible2 = IUPAC_CODES.get(nuc2, {nuc2})
    
    # Check if there's any overlap
    has_overlap = bool(possible1 & possible2)
    both_codes_sets = len(possible1) > 1 and len(possible2) > 1

    if enable_iupac_intersection:
        # Any overlap counts as a match, and it's ambiguous since nuc1 != nuc2
        return (has_overlap, has_overlap)
    else:
        # When disabled, only allow exact code match or standard nucleotide vs ambiguity code
        is_match = has_overlap and not both_codes_sets
        # It's ambiguous if there's a match but not exact (nuc1 != nuc2)
        return (is_match, is_match)



def _parse_suffix_gap_from_cigar(cigar_string):
    """
    Parse CIGAR string from the end to find suffix gaps.
    
    Args:
        cigar_string: CIGAR string from edlib (e.g., "6=3D")

    Returns:
        tuple: (gap_length, gap_in_query) or None if no suffix gap
        gap_in_query=True means query has gaps (target extends beyond query)
    """
    if not cigar_string:
        return None
    
    # Parse CIGAR operations (length, operation_type)
    cigar_ops = re.findall(r'(\d+)([=XIDM])', cigar_string)
    if not cigar_ops:
        return None
    
    # Scan from end looking for contiguous gaps
    gap_length = 0
    gap_type = None  # 'D' for deletions (gaps in query), 'I' for insertions (gaps in target)
    
    for length_str, op in reversed(cigar_ops):
        if op in ('D', 'I'):
            # First gap operation sets the type
            if gap_type is None:
                gap_type = op
            # Mixed gap types means we've hit the end of contiguous gaps
            elif gap_type != op:
                # Reset since we found mixed types - not a clean suffix gap
                gap_length = 0
                break
            gap_length += int(length_str)
        else:
            # Match/mismatch operations end the gap region
            break
    
    if gap_length > 0:
        gap_in_query = (gap_type == 'D')  # D = deletion = gap in query
        return (gap_length, gap_in_query)
    
    return None

def _find_scoring_region(seq1_aligned, seq2_aligned, end_skip_distance):
    """
    Find the [start, end] region of the alignment where mismatches should be counted.
    
    Implements MycoBLAST "digital end trimming" by skipping the first/last end_skip_distance
    nucleotides (not alignment positions) from each sequence to avoid counting sequencing 
    artifacts near read ends.
    
    Special case: When end_skip_distance=0, only score positions where both sequences 
    have non-gap characters (no overhang scoring).
    
    IMPORTANT: This function counts NUCLEOTIDES (non-gap characters), not alignment positions.
    End trimming only activates when both sequences have >= end_skip_distance nucleotides
    available to skip from each end.
    
    Behavior:
    - end_skip_distance=0: Score only overlap region (both sequences have content)
    - Short sequences (< 2×end_skip_distance nucleotides): Returns full range [0, len-1]
    - Long sequences (≥ 2×end_skip_distance nucleotides): Returns trimmed range excluding ends
    - Gap characters ('-') are ignored when counting nucleotides
    
    Example:
        seq1_aligned = "AAAA-TCGX-TTTT"  # 12 nucleotides, 14 alignment positions
        seq2_aligned = "-AAAATCGA-TTTT"  # 12 nucleotides, 14 alignment positions
        end_skip_distance = 3
        
        Result: Skip first 3 and last 3 nucleotides from each sequence
        → scoring_start=4, scoring_end=9 (positions where middle nucleotides align)
    
    Args:
        seq1_aligned, seq2_aligned: Aligned sequences with gaps (must be same length)
        end_skip_distance: Number of nucleotides to skip from each end (typically 20)
        
    Returns:
        tuple: (scoring_start, scoring_end) - inclusive range of alignment positions to score
    """
    alignment_length = len(seq1_aligned)
    
    # Special case: end_skip_distance=0 means score only overlap region
    if end_skip_distance == 0:
        # Find first position where both sequences have content
        scoring_start = 0
        for pos in range(alignment_length):
            if seq1_aligned[pos] != '-' and seq2_aligned[pos] != '-':
                scoring_start = pos
                break
        
        # Find last position where both sequences have content
        scoring_end = alignment_length - 1
        for pos in range(alignment_length - 1, -1, -1):
            if seq1_aligned[pos] != '-' and seq2_aligned[pos] != '-':
                scoring_end = pos
                break
        
        return scoring_start, scoring_end
    
    # General case: skip end_skip_distance nucleotides from each end
    # Find scoring start: first position where both sequences have >= end_skip_distance bp
    seq1_count = seq2_count = 0
    scoring_start = 0
    for pos in range(alignment_length):
        if seq1_aligned[pos] != '-':
            seq1_count += 1
        if seq2_aligned[pos] != '-':
            seq2_count += 1
        if seq1_count >= end_skip_distance and seq2_count >= end_skip_distance:
            scoring_start = pos
            break
    
    # Find scoring end: last position where both sequences have >= end_skip_distance bp remaining
    seq1_count = seq2_count = 0
    scoring_end = alignment_length - 1
    for pos in range(alignment_length - 1, -1, -1):
        if seq1_aligned[pos] != '-':
            seq1_count += 1
        if seq2_aligned[pos] != '-':
            seq2_count += 1
        if seq1_count >= end_skip_distance and seq2_count >= end_skip_distance:
            scoring_end = pos
            break
    
    return scoring_start, scoring_end


def _extract_left_context(seq1_aligned, seq2_aligned, position, length):
    """
    Extract up to 'length' nucleotide characters from positions before 'position',
    working backwards, enforcing consensus between sequences for MSA support.

    Context extraction rules (MSA-compatible):
    - Both sequences agree (same non-gap char) → use it (valid consensus)
    - One has gap, other has char → use the character (unambiguous)
    - Both have gaps (dual-gap) → skip position entirely
    - Both have different non-gap chars → conflict, return None (no consensus)

    Args:
        seq1_aligned, seq2_aligned: Aligned sequences with gaps
        position: Starting position (exclusive, work backward from here)
        length: Number of nucleotide characters to extract

    Returns:
        str: Context string in left-to-right order (e.g., "AAA"), or None if:
            - Insufficient context available (fewer than 'length' chars)
            - Sequences disagree at context position (conflicting characters)

    Examples:
        >>> _extract_left_context("AGG-AC", "AG-GAC", 2, 1)
        'G'  # Both sequences agree at position 1
        >>> _extract_left_context("AGT-AC", "AX-GAC", 2, 1)
        None  # Sequences disagree at position 1 (G vs X)
    """
    context_chars = []
    pos = position - 1

    while len(context_chars) < length and pos >= 0:
        char1 = seq1_aligned[pos]
        char2 = seq2_aligned[pos]

        # Skip dual-gap positions (common in MSA-derived alignments)
        if char1 == '-' and char2 == '-':
            pos -= 1
            continue

        # Both sequences have content - they must AGREE for valid consensus
        if char1 != '-' and char2 != '-':
            if char1.upper() == char2.upper():
                context_chars.append(char1)
            else:
                # Disagreement = no clear consensus context
                # Cannot reliably detect homopolymer extension
                return None
        else:
            # One has gap, other has character - use the non-gap character
            context_chars.append(char1 if char1 != '-' else char2)

        pos -= 1

    # Check if we collected enough context characters
    if len(context_chars) < length:
        return None  # Insufficient context available

    # Reverse to get left-to-right order
    return ''.join(reversed(context_chars))


def _extract_right_context(seq1_aligned, seq2_aligned, position, length):
    """
    Extract up to 'length' nucleotide characters from positions after 'position',
    working forwards, enforcing consensus between sequences for MSA support.

    Context extraction rules (MSA-compatible):
    - Both sequences agree (same non-gap char) → use it (valid consensus)
    - One has gap, other has char → use the character (unambiguous)
    - Both have gaps (dual-gap) → skip position entirely
    - Both have different non-gap chars → conflict, return None (no consensus)

    Args:
        seq1_aligned, seq2_aligned: Aligned sequences with gaps
        position: Starting position (exclusive, work forward from here)
        length: Number of nucleotide characters to extract

    Returns:
        str: Context string in left-to-right order (e.g., "TTT"), or None if:
            - Insufficient context available (fewer than 'length' chars)
            - Sequences disagree at context position (conflicting characters)

    Examples:
        >>> _extract_right_context("AGA--TT", "AGAT-TT", 4, 2)
        'TT'  # Both sequences agree at positions 5-6
        >>> _extract_right_context("AGA--TC", "AGAT-TG", 4, 2)
        None  # Sequences disagree at position 6 (C vs G)
    """
    context_chars = []
    pos = position + 1
    max_pos = len(seq1_aligned)

    while len(context_chars) < length and pos < max_pos:
        char1 = seq1_aligned[pos]
        char2 = seq2_aligned[pos]

        # Skip dual-gap positions (common in MSA-derived alignments)
        if char1 == '-' and char2 == '-':
            pos += 1
            continue

        # Both sequences have content - they must AGREE for valid consensus
        if char1 != '-' and char2 != '-':
            if char1.upper() == char2.upper():
                context_chars.append(char1)
            else:
                # Disagreement = no clear consensus context
                # Cannot reliably detect homopolymer extension
                return None
        else:
            # One has gap, other has character - use the non-gap character
            context_chars.append(char1 if char1 != '-' else char2)

        pos += 1

    # Check if we collected enough context characters
    if len(context_chars) < length:
        return None  # Insufficient context available

    return ''.join(context_chars)


# =============================================================================
# Variant Range Algorithm Functions (v0.2.0)
# =============================================================================

def _extract_allele(seq_aligned, start, end):
    """
    Extract non-gap characters from aligned sequence within range [start, end].

    Args:
        seq_aligned: Aligned sequence with gaps
        start: Start position (inclusive)
        end: End position (inclusive)

    Returns:
        tuple: (allele_string, list_of_source_positions)
    """
    chars = []
    positions = []
    for i in range(start, end + 1):
        if seq_aligned[i] != '-':
            chars.append(seq_aligned[i])
            positions.append(i)
    return (''.join(chars), positions)


def _motif_matches(chunk, motif, handle_iupac):
    """
    Check if chunk matches motif using IUPAC equivalence.

    Args:
        chunk: String to check
        motif: Motif pattern to match against
        handle_iupac: Whether to use IUPAC intersection matching

    Returns:
        bool: True if chunk matches motif
    """
    if len(chunk) != len(motif):
        return False
    for c1, c2 in zip(chunk, motif):
        is_match, _ = _are_nucleotides_equivalent(c1, c2, handle_iupac)
        if not is_match:
            return False
    return True


def _analyze_allele(allele, left_context, right_context, max_motif_length, handle_iupac):
    """
    Analyze an allele to determine what portions are repeat extensions.

    Uses split scoring: portions matching context are extensions, remainder is core.
    Supports IUPAC equivalence for extension matching.

    Args:
        allele: The allele string to analyze
        left_context: Context from left of variant range (or None)
        right_context: Context from right of variant range (or None)
        max_motif_length: Maximum motif length to try
        handle_iupac: Whether to use IUPAC intersection matching

    Returns:
        AlleleAnalysis: Analysis result with extension counts and core content
    """
    if not allele:
        return AlleleAnalysis(0, 0, '', True)  # Empty allele = pure extension

    chars = list(allele)
    n = len(chars)
    left_consumed = 0
    right_consumed = 0

    # LEFT EXTENSION: Try different motif lengths (largest first)
    if left_context:
        for motif_len in range(min(max_motif_length, len(left_context)), 0, -1):
            motif = left_context[-motif_len:]  # Last motif_len chars of left context

            # Check for degenerate case (homopolymer disguised as longer motif)
            if len(set(motif.upper())) == 1:
                motif = motif[0]
                motif_len = 1

            # Count complete motif matches from left
            consumed = 0
            pos = 0
            while pos + motif_len <= n:
                chunk = ''.join(chars[pos:pos + motif_len])
                if _motif_matches(chunk, motif, handle_iupac):
                    consumed += motif_len
                    pos += motif_len
                else:
                    break

            if consumed > 0:
                left_consumed = consumed
                break

    # RIGHT EXTENSION: Try different motif lengths (largest first)
    if right_context:
        remaining_start = left_consumed
        remaining_chars = chars[remaining_start:]

        for motif_len in range(min(max_motif_length, len(right_context)), 0, -1):
            motif = right_context[:motif_len]  # First motif_len chars of right context

            # Check for degenerate case
            if len(set(motif.upper())) == 1:
                motif = motif[0]
                motif_len = 1

            # Count complete motif matches from right
            consumed = 0
            pos = len(remaining_chars)

            while pos >= motif_len:
                chunk = ''.join(remaining_chars[pos - motif_len:pos])
                if _motif_matches(chunk, motif, handle_iupac):
                    consumed += motif_len
                    pos -= motif_len
                else:
                    break

            if consumed > 0:
                right_consumed = consumed
                break

    # Core content is what remains
    core_start = left_consumed
    core_end = n - right_consumed
    core_content = ''.join(chars[core_start:core_end]) if core_end > core_start else ''

    return AlleleAnalysis(
        left_extension_count=left_consumed,
        right_extension_count=right_consumed,
        core_content=core_content,
        is_pure_extension=(len(core_content) == 0)
    )


def _find_variant_ranges(seq1_aligned, seq2_aligned, scoring_start, scoring_end, handle_iupac):
    """
    Find all variant ranges within the scoring region.

    A variant range is a maximal contiguous region where at least one position
    is NOT a non-gap match or dual-gap. Bounded by match positions on left/right.

    Args:
        seq1_aligned, seq2_aligned: Aligned sequences with gaps
        scoring_start, scoring_end: Scoring region boundaries (inclusive)
        handle_iupac: Whether to use IUPAC intersection matching

    Returns:
        List of tuples: [(start, end, left_bound_pos, right_bound_pos), ...]
        where start/end are inclusive positions of the variant range,
        and bound positions point to the bounding match positions (-1 if none).
    """
    variant_ranges = []
    i = scoring_start

    def is_match_position(pos):
        """Check if position is a match for variant range boundary detection.

        Only non-gap matches count as boundaries. Dual-gaps are NOT matches
        because they should be included within variant ranges, not split them.
        """
        c1, c2 = seq1_aligned[pos], seq2_aligned[pos]
        # Both must be non-gap and equivalent to be a match
        if c1 != '-' and c2 != '-':
            is_match, _ = _are_nucleotides_equivalent(c1, c2, handle_iupac)
            return is_match
        # Any gap (single-sided or dual) is not a match for boundary purposes
        return False

    while i <= scoring_end:
        if is_match_position(i):
            i += 1
            continue

        # Found start of variant range
        variant_start = i
        # Left bound is the position just before variant start (if within scoring region)
        left_bound_pos = i - 1 if i > scoring_start else -1

        # Scan to find end of variant range
        while i <= scoring_end:
            if is_match_position(i):
                break
            i += 1

        variant_end = i - 1
        # Right bound is the position just after variant end (if within scoring region)
        right_bound_pos = i if i <= scoring_end else -1

        variant_ranges.append((variant_start, variant_end, left_bound_pos, right_bound_pos))

    return variant_ranges


def _score_variant_range(allele1, analysis1, allele2, analysis2, adjustment_params):
    """
    Score a variant range given two allele analyses using Occam's razor.

    Scoring rules (when normalize_homopolymers=True):
    - Both pure extensions → 0 edits (homopolymer equivalent)
    - One pure extension, other has core → score the core as edits
    - Both have core → compare cores

    When normalize_homopolymers=False:
    - Extensions are treated as regular indels

    Args:
        allele1, allele2: The allele strings
        analysis1, analysis2: AlleleAnalysis for each allele
        adjustment_params: AdjustmentParams for scoring behavior

    Returns:
        dict: {
            'edits': int,
            'scored_positions': int,
            'both_pure_extension': bool
        }
    """
    # When homopolymer normalization is disabled, treat extensions as indels
    if not adjustment_params.normalize_homopolymers:
        # Total content from both alleles (treating as a regular indel region)
        total_len = max(len(allele1), len(allele2))
        if total_len == 0:
            return {'edits': 0, 'scored_positions': 0, 'both_pure_extension': False}

        if adjustment_params.normalize_indels:
            return {
                'edits': 1,
                'scored_positions': 1,
                'both_pure_extension': False
            }
        else:
            return {
                'edits': total_len,
                'scored_positions': total_len,
                'both_pure_extension': False
            }

    # Both pure extensions -> homopolymer equivalent
    if analysis1.is_pure_extension and analysis2.is_pure_extension:
        return {
            'edits': 0,
            'scored_positions': 0,
            'both_pure_extension': True
        }

    # One pure extension, other has core -> count core as edits
    if analysis1.is_pure_extension:
        core = analysis2.core_content
        if adjustment_params.normalize_indels:
            return {
                'edits': 1 if core else 0,
                'scored_positions': 1 if core else 0,
                'both_pure_extension': False
            }
        else:
            return {
                'edits': len(core),
                'scored_positions': len(core),
                'both_pure_extension': False
            }

    if analysis2.is_pure_extension:
        core = analysis1.core_content
        if adjustment_params.normalize_indels:
            return {
                'edits': 1 if core else 0,
                'scored_positions': 1 if core else 0,
                'both_pure_extension': False
            }
        else:
            return {
                'edits': len(core),
                'scored_positions': len(core),
                'both_pure_extension': False
            }

    # Both have core -> compare cores
    core1, core2 = analysis1.core_content, analysis2.core_content

    if core1 == core2:
        # Cores are identical - just extension differences
        return {
            'edits': 0,
            'scored_positions': len(core1),
            'both_pure_extension': False
        }

    # Cores differ - compute edit count
    min_len = min(len(core1), len(core2))
    max_len = max(len(core1), len(core2))

    # Count substitutions in overlapping region
    substitutions = sum(1 for i in range(min_len)
                        if not _are_nucleotides_equivalent(core1[i], core2[i],
                                                          adjustment_params.handle_iupac_overlap)[0])

    # Handle length difference as indel
    if max_len > min_len:
        if adjustment_params.normalize_indels:
            indel_edits = 1
            indel_scored = 1
        else:
            indel_edits = max_len - min_len
            indel_scored = max_len - min_len
    else:
        indel_edits = 0
        indel_scored = 0

    return {
        'edits': substitutions + indel_edits,
        'scored_positions': min_len + indel_scored,
        'both_pure_extension': False
    }


def _generate_variant_score_strings(seq1_aligned, seq2_aligned, start, end,
                                     analysis1, analysis2, allele1_positions, allele2_positions,
                                     scoring_format, adjustment_params):
    """
    Generate score_aligned strings for both sequences in a variant range.

    Each sequence gets its own visualization string. The visualization reflects
    how positions were scored:
    - Extension positions show extension marker (=)
    - Core positions show match (|) if cores match, mismatch ( ) if they differ
    - Gap positions mirror the other sequence's marker (borrowing for extensions)

    Args:
        seq1_aligned, seq2_aligned: Aligned sequences
        start, end: Variant range boundaries (inclusive)
        analysis1, analysis2: AlleleAnalysis for each allele
        allele1_positions, allele2_positions: Source positions for each allele's chars
        scoring_format: ScoringFormat for visualization
        adjustment_params: AdjustmentParams for scoring behavior

    Returns:
        tuple: (score_string_seq1, score_string_seq2)
    """
    # Build position classification sets for seq1
    seq1_left_ext_positions = set(allele1_positions[:analysis1.left_extension_count])
    seq1_right_ext_positions = set(allele1_positions[-analysis1.right_extension_count:]
                                   if analysis1.right_extension_count > 0 else [])
    seq1_ext_positions = seq1_left_ext_positions | seq1_right_ext_positions

    # Build position classification sets for seq2
    seq2_left_ext_positions = set(allele2_positions[:analysis2.left_extension_count])
    seq2_right_ext_positions = set(allele2_positions[-analysis2.right_extension_count:]
                                   if analysis2.right_extension_count > 0 else [])
    seq2_ext_positions = seq2_left_ext_positions | seq2_right_ext_positions

    # Calculate core positions for each allele
    core1_start = analysis1.left_extension_count
    core1_end = len(allele1_positions) - analysis1.right_extension_count
    seq1_core_positions = set(allele1_positions[core1_start:core1_end]) if core1_end > core1_start else set()

    core2_start = analysis2.left_extension_count
    core2_end = len(allele2_positions) - analysis2.right_extension_count
    seq2_core_positions = set(allele2_positions[core2_start:core2_end]) if core2_end > core2_start else set()

    # Determine if cores match (for visualization - matched cores show as |)
    cores_match = analysis1.core_content == analysis2.core_content

    # Choose extension marker based on normalization settings
    ext_marker = (scoring_format.homopolymer_extension if adjustment_params.normalize_homopolymers
                  else scoring_format.indel_extension)

    # Generate visualization for each sequence
    score_chars_seq1 = []
    score_chars_seq2 = []
    seen_core_start = False  # Track if we've seen the first core position (for indel normalization)

    for pos in range(start, end + 1):
        char1 = seq1_aligned[pos]
        char2 = seq2_aligned[pos]

        # Case 1: Dual-gap - both show match
        if char1 == '-' and char2 == '-':
            score_chars_seq1.append(scoring_format.match)
            score_chars_seq2.append(scoring_format.match)
            continue

        # Case 2: Both have content
        if char1 != '-' and char2 != '-':
            # When homopolymer normalization is disabled, treat all content as core (no extensions)
            is_ext1 = adjustment_params.normalize_homopolymers and pos in seq1_ext_positions
            is_ext2 = adjustment_params.normalize_homopolymers and pos in seq2_ext_positions
            is_core1 = pos in seq1_core_positions or not adjustment_params.normalize_homopolymers
            is_core2 = pos in seq2_core_positions or not adjustment_params.normalize_homopolymers

            # If one is extension and other is core, show based on whether cores match
            if (is_ext1 and is_core2) or (is_ext2 and is_core1):
                if cores_match:
                    # Cores match - no mismatch counted
                    score_chars_seq1.append(scoring_format.match)
                    score_chars_seq2.append(scoring_format.match)
                else:
                    # Cores differ - mismatch counted
                    score_chars_seq1.append(scoring_format.substitution)
                    score_chars_seq2.append(scoring_format.substitution)
            elif is_ext1 and is_ext2:
                # Both extensions - show extension marker
                score_chars_seq1.append(ext_marker)
                score_chars_seq2.append(ext_marker)
            elif is_core1 and is_core2:
                # Both core - show match or mismatch based on core comparison
                marker = scoring_format.match if cores_match else scoring_format.substitution
                score_chars_seq1.append(marker)
                score_chars_seq2.append(marker)
            else:
                # Fallback (shouldn't happen normally)
                score_chars_seq1.append(scoring_format.substitution)
                score_chars_seq2.append(scoring_format.substitution)
            continue

        # Case 3: seq1 has gap, seq2 has content
        if char1 == '-':
            # When homopolymer normalization is disabled, treat all content as indels (no extensions)
            is_ext2 = adjustment_params.normalize_homopolymers and pos in seq2_ext_positions
            is_core2 = pos in seq2_core_positions or not adjustment_params.normalize_homopolymers

            if is_ext2:
                # seq2 is extension - seq1 borrows extension marker
                score_chars_seq1.append(ext_marker)
                score_chars_seq2.append(ext_marker)
            elif is_core2:
                # seq2 is core - show indel markers (start for first, extension for rest)
                if adjustment_params.normalize_indels and seen_core_start:
                    # Subsequent core position - show as indel extension
                    score_chars_seq1.append(scoring_format.indel_extension)
                    score_chars_seq2.append(scoring_format.indel_extension)
                else:
                    # First core position - show as indel start
                    score_chars_seq1.append(scoring_format.indel_start)
                    score_chars_seq2.append(scoring_format.indel_start)
                    seen_core_start = True
            else:
                # Fallback (shouldn't happen)
                score_chars_seq1.append(scoring_format.substitution)
                score_chars_seq2.append(scoring_format.substitution)
            continue

        # Case 4: seq1 has content, seq2 has gap
        if char2 == '-':
            # When homopolymer normalization is disabled, treat all content as indels (no extensions)
            is_ext1 = adjustment_params.normalize_homopolymers and pos in seq1_ext_positions
            is_core1 = pos in seq1_core_positions or not adjustment_params.normalize_homopolymers

            if is_ext1:
                # seq1 is extension - seq2 borrows extension marker
                score_chars_seq1.append(ext_marker)
                score_chars_seq2.append(ext_marker)
            elif is_core1:
                # seq1 is core - show indel markers (start for first, extension for rest)
                if adjustment_params.normalize_indels and seen_core_start:
                    # Subsequent core position - show as indel extension
                    score_chars_seq1.append(scoring_format.indel_extension)
                    score_chars_seq2.append(scoring_format.indel_extension)
                else:
                    # First core position - show as indel start
                    score_chars_seq1.append(scoring_format.indel_start)
                    score_chars_seq2.append(scoring_format.indel_start)
                    seen_core_start = True
            else:
                # Fallback (shouldn't happen)
                score_chars_seq1.append(scoring_format.substitution)
                score_chars_seq2.append(scoring_format.substitution)
            continue

    return ''.join(score_chars_seq1), ''.join(score_chars_seq2)


def score_alignment(seq1_aligned, seq2_aligned, adjustment_params=None, scoring_format=None):
    """
    Score alignment and count edits with configurable MycoBLAST-style adjustments.

    Applies various preprocessing adjustments based on adjustment_params:
    - End trimming: Skip mismatches within end_skip_distance bp from either end (set 0 to disable)
    - Homopolymer adjustment: Ignore differences in homopolymer run lengths
    - IUPAC handling: Allow different ambiguity codes to match via intersection
    - Indel normalization: Count contiguous indels as single evolutionary events

    Args:
        seq1_aligned (str): First sequence with gaps ('-') inserted
        seq2_aligned (str): Second sequence with gaps ('-') inserted
        adjustment_params (AdjustmentParams, optional): Parameters controlling which adjustments to apply.
                                                       Defaults to DEFAULT_ADJUSTMENT_PARAMS.
        scoring_format (ScoringFormat, optional): Format codes for alignment visualization.
                                                 Defaults to DEFAULT_SCORING_FORMAT.

    Returns:
        AlignmentResult: Dataclass containing:
            - identity (float): Identity score based on adjustment parameters
            - mismatches (int): Number of mismatches/edits counted
            - scored_positions (int): Number of positions used for identity calculation
            - seq1_coverage (float): Fraction of seq1 used in scoring region
            - seq2_coverage (float): Fraction of seq2 used in scoring region
            - seq1_aligned (str): Input seq1_aligned (passed through for consistency)
            - seq2_aligned (str): Input seq2_aligned (passed through for consistency)
            - score_aligned (str): Scoring codes string for visualization
    """
    # Use default parameters if none provided
    if adjustment_params is None:
        adjustment_params = DEFAULT_ADJUSTMENT_PARAMS
    if scoring_format is None:
        scoring_format = DEFAULT_SCORING_FORMAT
    
    # Input validation: aligned sequences must have the same length
    if len(seq1_aligned) != len(seq2_aligned):
        raise ValueError(f"Aligned sequences must have same length: seq1={len(seq1_aligned)}, seq2={len(seq2_aligned)}")
    
    edits = 0
    scored_positions = 0  # Single counter for identity denominator
    
    # Total length for iteration
    total_alignment_length = len(seq1_aligned)


    # Calculate coverage independently (simple counting within alignment region)
    seq1_coverage_positions = 0
    seq2_coverage_positions = 0
    seq1_total_length = 0
    seq2_total_length = 0
    
    # Find alignment bounds inline
    alignment_start = 0
    alignment_end = total_alignment_length - 1
    
    # Find first position where both sequences have content
    for pos in range(total_alignment_length):
        if seq1_aligned[pos] != '-' and seq2_aligned[pos] != '-':
            alignment_start = pos
            break
    
    # Find last position where both sequences have content
    for pos in range(total_alignment_length - 1, -1, -1):
        if seq1_aligned[pos] != '-' and seq2_aligned[pos] != '-':
            alignment_end = pos
            break
    
    for pos in range(total_alignment_length):
        # Count total sequence lengths (all non-gap characters)
        if seq1_aligned[pos] != '-':
            seq1_total_length += 1
        if seq2_aligned[pos] != '-':
            seq2_total_length += 1
            
        # Count coverage positions (non-gap characters within alignment region)
        if alignment_start <= pos <= alignment_end:
            if seq1_aligned[pos] != '-':
                seq1_coverage_positions += 1
            if seq2_aligned[pos] != '-':
                seq2_coverage_positions += 1

    # Always build alignment scoring codes for visualization
    score_aligned_seq1 = []
    score_aligned_seq2 = []

    # Find the scoring region boundaries (end trimming controlled by end_skip_distance)
    scoring_start, scoring_end = _find_scoring_region(
        seq1_aligned, seq2_aligned, adjustment_params.end_skip_distance
    )

    # =========================================================================
    # Variant Range Algorithm (v0.2.0)
    # =========================================================================
    # Process alignment by:
    # 1. Add end-trimmed markers for positions outside scoring region
    # 2. Find variant ranges (contiguous non-match regions)
    # 3. Process match positions and variant ranges

    # Add end-trimmed markers for positions before scoring region
    for i in range(scoring_start):
        score_aligned_seq1.append(scoring_format.end_trimmed)
        score_aligned_seq2.append(scoring_format.end_trimmed)

    # Find all variant ranges within the scoring region
    variant_ranges = _find_variant_ranges(
        seq1_aligned, seq2_aligned, scoring_start, scoring_end,
        adjustment_params.handle_iupac_overlap
    )

    # Build a set of positions that are in variant ranges for quick lookup
    variant_positions = set()
    for vr_start, vr_end, _, _ in variant_ranges:
        for pos in range(vr_start, vr_end + 1):
            variant_positions.add(pos)

    # Process each position in the scoring region
    pos = scoring_start
    vr_index = 0  # Current variant range index

    while pos <= scoring_end:
        # Check if this position starts a variant range
        if vr_index < len(variant_ranges) and pos == variant_ranges[vr_index][0]:
            vr_start, vr_end, left_bound, right_bound = variant_ranges[vr_index]

            # Extract alleles from both sequences
            allele1, allele1_positions = _extract_allele(seq1_aligned, vr_start, vr_end)
            allele2, allele2_positions = _extract_allele(seq2_aligned, vr_start, vr_end)

            # Get context for homopolymer detection
            # Try different context lengths (largest to smallest) until we get valid context
            max_ctx_len = adjustment_params.max_repeat_motif_length
            left_context = None
            right_context = None

            if left_bound >= 0:
                for ctx_len in range(max_ctx_len, 0, -1):
                    left_context = _extract_left_context(
                        seq1_aligned, seq2_aligned, vr_start, ctx_len)
                    if left_context is not None:
                        break

            if right_bound >= 0:
                for ctx_len in range(max_ctx_len, 0, -1):
                    right_context = _extract_right_context(
                        seq1_aligned, seq2_aligned, vr_end, ctx_len)
                    if right_context is not None:
                        break

            # Analyze alleles
            analysis1 = _analyze_allele(
                allele1, left_context, right_context,
                adjustment_params.max_repeat_motif_length,
                adjustment_params.handle_iupac_overlap
            )
            analysis2 = _analyze_allele(
                allele2, left_context, right_context,
                adjustment_params.max_repeat_motif_length,
                adjustment_params.handle_iupac_overlap
            )

            # Score the variant range using Occam's razor
            vr_score = _score_variant_range(
                allele1, analysis1, allele2, analysis2, adjustment_params
            )

            edits += vr_score['edits']
            scored_positions += vr_score['scored_positions']

            # Generate score strings for this variant range (one for each sequence)
            vr_score_seq1, vr_score_seq2 = _generate_variant_score_strings(
                seq1_aligned, seq2_aligned, vr_start, vr_end,
                analysis1, analysis2, allele1_positions, allele2_positions,
                scoring_format, adjustment_params
            )
            score_aligned_seq1.append(vr_score_seq1)
            score_aligned_seq2.append(vr_score_seq2)

            # Move past this variant range
            pos = vr_end + 1
            vr_index += 1

        else:
            # This is a match position (both non-gap and equivalent, or dual-gap)
            char1, char2 = seq1_aligned[pos], seq2_aligned[pos]

            # Check for dual-gap (treated as match)
            if char1 == '-' and char2 == '-':
                scored_positions += 1
                score_aligned_seq1.append(scoring_format.match)
                score_aligned_seq2.append(scoring_format.match)
                pos += 1
                continue

            is_match, is_ambiguous = _are_nucleotides_equivalent(
                char1, char2, adjustment_params.handle_iupac_overlap)

            scored_positions += 1
            if is_ambiguous:
                score_aligned_seq1.append(scoring_format.ambiguous_match)
                score_aligned_seq2.append(scoring_format.ambiguous_match)
            else:
                score_aligned_seq1.append(scoring_format.match)
                score_aligned_seq2.append(scoring_format.match)
            pos += 1

    # Add end-trimmed markers for positions after scoring region
    for i in range(scoring_end + 1, total_alignment_length):
        score_aligned_seq1.append(scoring_format.end_trimmed)
        score_aligned_seq2.append(scoring_format.end_trimmed)

    # Calculate coverage as fraction of sequence used in alignment region
    seq1_coverage = seq1_coverage_positions / seq1_total_length if seq1_total_length > 0 else 0.0
    seq2_coverage = seq2_coverage_positions / seq2_total_length if seq2_total_length > 0 else 0.0

    # Create scoring codes strings for visualization
    score_aligned_str = ''.join(score_aligned_seq1)
    score_aligned_str_seq2 = ''.join(score_aligned_seq2)

    # Calculate identity metric: identity = 1 - (edits / scored_positions)
    identity = 1.0 - (edits / scored_positions) if scored_positions > 0 else 1.0

    # Return AlignmentResult with calculated identity value
    return AlignmentResult(
        identity=identity,
        mismatches=edits,
        scored_positions=scored_positions,
        seq1_coverage=seq1_coverage,
        seq2_coverage=seq2_coverage,
        seq1_aligned=seq1_aligned,
        seq2_aligned=seq2_aligned,
        score_aligned=score_aligned_str,
        score_aligned_seq2=score_aligned_str_seq2
    )


def align_edlib_bidirectional(seq1, seq2):
    """
    Multi-stage alignment optimization using CIGAR-based suffix detection.

    Process:
    1. Reverse complement both sequences
    2. Global alignment with task=locations, parse CIGAR for suffix gaps
    3. Trim sequences based on gap info
    4. Reverse complement back to forward orientation
    5. Global alignment with task=path for final result
    6. Parse CIGAR for suffix gaps and trim final alignment

    Args:
        seq1, seq2: Original DNA sequences

    Returns:
        dict: final_alignment with 'aligned_seq1' and 'aligned_seq2' keys
        Or None if alignment fails
    """

    # Safety check: ensure sequences are non-empty
    if len(seq1) == 0 or len(seq2) == 0:
        return None  # Sentinel value for failed alignment

    current_seq1, current_seq2 = seq1, seq2

    # Track trimming information with local variables
    seq1_prefix_trimmed = 0
    seq1_suffix_trimmed = 0
    seq2_prefix_trimmed = 0
    seq2_suffix_trimmed = 0

    # Step 1: Reverse complement both sequences
    rc_seq1 = _reverse_complement(current_seq1)
    rc_seq2 = _reverse_complement(current_seq2)

    # Step 2: RC alignment with task=path for CIGAR-based gap detection
    result = edlib.align(rc_seq1, rc_seq2, mode="HW", task="path")

    # Step 3: Check for alignment failure
    if result['editDistance'] == -1:
        return None  # Sentinel value for failed alignment

    # Step 4: Parse CIGAR for suffix gaps (original prefix gaps)
    cigar_string = result.get('cigar')
    if cigar_string:
        gap_info = _parse_suffix_gap_from_cigar(cigar_string)

        if gap_info:
            gap_length, gap_in_query = gap_info

            # Trim the RC sequences
            if gap_in_query:
                # Query (rc_seq1) has gaps, trim rc_seq2
                rc_seq2 = rc_seq2[:-gap_length] if gap_length > 0 else rc_seq2
                seq2_prefix_trimmed = gap_length  # suffix in RC = prefix in forward
            else:
                # Target (rc_seq2) has gaps, trim rc_seq1
                rc_seq1 = rc_seq1[:-gap_length] if gap_length > 0 else rc_seq1
                seq1_prefix_trimmed = gap_length  # suffix in RC = prefix in forward

    # Step 5: Reverse complement back to forward orientation
    current_seq1 = _reverse_complement(rc_seq1)
    current_seq2 = _reverse_complement(rc_seq2)

    # Step 6: Final forward alignment with task=path
    result = edlib.align(current_seq1, current_seq2, mode="HW", task="path")

    if result['editDistance'] == -1:
        return None  # Sentinel value for failed alignment

    # Step 8: Get nice alignment 
    alignment = edlib.getNiceAlignment(result, current_seq1, current_seq2)
    seq1_aligned = alignment['query_aligned']
    seq2_aligned = alignment['target_aligned']

    # Step 9: Re-attach removed prefix and suffix regions with gap padding
    # Only one sequence can have prefix trimmed, only one can have suffix trimmed

    # Handle prefix: one sequence has actual prefix, other gets gap padding
    if seq1_prefix_trimmed > 0:
        seq1_prefix_part = seq1[:seq1_prefix_trimmed]
        seq2_prefix_part = '-' * seq1_prefix_trimmed
    elif seq2_prefix_trimmed > 0:
        seq1_prefix_part = '-' * seq2_prefix_trimmed
        seq2_prefix_part = seq2[:seq2_prefix_trimmed]
    else:
        seq1_prefix_part = ""
        seq2_prefix_part = ""

    # Re-attach prefix and suffix to alignment
    seq1_aligned = seq1_prefix_part + seq1_aligned
    seq2_aligned = seq2_prefix_part + seq2_aligned

    # Return final alignment with re-attached sequences
    final_alignment = {
        'aligned_seq1': seq1_aligned,
        'aligned_seq2': seq2_aligned
    }

    return final_alignment


def align_and_score(seq1, seq2, adjustment_params=None, scoring_format=None):
    """
    Calculate adjusted and full identity between two DNA sequences.

    Implements the MycoBLAST preprocessing approach with configurable adjustments:
    - Homopolymer length normalization: ignore differences in homopolymer run lengths
    - IUPAC ambiguity code handling: allow different ambiguity codes to match via intersection
    - End trimming: skip mismatches in end regions (set end_skip_distance=0 to disable)
    - Indel normalization: count contiguous indels as single evolutionary events

    This is particularly useful for mycological DNA barcoding where technical artifacts
    can obscure true phylogenetic signal in sequence-based identifications.

    Args:
        seq1 (str): First DNA sequence
        seq2 (str): Second DNA sequence
        adjustment_params (AdjustmentParams, optional): Parameters controlling which adjustments to apply.
                                                       Defaults to DEFAULT_ADJUSTMENT_PARAMS.
        scoring_format (ScoringFormat, optional): Format codes for alignment visualization.
                                                 Defaults to DEFAULT_SCORING_FORMAT.

    Returns:
        AlignmentResult: Dataclass containing:
            - identity (float): Identity score based on adjustment parameters
            - mismatches (int): Number of mismatches/edits counted
            - scored_positions (int): Number of positions used for identity calculation
            - seq1_coverage (float): Fraction of seq1 used in scoring region
            - seq2_coverage (float): Fraction of seq2 used in scoring region
            - seq1_aligned (str): First aligned sequence with gaps
            - seq2_aligned (str): Second aligned sequence with gaps
            - score_aligned (str): Scoring codes string for visualization

    Example:
        >>> result = align_and_score("AAATTTGGG","AAAATTTGGG")
        >>> print(f"Identity: {result.identity:.3f}")
        >>> print(f"Coverage: {result.seq1_coverage:.3f}")
    """

    # Use default parameters if none provided
    if adjustment_params is None:
        adjustment_params = DEFAULT_ADJUSTMENT_PARAMS
    if scoring_format is None:
        scoring_format = DEFAULT_SCORING_FORMAT
    
    # Safety check: ensure sequences are non-empty
    if len(seq1) == 0 or len(seq2) == 0:
        return AlignmentResult(
            identity=0.0,
            mismatches=0,
            scored_positions=0,
            seq1_coverage=0.0,
            seq2_coverage=0.0,
            seq1_aligned='',
            seq2_aligned='',
            score_aligned=''
        )

    # Perform multi-stage bidirectional alignment with suffix trimming
    align_result = align_edlib_bidirectional(seq1, seq2)

    # Check for alignment failure (sentinel value)
    if align_result is None:
        # Alignment failed - return zero identity
        return AlignmentResult(
            identity=0.0,
            mismatches=-1,
            scored_positions=0,
            seq1_coverage=0.0,
            seq2_coverage=0.0,
            seq1_aligned='',
            seq2_aligned='',
            score_aligned=''
        )

    # Get successful alignment result
    final_alignment = align_result

    # Use the alignment returned from optimization
    seq1_aligned = final_alignment['aligned_seq1']
    seq2_aligned = final_alignment['aligned_seq2']

    # Score the alignment and calculate identity metrics
    # score_alignment now calculates everything including identity values
    return score_alignment(
        seq1_aligned, seq2_aligned, adjustment_params, scoring_format
    )

