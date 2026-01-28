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
        score_aligned: Scoring visualization string
    """
    identity: float
    mismatches: int
    scored_positions: int
    seq1_coverage: float
    seq2_coverage: float
    seq1_aligned: str
    seq2_aligned: str
    score_aligned: str
    


@dataclass(frozen=True)
class ScoringFormat:
    """Format codes for alignment scoring visualization."""
    match: str = '|'                    # Exact match (A=A, C=C, G=G, T=T)
    ambiguous_match: str = '='          # Ambiguous nucleotide match (different IUPAC codes with intersection)
    substitution: str = ' '             # Nucleotide substitution
    indel_start: str = ' '              # First position of indel (scored)
    indel_extension: str = '-'          # Indel positions skipped due to normalization
    homopolymer_extension: str = '='    # Homopolymer length difference
    end_trimmed: str = '.'              # Position outside scoring region (end trimmed)
    dual_gap: str = '.'                 # Both sequences have gap (MSA artifact, not scored)
    
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


@dataclass(frozen=True)
class VariantRangeInfo:
    """Information about a single variant range within an alignment.

    Used by the gap adjustment algorithm to track variant range boundaries
    and their allele analyses for reconstructing adjusted alignments.

    Attributes:
        start: Start position of variant range (inclusive)
        end: End position of variant range (inclusive)
        left_bound_pos: Position of left bounding match (-1 if none)
        right_bound_pos: Position of right bounding match (-1 if none)
        allele1: Non-gap content from seq1 in this range
        allele2: Non-gap content from seq2 in this range
        allele1_positions: Source positions for allele1 characters
        allele2_positions: Source positions for allele2 characters
        analysis1: AlleleAnalysis for allele1
        analysis2: AlleleAnalysis for allele2
        score_result: Scoring result dict from _score_variant_range
    """
    start: int
    end: int
    left_bound_pos: int
    right_bound_pos: int
    allele1: str
    allele2: str
    allele1_positions: tuple  # Use tuple for frozen dataclass (immutable)
    allele2_positions: tuple
    analysis1: AlleleAnalysis
    analysis2: AlleleAnalysis
    score_result: dict


@dataclass(frozen=True)
class AlignmentAnalysis:
    """Complete analysis of an alignment for gap adjustment.

    Contains all computed metrics and variant range information needed
    to reconstruct an adjusted alignment where gaps match the scoring analysis.

    Attributes:
        identity: Calculated identity score (0.0-1.0)
        mismatches: Number of edits/mismatches counted
        scored_positions: Number of positions used for identity calculation
        seq1_coverage: Fraction of seq1 used in alignment region
        seq2_coverage: Fraction of seq2 used in alignment region
        scoring_start: Start of scoring region (inclusive)
        scoring_end: End of scoring region (inclusive)
        variant_ranges: List of VariantRangeInfo for each variant range
    """
    identity: float
    mismatches: int
    scored_positions: int
    seq1_coverage: float
    seq2_coverage: float
    scoring_start: int
    scoring_end: int
    variant_ranges: tuple  # Tuple of VariantRangeInfo for frozen dataclass


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


# Fast lookup sets for common cases (avoid repeated set creation)
_STANDARD_NUCS = {'A', 'T', 'C', 'G', 'a', 't', 'c', 'g'}
_STANDARD_NUCS_AND_GAP = {'A', 'T', 'C', 'G', 'a', 't', 'c', 'g', '-'}


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
    # Fast path: exact match (case-sensitive) - most common case
    if nuc1 == nuc2:
        return (True, nuc1 not in _STANDARD_NUCS_AND_GAP)

    # Fast path: case-insensitive match for standard nucleotides
    if nuc1 in _STANDARD_NUCS and nuc2 in _STANDARD_NUCS:
        # Both are standard nucleotides - check case-insensitive match
        if nuc1.upper() == nuc2.upper():
            return (True, False)
        else:
            return (False, False)

    # Slow path: handle IUPAC codes - convert to uppercase
    nuc1_upper = nuc1.upper()
    nuc2_upper = nuc2.upper()

    # Check for exact match after uppercasing
    if nuc1_upper == nuc2_upper:
        is_standard_or_gap = nuc1_upper in {'A', 'T', 'C', 'G', '-'}
        return (True, not is_standard_or_gap)

    # Get possible nucleotides for each code
    possible1 = IUPAC_CODES.get(nuc1_upper, {nuc1_upper})
    possible2 = IUPAC_CODES.get(nuc2_upper, {nuc2_upper})

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
    Extract up to 'length' nucleotide characters from match positions before 'position'.

    Context is only valid from positions where both sequences agree (match positions).
    Since seq1 == seq2 at match positions, we read from seq1 and verify agreement.

    Rules:
    - Match (seq1 == seq2, non-gap) → use the character
    - Dual-gap (both '-') → skip (not a match, but not a conflict)
    - Any other case (gap or mismatch) → hit another variant range, return None

    Args:
        seq1_aligned, seq2_aligned: Aligned sequences with gaps
        position: Starting position (exclusive, work backward from here)
        length: Number of nucleotide characters to extract

    Returns:
        str: Context string in left-to-right order (e.g., "AAA"), or None if:
            - Insufficient context available (fewer than 'length' chars)
            - Hit another variant range (non-match position)

    Examples:
        >>> _extract_left_context("AGG-AC", "AG-GAC", 2, 1)
        'G'  # Position 1 is a match
        >>> _extract_left_context("AGT-AC", "AX-GAC", 2, 1)
        None  # Position 1 is not a match (G vs X)
    """
    context_chars = []
    pos = position - 1
    collected = 0

    while collected < length and pos >= 0:
        char1 = seq1_aligned[pos]
        char2 = seq2_aligned[pos]

        # Dual-gap: skip (not a match position, but doesn't end context search)
        if char1 == '-' and char2 == '-':
            pos -= 1
            continue

        # Check for match: must be non-gap and equal (case-insensitive)
        if char1 != '-' and char2 != '-' and (char1 == char2 or char1.upper() == char2.upper()):
            context_chars.append(char1)
            collected += 1
            pos -= 1
            continue

        # Any other case: we've hit a variant range (gap or mismatch)
        return None

    if collected < length:
        return None

    return ''.join(reversed(context_chars))


def _extract_right_context(seq1_aligned, seq2_aligned, position, length):
    """
    Extract up to 'length' nucleotide characters from match positions after 'position'.

    Context is only valid from positions where both sequences agree (match positions).
    Since seq1 == seq2 at match positions, we read from seq1 and verify agreement.

    Rules:
    - Match (seq1 == seq2, non-gap) → use the character
    - Dual-gap (both '-') → skip (not a match, but not a conflict)
    - Any other case (gap or mismatch) → hit another variant range, return None

    Args:
        seq1_aligned, seq2_aligned: Aligned sequences with gaps
        position: Starting position (exclusive, work forward from here)
        length: Number of nucleotide characters to extract

    Returns:
        str: Context string in left-to-right order (e.g., "TTT"), or None if:
            - Insufficient context available (fewer than 'length' chars)
            - Hit another variant range (non-match position)

    Examples:
        >>> _extract_right_context("AGA--TT", "AGAT-TT", 4, 2)
        'TT'  # Positions 5-6 are matches
        >>> _extract_right_context("AGA--TC", "AGAT-TG", 4, 2)
        None  # Position 6 is not a match (C vs G)
    """
    context_chars = []
    pos = position + 1
    max_pos = len(seq1_aligned)
    collected = 0

    while collected < length and pos < max_pos:
        char1 = seq1_aligned[pos]
        char2 = seq2_aligned[pos]

        # Dual-gap: skip (not a match position, but doesn't end context search)
        if char1 == '-' and char2 == '-':
            pos += 1
            continue

        # Check for match: must be non-gap and equal (case-insensitive)
        if char1 != '-' and char2 != '-' and (char1 == char2 or char1.upper() == char2.upper()):
            context_chars.append(char1)
            collected += 1
            pos += 1
            continue

        # Any other case: we've hit a variant range (gap or mismatch)
        return None

    if collected < length:
        return None

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

    n = len(allele)
    left_consumed = 0
    right_consumed = 0

    # LEFT EXTENSION: Try different motif lengths (largest first)
    if left_context:
        for motif_len in range(min(max_motif_length, len(left_context)), 0, -1):
            motif = left_context[-motif_len:]  # Last motif_len chars of left context

            # Check for degenerate case (homopolymer disguised as longer motif)
            # Fast check: compare first and last char instead of building a set
            if motif_len > 1 and motif[0].upper() == motif[-1].upper():
                # Check all chars are same
                first_upper = motif[0].upper()
                if all(c.upper() == first_upper for c in motif):
                    motif = motif[0]
                    motif_len = 1

            # Count complete motif matches from left (use string slicing, not list)
            consumed = 0
            pos = 0
            while pos + motif_len <= n:
                chunk = allele[pos:pos + motif_len]
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
        remaining_len = n - remaining_start

        for motif_len in range(min(max_motif_length, len(right_context)), 0, -1):
            motif = right_context[:motif_len]  # First motif_len chars of right context

            # Check for degenerate case (homopolymer disguised as longer motif)
            # Fast check: compare first and last char instead of building a set
            if motif_len > 1 and motif[0].upper() == motif[-1].upper():
                # Check all chars are same
                first_upper = motif[0].upper()
                if all(c.upper() == first_upper for c in motif):
                    motif = motif[0]
                    motif_len = 1

            # Count complete motif matches from right (use string slicing, not list)
            consumed = 0
            pos = remaining_len

            while pos >= motif_len:
                # Slice from the remaining portion of allele
                chunk_start = remaining_start + pos - motif_len
                chunk_end = remaining_start + pos
                chunk = allele[chunk_start:chunk_end]
                if _motif_matches(chunk, motif, handle_iupac):
                    consumed += motif_len
                    pos -= motif_len
                else:
                    break

            if consumed > 0:
                right_consumed = consumed
                break

    # Core content is what remains (string slicing, no list conversion)
    core_start = left_consumed
    core_end = n - right_consumed
    core_content = allele[core_start:core_end] if core_end > core_start else ''

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

    # Inline match checking for performance - avoid function call overhead
    while i <= scoring_end:
        c1, c2 = seq1_aligned[i], seq2_aligned[i]
        # Fast path: exact match of non-gap characters (most common case)
        if c1 == c2 and c1 != '-':
            i += 1
            continue
        # Check for gap - any gap means not a match
        if c1 == '-' or c2 == '-':
            pass  # Not a match, fall through to variant range handling
        # Slow path: check IUPAC equivalence for different characters
        elif _are_nucleotides_equivalent(c1, c2, handle_iupac)[0]:
            i += 1
            continue

        # Found start of variant range
        variant_start = i
        # Left bound is the position just before variant start (if within scoring region)
        left_bound_pos = i - 1 if i > scoring_start else -1

        # Scan to find end of variant range (inline match check for performance)
        i += 1
        while i <= scoring_end:
            c1, c2 = seq1_aligned[i], seq2_aligned[i]
            # Fast path: exact match of non-gap characters
            if c1 == c2 and c1 != '-':
                break
            # Gap means not a match - continue in variant range
            if c1 == '-' or c2 == '-':
                i += 1
                continue
            # Slow path: check IUPAC equivalence
            if _are_nucleotides_equivalent(c1, c2, handle_iupac)[0]:
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


def _generate_variant_score_string(seq1_aligned, seq2_aligned, start, end,
                                    analysis1, analysis2, allele1_positions, allele2_positions,
                                    scoring_format, adjustment_params):
    """
    Generate score_aligned string for a variant range.

    The visualization reflects how positions were scored:
    - Extension positions show extension marker (=)
    - Core positions show match (|) if cores match, mismatch ( ) if they differ
    - Gap positions show extension marker if the other sequence is an extension

    Args:
        seq1_aligned, seq2_aligned: Aligned sequences
        start, end: Variant range boundaries (inclusive)
        analysis1, analysis2: AlleleAnalysis for each allele
        allele1_positions, allele2_positions: Source positions for each allele's chars
        scoring_format: ScoringFormat for visualization
        adjustment_params: AdjustmentParams for scoring behavior

    Returns:
        str: Score visualization string for this variant range
    """
    # Build extension position sets for each sequence
    # Note: Extension and core positions form a partition of all content positions,
    # so is_core = not is_ext (no need to build separate core position sets)
    left1, right1 = analysis1.left_extension_count, analysis1.right_extension_count
    left2, right2 = analysis2.left_extension_count, analysis2.right_extension_count

    seq1_ext_positions = (set(allele1_positions[:left1]) |
                          set(allele1_positions[-right1:] if right1 > 0 else []))
    seq2_ext_positions = (set(allele2_positions[:left2]) |
                          set(allele2_positions[-right2:] if right2 > 0 else []))

    # Determine if cores match (for visualization - matched cores show as |)
    cores_match = analysis1.core_content == analysis2.core_content

    # Choose extension marker based on normalization settings
    ext_marker = (scoring_format.homopolymer_extension if adjustment_params.normalize_homopolymers
                  else scoring_format.indel_extension)

    # Generate visualization
    score_chars = []
    seen_core_start = False  # Track if we've seen the first core position (for indel normalization)

    for pos in range(start, end + 1):
        char1 = seq1_aligned[pos]
        char2 = seq2_aligned[pos]

        # Case 1: Dual-gap - both sequences have gap (MSA artifact, not scored)
        if char1 == '-' and char2 == '-':
            score_chars.append(scoring_format.dual_gap)
            continue

        # Case 2: Both have content
        if char1 != '-' and char2 != '-':
            # When homopolymer normalization is disabled, is_ext is always False (all content is core)
            is_ext1 = adjustment_params.normalize_homopolymers and pos in seq1_ext_positions
            is_ext2 = adjustment_params.normalize_homopolymers and pos in seq2_ext_positions
            # Extension and core partition content positions, so is_core = not is_ext
            is_core1 = not is_ext1
            is_core2 = not is_ext2

            if is_ext1 and is_ext2:
                # Both extensions - show extension marker
                score_chars.append(ext_marker)
            elif is_core1 and is_core2:
                # Both core - show match or mismatch based on core comparison
                score_chars.append(scoring_format.match if cores_match else scoring_format.substitution)
            else:
                # One is extension, one is core - show based on whether cores match
                if cores_match:
                    # Cores match - show extension marker for seq1's extension, match for seq1's core
                    score_chars.append(ext_marker if is_ext1 else scoring_format.match)
                else:
                    # Cores differ - mismatch counted
                    score_chars.append(scoring_format.substitution)
            continue

        # Case 3: seq1 has gap, seq2 has content
        if char1 == '-':
            is_ext2 = adjustment_params.normalize_homopolymers and pos in seq2_ext_positions

            if is_ext2:
                # seq2 is extension - show extension marker
                score_chars.append(ext_marker)
            else:
                # seq2 is core - check if cores match (only when HP normalization enabled)
                if adjustment_params.normalize_homopolymers and cores_match:
                    # Cores match - seq1 (gap) absorbs extension → show extension marker
                    score_chars.append(ext_marker)
                else:
                    # Cores differ or HP normalization disabled - show indel markers
                    if adjustment_params.normalize_indels and seen_core_start:
                        score_chars.append(scoring_format.indel_extension)
                    else:
                        score_chars.append(scoring_format.indel_start)
                        seen_core_start = True
            continue

        # Case 4: seq1 has content, seq2 has gap
        if char2 == '-':
            is_ext1 = adjustment_params.normalize_homopolymers and pos in seq1_ext_positions

            if is_ext1:
                # seq1 is extension - show extension marker
                score_chars.append(ext_marker)
            else:
                # seq1 is core - check if cores match (only when HP normalization enabled)
                if adjustment_params.normalize_homopolymers and cores_match:
                    # Cores match - seq1 (core) is the matching core → show match marker
                    score_chars.append(scoring_format.match)
                else:
                    # Cores differ or HP normalization disabled - show indel markers
                    if adjustment_params.normalize_indels and seen_core_start:
                        score_chars.append(scoring_format.indel_extension)
                    else:
                        score_chars.append(scoring_format.indel_start)
                        seen_core_start = True
            continue

    return ''.join(score_chars)


def _unified_analysis(seq1_aligned, seq2_aligned, adjustment_params=None):
    """
    Unified analysis of alignment extracting variant range information and computing metrics.

    This is the single source of truth for alignment metrics. Both adjust_gaps=True
    and adjust_gaps=False code paths use this function, ensuring consistent metrics
    regardless of output mode.

    Returns an AlignmentAnalysis object with all computed metrics and detailed
    variant range information needed to generate either "stuffed" or "adjusted"
    output alignments.

    Args:
        seq1_aligned (str): First sequence with gaps ('-') inserted
        seq2_aligned (str): Second sequence with gaps ('-') inserted
        adjustment_params (AdjustmentParams, optional): Parameters controlling analysis.
                                                       Defaults to DEFAULT_ADJUSTMENT_PARAMS.

    Returns:
        AlignmentAnalysis: Analysis containing metrics and variant range information
    """
    if adjustment_params is None:
        adjustment_params = DEFAULT_ADJUSTMENT_PARAMS

    # Input validation
    if len(seq1_aligned) != len(seq2_aligned):
        raise ValueError(f"Aligned sequences must have same length: seq1={len(seq1_aligned)}, seq2={len(seq2_aligned)}")

    total_alignment_length = len(seq1_aligned)
    edits = 0
    scored_positions = 0

    # Calculate coverage (same as _score_alignment_impl)
    seq1_coverage_positions = 0
    seq2_coverage_positions = 0
    seq1_total_length = 0
    seq2_total_length = 0

    # Find alignment bounds
    alignment_start = 0
    alignment_end = total_alignment_length - 1

    for pos in range(total_alignment_length):
        if seq1_aligned[pos] != '-' and seq2_aligned[pos] != '-':
            alignment_start = pos
            break

    for pos in range(total_alignment_length - 1, -1, -1):
        if seq1_aligned[pos] != '-' and seq2_aligned[pos] != '-':
            alignment_end = pos
            break

    for pos in range(total_alignment_length):
        if seq1_aligned[pos] != '-':
            seq1_total_length += 1
        if seq2_aligned[pos] != '-':
            seq2_total_length += 1
        if alignment_start <= pos <= alignment_end:
            if seq1_aligned[pos] != '-':
                seq1_coverage_positions += 1
            if seq2_aligned[pos] != '-':
                seq2_coverage_positions += 1

    # Find scoring region
    scoring_start, scoring_end = _find_scoring_region(
        seq1_aligned, seq2_aligned, adjustment_params.end_skip_distance
    )

    # Find variant ranges
    raw_variant_ranges = _find_variant_ranges(
        seq1_aligned, seq2_aligned, scoring_start, scoring_end,
        adjustment_params.handle_iupac_overlap
    )

    # Process variant ranges and collect detailed info
    variant_range_infos = []

    for vr_start, vr_end, left_bound, right_bound in raw_variant_ranges:
        # Extract alleles
        allele1, allele1_positions = _extract_allele(seq1_aligned, vr_start, vr_end)
        allele2, allele2_positions = _extract_allele(seq2_aligned, vr_start, vr_end)

        # Get context
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

        # Score the variant range
        vr_score = _score_variant_range(
            allele1, analysis1, allele2, analysis2, adjustment_params
        )

        edits += vr_score['edits']
        scored_positions += vr_score['scored_positions']

        # Create VariantRangeInfo
        vr_info = VariantRangeInfo(
            start=vr_start,
            end=vr_end,
            left_bound_pos=left_bound,
            right_bound_pos=right_bound,
            allele1=allele1,
            allele2=allele2,
            allele1_positions=tuple(allele1_positions),
            allele2_positions=tuple(allele2_positions),
            analysis1=analysis1,
            analysis2=analysis2,
            score_result=vr_score
        )
        variant_range_infos.append(vr_info)

    # Count match positions (same logic as _score_alignment_impl)
    pos = scoring_start
    vr_index = 0
    num_variant_ranges = len(raw_variant_ranges)

    while pos <= scoring_end:
        if vr_index < num_variant_ranges and pos == raw_variant_ranges[vr_index][0]:
            # Skip variant range (already counted above)
            pos = raw_variant_ranges[vr_index][1] + 1
            vr_index += 1
        else:
            char1, char2 = seq1_aligned[pos], seq2_aligned[pos]
            if char1 == '-' and char2 == '-':
                # Dual-gap - not scored
                pass
            else:
                scored_positions += 1
            pos += 1

    # Calculate coverage
    seq1_coverage = seq1_coverage_positions / seq1_total_length if seq1_total_length > 0 else 0.0
    seq2_coverage = seq2_coverage_positions / seq2_total_length if seq2_total_length > 0 else 0.0

    # Calculate identity
    identity = 1.0 - (edits / scored_positions) if scored_positions > 0 else 1.0

    return AlignmentAnalysis(
        identity=identity,
        mismatches=edits,
        scored_positions=scored_positions,
        seq1_coverage=seq1_coverage,
        seq2_coverage=seq2_coverage,
        scoring_start=scoring_start,
        scoring_end=scoring_end,
        variant_ranges=tuple(variant_range_infos)
    )


def _adjust_variant_range(vr_info, adjustment_params):
    """
    Adjust a single variant range so gaps match the scoring analysis.

    Uses middle-padded alignment:
    - Left extensions: pad shorter with gaps at END (closest to core)
    - Core: left-align, pad shorter with gaps at END
    - Right extensions: pad shorter with gaps at START (closest to core)

    Args:
        vr_info: VariantRangeInfo with allele analyses
        adjustment_params: AdjustmentParams for scoring behavior

    Returns:
        tuple: (adjusted_seq1_part, adjusted_seq2_part) strings for this variant range
    """
    analysis1 = vr_info.analysis1
    analysis2 = vr_info.analysis2
    allele1 = vr_info.allele1
    allele2 = vr_info.allele2

    # Extract extension and core portions
    left1 = analysis1.left_extension_count
    right1 = analysis1.right_extension_count
    left2 = analysis2.left_extension_count
    right2 = analysis2.right_extension_count

    ext1_left = allele1[:left1] if left1 > 0 else ""
    ext1_right = allele1[len(allele1) - right1:] if right1 > 0 else ""
    core1 = analysis1.core_content

    ext2_left = allele2[:left2] if left2 > 0 else ""
    ext2_right = allele2[len(allele2) - right2:] if right2 > 0 else ""
    core2 = analysis2.core_content

    result_seq1 = []
    result_seq2 = []

    # Left extensions (middle-padded: gaps at END, closest to core)
    max_left = max(len(ext1_left), len(ext2_left))
    for i in range(max_left):
        c1 = ext1_left[i] if i < len(ext1_left) else '-'
        c2 = ext2_left[i] if i < len(ext2_left) else '-'
        result_seq1.append(c1)
        result_seq2.append(c2)

    # Core content (left-aligned, shorter padded with gaps at end)
    max_core = max(len(core1), len(core2))
    for i in range(max_core):
        c1 = core1[i] if i < len(core1) else '-'
        c2 = core2[i] if i < len(core2) else '-'
        result_seq1.append(c1)
        result_seq2.append(c2)

    # Right extensions (middle-padded: gaps at START, closest to core)
    max_right = max(len(ext1_right), len(ext2_right))
    pad1 = max_right - len(ext1_right)
    pad2 = max_right - len(ext2_right)
    for i in range(max_right):
        c1 = '-' if i < pad1 else ext1_right[i - pad1]
        c2 = '-' if i < pad2 else ext2_right[i - pad2]
        result_seq1.append(c1)
        result_seq2.append(c2)

    return ''.join(result_seq1), ''.join(result_seq2)


def _adjust_alignment_gaps(seq1_aligned, seq2_aligned, analysis, adjustment_params):
    """
    Adjust alignment so gap positions match the variant range scoring analysis.

    This produces a "canonical" alignment where the scoring string is intuitive
    to interpret position-by-position.

    Algorithm:
    1. For positions before scoring region: emit as-is (end-trimmed)
    2. For positions in scoring region:
       - Match positions: emit both characters as-is
       - Dual-gap positions: emit both gaps as-is
       - Variant range: emit adjusted variant range
    3. For positions after scoring region: emit as-is (end-trimmed)

    Args:
        seq1_aligned, seq2_aligned: Original aligned sequences
        analysis: AlignmentAnalysis from _unified_analysis
        adjustment_params: AdjustmentParams for scoring behavior

    Returns:
        tuple: (adjusted_seq1, adjusted_seq2) strings
    """
    result_seq1 = []
    result_seq2 = []

    # Build lookup for variant range start positions
    vr_by_start = {vr.start: vr for vr in analysis.variant_ranges}
    vr_ends = {vr.start: vr.end for vr in analysis.variant_ranges}

    pos = 0
    alignment_length = len(seq1_aligned)

    while pos < alignment_length:
        # Outside scoring region: emit as-is
        if pos < analysis.scoring_start or pos > analysis.scoring_end:
            result_seq1.append(seq1_aligned[pos])
            result_seq2.append(seq2_aligned[pos])
            pos += 1
            continue

        # Check if this position starts a variant range
        if pos in vr_by_start:
            vr_info = vr_by_start[pos]
            vr_end = vr_ends[pos]

            # Emit adjusted variant range
            adj_seq1, adj_seq2 = _adjust_variant_range(vr_info, adjustment_params)
            result_seq1.append(adj_seq1)
            result_seq2.append(adj_seq2)

            # Skip to end of variant range
            pos = vr_end + 1
            continue

        # Match position or dual-gap: emit as-is
        result_seq1.append(seq1_aligned[pos])
        result_seq2.append(seq2_aligned[pos])
        pos += 1

    return ''.join(result_seq1), ''.join(result_seq2)


def _generate_stuffed_score_string(seq1_aligned, seq2_aligned, analysis, adjustment_params, scoring_format):
    """
    Generate "stuffed" score string for original alignment (adjust_gaps=False).

    This generates a scoring visualization where gap positions are "stuffed" with
    markers based on the variant range analysis. The visualization reflects how
    positions were scored using the unified analysis.

    Args:
        seq1_aligned, seq2_aligned: Original aligned sequences (unchanged)
        analysis: AlignmentAnalysis from _unified_analysis
        adjustment_params: AdjustmentParams for scoring behavior
        scoring_format: ScoringFormat for visualization

    Returns:
        str: Scoring visualization string matching original alignment length
    """
    total_length = len(seq1_aligned)
    score_chars = []

    # Add end-trimmed markers for positions before scoring region
    for _ in range(analysis.scoring_start):
        score_chars.append(scoring_format.end_trimmed)

    # Build lookup for variant ranges by start position
    vr_by_start = {vr.start: vr for vr in analysis.variant_ranges}

    # Process positions in scoring region
    pos = analysis.scoring_start
    while pos <= analysis.scoring_end:
        # Check if this position starts a variant range
        if pos in vr_by_start:
            vr_info = vr_by_start[pos]

            # Generate score string for this variant range using existing function
            vr_score_str = _generate_variant_score_string(
                seq1_aligned, seq2_aligned, vr_info.start, vr_info.end,
                vr_info.analysis1, vr_info.analysis2,
                vr_info.allele1_positions, vr_info.allele2_positions,
                scoring_format, adjustment_params
            )
            score_chars.append(vr_score_str)

            # Move past this variant range
            pos = vr_info.end + 1
        else:
            # Match position or dual-gap
            char1, char2 = seq1_aligned[pos], seq2_aligned[pos]

            if char1 == '-' and char2 == '-':
                # Dual-gap (MSA artifact - not scored)
                score_chars.append(scoring_format.dual_gap)
            else:
                # Match position - check for ambiguous match
                is_match, is_ambiguous = _are_nucleotides_equivalent(
                    char1, char2, adjustment_params.handle_iupac_overlap)
                if is_ambiguous:
                    score_chars.append(scoring_format.ambiguous_match)
                else:
                    score_chars.append(scoring_format.match)
            pos += 1

    # Add end-trimmed markers for positions after scoring region
    for _ in range(analysis.scoring_end + 1, total_length):
        score_chars.append(scoring_format.end_trimmed)

    return ''.join(score_chars)


def _generate_stuffed_output(seq1_aligned, seq2_aligned, analysis, adjustment_params, scoring_format):
    """
    Generate AlignmentResult for adjust_gaps=False using original alignment strings.

    The alignment strings are preserved as-is. The score string is generated to
    reflect how positions were scored using the unified analysis.

    Args:
        seq1_aligned, seq2_aligned: Original aligned sequences
        analysis: AlignmentAnalysis from _unified_analysis
        adjustment_params: AdjustmentParams for scoring behavior
        scoring_format: ScoringFormat for visualization

    Returns:
        AlignmentResult: Result with original alignment strings and stuffed score string
    """
    score_aligned = _generate_stuffed_score_string(
        seq1_aligned, seq2_aligned, analysis, adjustment_params, scoring_format
    )

    return AlignmentResult(
        identity=analysis.identity,
        mismatches=analysis.mismatches,
        scored_positions=analysis.scored_positions,
        seq1_coverage=analysis.seq1_coverage,
        seq2_coverage=analysis.seq2_coverage,
        seq1_aligned=seq1_aligned,
        seq2_aligned=seq2_aligned,
        score_aligned=score_aligned
    )


def _generate_adjusted_score_string(adj_seq1, adj_seq2, analysis, adjustment_params, scoring_format):
    """
    Generate intuitive score string for gap-adjusted alignment.

    In the adjusted alignment, gap positions naturally align with the scoring
    analysis, making the visualization more intuitive. This function generates
    a score string where:
    - Extension positions show extension marker (=)
    - Core positions show match (|) or mismatch ( )
    - Match positions show match (|) or ambiguous match (=)

    Args:
        adj_seq1, adj_seq2: Gap-adjusted aligned sequences
        analysis: AlignmentAnalysis from _unified_analysis
        adjustment_params: AdjustmentParams for scoring behavior
        scoring_format: ScoringFormat for visualization

    Returns:
        str: Scoring visualization string for adjusted alignment
    """
    adj_length = len(adj_seq1)
    score_chars = []

    # Adjusted alignment may have different length than original
    # We need to track position mapping through the adjustment

    # For the adjusted alignment, we iterate through positions and generate
    # appropriate markers based on content
    for pos in range(adj_length):
        char1 = adj_seq1[pos]
        char2 = adj_seq2[pos]

        # Dual-gap (should be rare in adjusted alignment within scoring region)
        if char1 == '-' and char2 == '-':
            score_chars.append(scoring_format.dual_gap)
            continue

        # Both have content
        if char1 != '-' and char2 != '-':
            # Check if they match
            is_match, is_ambiguous = _are_nucleotides_equivalent(
                char1, char2, adjustment_params.handle_iupac_overlap)
            if is_match:
                if is_ambiguous:
                    score_chars.append(scoring_format.ambiguous_match)
                else:
                    score_chars.append(scoring_format.match)
            else:
                score_chars.append(scoring_format.substitution)
            continue

        # One has content, one has gap - this is an extension or core indel position
        # In the adjusted alignment, gaps adjacent to matching context are extensions
        if adjustment_params.normalize_homopolymers:
            score_chars.append(scoring_format.homopolymer_extension)
        else:
            score_chars.append(scoring_format.indel_extension)

    return ''.join(score_chars)


def _generate_adjusted_output(seq1_aligned, seq2_aligned, analysis, adjustment_params, scoring_format):
    """
    Generate AlignmentResult for adjust_gaps=True with gap-adjusted alignment.

    The alignment is adjusted so gap positions match the scoring analysis,
    making the visualization more intuitive to interpret position-by-position.

    Args:
        seq1_aligned, seq2_aligned: Original aligned sequences
        analysis: AlignmentAnalysis from _unified_analysis
        adjustment_params: AdjustmentParams for scoring behavior
        scoring_format: ScoringFormat for visualization

    Returns:
        AlignmentResult: Result with adjusted alignment strings and intuitive score string
    """
    # Adjust the alignment gaps using the analysis
    adj_seq1, adj_seq2 = _adjust_alignment_gaps(
        seq1_aligned, seq2_aligned, analysis, adjustment_params
    )

    # Generate score string for the adjusted alignment
    score_aligned = _generate_adjusted_score_string(
        adj_seq1, adj_seq2, analysis, adjustment_params, scoring_format
    )

    return AlignmentResult(
        identity=analysis.identity,
        mismatches=analysis.mismatches,
        scored_positions=analysis.scored_positions,
        seq1_coverage=analysis.seq1_coverage,
        seq2_coverage=analysis.seq2_coverage,
        seq1_aligned=adj_seq1,
        seq2_aligned=adj_seq2,
        score_aligned=score_aligned
    )


def _score_alignment_impl(seq1_aligned, seq2_aligned, adjustment_params=None, scoring_format=None):
    """
    DEPRECATED: Legacy scoring implementation kept for reference.

    This function is no longer used by score_alignment(). The unified architecture
    now uses _unified_analysis() with _generate_stuffed_output() or _generate_adjusted_output()
    to ensure both adjust_gaps=True and adjust_gaps=False produce identical metrics.

    This implementation is preserved for reference and comparison testing.

    ---
    Original docstring:

    Internal implementation: Score alignment and count edits with configurable MycoBLAST-style adjustments.

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

    # Find all variant ranges within the scoring region
    variant_ranges = _find_variant_ranges(
        seq1_aligned, seq2_aligned, scoring_start, scoring_end,
        adjustment_params.handle_iupac_overlap
    )
    num_variant_ranges = len(variant_ranges)  # Cache length

    # Process each position in the scoring region
    pos = scoring_start
    vr_index = 0  # Current variant range index

    while pos <= scoring_end:
        # Check if this position starts a variant range
        if vr_index < num_variant_ranges and pos == variant_ranges[vr_index][0]:
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

            # Generate score string for this variant range
            vr_score_str = _generate_variant_score_string(
                seq1_aligned, seq2_aligned, vr_start, vr_end,
                analysis1, analysis2, allele1_positions, allele2_positions,
                scoring_format, adjustment_params
            )
            score_aligned_seq1.append(vr_score_str)

            # Move past this variant range
            pos = vr_end + 1
            vr_index += 1

        else:
            # This is a match position (both non-gap and equivalent, or dual-gap)
            char1, char2 = seq1_aligned[pos], seq2_aligned[pos]

            # Check for dual-gap (MSA artifact - not scored, just visualized)
            if char1 == '-' and char2 == '-':
                # Dual-gaps are NOT counted in scored_positions
                score_aligned_seq1.append(scoring_format.dual_gap)
                pos += 1
                continue

            is_match, is_ambiguous = _are_nucleotides_equivalent(
                char1, char2, adjustment_params.handle_iupac_overlap)

            scored_positions += 1
            if is_ambiguous:
                score_aligned_seq1.append(scoring_format.ambiguous_match)
            else:
                score_aligned_seq1.append(scoring_format.match)
            pos += 1

    # Add end-trimmed markers for positions after scoring region
    for i in range(scoring_end + 1, total_alignment_length):
        score_aligned_seq1.append(scoring_format.end_trimmed)

    # Calculate coverage as fraction of sequence used in alignment region
    seq1_coverage = seq1_coverage_positions / seq1_total_length if seq1_total_length > 0 else 0.0
    seq2_coverage = seq2_coverage_positions / seq2_total_length if seq2_total_length > 0 else 0.0

    # Create scoring codes strings for visualization
    score_aligned_str = ''.join(score_aligned_seq1)

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
        score_aligned=score_aligned_str
    )


def score_alignment(seq1_aligned, seq2_aligned, adjustment_params=None, scoring_format=None,
                    adjust_gaps=False):
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
        adjust_gaps (bool): If True, rewrite alignment so gap positions match the
                           scoring analysis. Output *_aligned strings may have different
                           length than input. Defaults to False for backward compatibility.

    Returns:
        AlignmentResult: Dataclass containing:
            - identity (float): Identity score based on adjustment parameters
            - mismatches (int): Number of mismatches/edits counted
            - scored_positions (int): Number of positions used for identity calculation
            - seq1_coverage (float): Fraction of seq1 used in scoring region
            - seq2_coverage (float): Fraction of seq2 used in scoring region
            - seq1_aligned (str): Input seq1_aligned (or adjusted if adjust_gaps=True)
            - seq2_aligned (str): Input seq2_aligned (or adjusted if adjust_gaps=True)
            - score_aligned (str): Scoring codes string for visualization

    Note:
        When adjust_gaps=True, the output aligned sequences may have different length
        than the input because variant ranges are rewritten to match the scoring analysis.
        The identity metrics will be identical whether adjust_gaps is True or False,
        but the visualization will be more intuitive with adjust_gaps=True.
    """
    # Use default parameters if none provided
    if adjustment_params is None:
        adjustment_params = DEFAULT_ADJUSTMENT_PARAMS
    if scoring_format is None:
        scoring_format = DEFAULT_SCORING_FORMAT

    # Single analysis pass - computes all metrics using unified logic
    analysis = _unified_analysis(seq1_aligned, seq2_aligned, adjustment_params)

    # Generate output based on adjust_gaps mode
    if adjust_gaps:
        # Adjusted output: gap positions match scoring analysis
        return _generate_adjusted_output(
            seq1_aligned, seq2_aligned, analysis, adjustment_params, scoring_format
        )
    else:
        # Stuffed output: original alignment with stuffed scoring visualization
        return _generate_stuffed_output(
            seq1_aligned, seq2_aligned, analysis, adjustment_params, scoring_format
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


def align_and_score(seq1, seq2, adjustment_params=None, scoring_format=None, adjust_gaps=False):
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
        adjust_gaps (bool): If True, rewrite alignment so gap positions match the
                           scoring analysis. Output *_aligned strings may have different
                           length than input. Defaults to False for backward compatibility.

    Returns:
        AlignmentResult: Dataclass containing:
            - identity (float): Identity score based on adjustment parameters
            - mismatches (int): Number of mismatches/edits counted
            - scored_positions (int): Number of positions used for identity calculation
            - seq1_coverage (float): Fraction of seq1 used in scoring region
            - seq2_coverage (float): Fraction of seq2 used in scoring region
            - seq1_aligned (str): First aligned sequence with gaps (or adjusted if adjust_gaps=True)
            - seq2_aligned (str): Second aligned sequence with gaps (or adjusted if adjust_gaps=True)
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
        seq1_aligned, seq2_aligned, adjustment_params, scoring_format, adjust_gaps
    )

