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
        score_aligned: Scoring visualization string showing match/mismatch patterns
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
    end_skip_distance: int = 20              # Nucleotides to skip from each end (0 to disable)
    max_repeat_motif_length: int = 2         # Maximum repeat motif length to detect (default: dinucleotides)

    def __post_init__(self):
        """Validate parameter combinations."""
        if self.normalize_homopolymers and self.max_repeat_motif_length < 1:
            raise ValueError(
                f"Contradictory configuration: normalize_homopolymers=True requires "
                f"max_repeat_motif_length >= 1 to detect homopolymers, got {self.max_repeat_motif_length}"
            )


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
        # Exact match - check if it's one of the standard nucleotides
        is_standard = nuc1 in {'A', 'T', 'C', 'G'}
        return (True, not is_standard)
    
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



def _process_indel_left_right(seq1_aligned, seq2_aligned, start, end, adjustment_params, scoring_format):
    """
    Process an indel region using left-right repeat motif/homopolymer extension algorithm.
    
    This implements the general case algorithm that can handle mixed indels
    by processing repeat motif or homopolymer extensions from both ends, leaving only the
    non-repeat content to be treated as regular indels.
    
    Algorithm:
    1. Extract the full indel content
    2. For left side: try different context lengths to find repeat motifs
    3. For right side: independently try different context lengths  
    4. Process extensions based on detected motif lengths
    5. Generate scoring string and counts based on adjustment parameters
    
    Args:
        seq1_aligned, seq2_aligned: Aligned sequences with gaps
        start, end: Start/end positions of indel region (inclusive)
        adjustment_params: AdjustmentParams controlling how to handle different components
        scoring_format: ScoringFormat for generating visualization symbols
        
    Returns:
        dict: {
            'score_string': str - scoring visualization for this indel region,
            'edits': int - number of edits to count,
            'scored_positions': int - number of positions to add to denominator,
            'all_homopolymer': bool - True if entire indel is homopolymer/repeat extension
        }
    """
    # Extract indel characters in order
    indel_chars = []
    indel_positions = []
    
    for i in range(start, end + 1):
        if seq1_aligned[i] != '-':
            indel_chars.append(seq1_aligned[i])
            indel_positions.append(i)
        elif seq2_aligned[i] != '-':
            indel_chars.append(seq2_aligned[i])
            indel_positions.append(i)
    
    if not indel_chars:
        return {
            'score_string': '',
            'edits': 0,
            'scored_positions': 0,
            'all_homopolymer': False
        }
    
    # Determine which sequence has gaps (the context sequence)
    gap_in_seq1 = (seq1_aligned[start] == '-')
    context_seq = seq2_aligned if gap_in_seq1 else seq1_aligned
    
    # LEFT SIDE PROCESSING - Try different context lengths
    left_repeat_count = 0
    left_context_length = 0
    max_motif_length = adjustment_params.max_repeat_motif_length
    
    for try_length in range(min(max_motif_length, start), 0, -1):
        if start >= try_length:
            left_context = context_seq[start - try_length:start]
            
            # Check for degenerate case (homopolymer disguised as longer motif)
            if len(set(left_context)) == 1:
                # It's actually a homopolymer, use length 1
                left_context = left_context[0]
                actual_length = 1
            else:
                actual_length = try_length
            
            # Count how many complete motifs we can consume from left
            consumed = 0
            while consumed + actual_length <= len(indel_chars):
                if actual_length == 1:
                    # Single character comparison
                    if indel_chars[consumed] == left_context:
                        consumed += 1
                    else:
                        break
                else:
                    # Multi-character motif comparison
                    chunk = ''.join(indel_chars[consumed:consumed + actual_length])
                    if chunk == left_context:
                        consumed += actual_length
                    else:
                        break
            
            if consumed > 0:
                left_repeat_count = consumed
                left_context_length = actual_length
                break  # Found a working context length
    
    # RIGHT SIDE PROCESSING - Independent of left side
    right_repeat_count = 0
    right_context_length = 0
    
    # Calculate remaining indel after left consumption
    remaining_start = left_repeat_count
    remaining_chars = indel_chars[remaining_start:]
    
    for try_length in range(min(max_motif_length, len(context_seq) - end - 1), 0, -1):
        if end + try_length < len(context_seq):
            right_context = context_seq[end + 1:end + 1 + try_length]
            
            # Check for degenerate case
            if len(set(right_context)) == 1:
                right_context = right_context[0]
                actual_length = 1
            else:
                actual_length = try_length
            
            # Count from right side of remaining indel
            consumed = 0
            pos = len(remaining_chars)
            
            while pos >= actual_length:
                if actual_length == 1:
                    # Single character comparison
                    if remaining_chars[pos - 1] == right_context:
                        consumed += 1
                        pos -= 1
                    else:
                        break
                else:
                    # Multi-character motif comparison
                    chunk = ''.join(remaining_chars[pos - actual_length:pos])
                    if chunk == right_context:
                        consumed += actual_length
                        pos -= actual_length
                    else:
                        break
            
            if consumed > 0:
                right_repeat_count = consumed
                right_context_length = actual_length
                break
    
    # Calculate middle (non-repeat) portion
    middle_indel_count = len(indel_chars) - left_repeat_count - right_repeat_count
    
    all_homopolymer = (middle_indel_count == 0)
    
    # Score each component independently and sum results
    total_edits = 0
    total_scored_positions = 0
    score_parts = []
    
    # 1. Score left repeat region (homopolymer extensions)
    if left_repeat_count > 0:
        if adjustment_params.normalize_homopolymers:
            # Homopolymer extensions: no edits, no scored positions
            left_edits = 0
            left_scored_positions = 0
            left_score = scoring_format.homopolymer_extension * left_repeat_count
        else:
            # Treat as regular indel
            if adjustment_params.normalize_indels:
                left_edits = 1
                left_scored_positions = 1
                left_score = scoring_format.indel_start + (scoring_format.indel_extension * (left_repeat_count - 1))
            else:
                left_edits = left_repeat_count
                left_scored_positions = left_repeat_count
                left_score = scoring_format.indel_start * left_repeat_count
        
        total_edits += left_edits
        total_scored_positions += left_scored_positions
        score_parts.append(left_score)
    
    # 2. Score middle indel region
    if middle_indel_count > 0:
        if adjustment_params.normalize_indels:
            middle_edits = 1
            middle_scored_positions = 1
            middle_score = scoring_format.indel_start + (scoring_format.indel_extension * (middle_indel_count - 1))
        else:
            middle_edits = middle_indel_count
            middle_scored_positions = middle_indel_count
            middle_score = scoring_format.indel_start * middle_indel_count
        
        total_edits += middle_edits
        total_scored_positions += middle_scored_positions
        score_parts.append(middle_score)
    
    # 3. Score right repeat region (homopolymer extensions)  
    if right_repeat_count > 0:
        if adjustment_params.normalize_homopolymers:
            # Homopolymer extensions: no edits, no scored positions
            right_edits = 0
            right_scored_positions = 0
            right_score = scoring_format.homopolymer_extension * right_repeat_count
        else:
            # Treat as regular indel
            if adjustment_params.normalize_indels:
                right_edits = 1
                right_scored_positions = 1
                right_score = scoring_format.indel_start + (scoring_format.indel_extension * (right_repeat_count - 1))
            else:
                right_edits = right_repeat_count
                right_scored_positions = right_repeat_count
                right_score = scoring_format.indel_start * right_repeat_count
        
        total_edits += right_edits
        total_scored_positions += right_scored_positions
        score_parts.append(right_score)
    
    return {
        'score_string': ''.join(score_parts),
        'edits': total_edits,
        'scored_positions': total_scored_positions,
        'all_homopolymer': all_homopolymer
    }


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
    score_aligned = []

    # Find the scoring region boundaries (end trimming controlled by end_skip_distance)
    scoring_start, scoring_end = _find_scoring_region(
        seq1_aligned, seq2_aligned, adjustment_params.end_skip_distance
    )

    i = 0

    while i < total_alignment_length:
        seq1_char = seq1_aligned[i]
        seq2_char = seq2_aligned[i]

        # Check if outside scoring region (end mismatch skipping)
        if i < scoring_start or i > scoring_end:
            # Add scoring code for skipped region
            score_aligned.append(scoring_format.end_trimmed)
            i += 1
            continue

        # We're in the scoring region

        # Check for exact match or IUPAC equivalence
        is_equivalent, is_ambiguous = _are_nucleotides_equivalent(seq1_char, seq2_char, adjustment_params.handle_iupac_overlap)

        if is_equivalent:
            # Match - no edit
            scored_positions += 1
            # Use appropriate scoring code based on match type
            if is_ambiguous:
                score_aligned.append(scoring_format.ambiguous_match)
            else:
                score_aligned.append(scoring_format.match)
            i += 1
        elif seq1_char == '-' or seq2_char == '-':
            # Start of an indel region - find end of contiguous gaps in the same sequence
            indel_start = i
            gap_in_seq1 = (seq1_char == '-')  # Track which sequence has the gap

            # Find the end of this indel region (gaps must stay in the same sequence)
            while i < total_alignment_length:
                seq1_has_gap = (seq1_aligned[i] == '-')
                seq2_has_gap = (seq2_aligned[i] == '-')
                
                # Must be an indel position and gap must be in the same sequence as when we started
                if (seq1_has_gap or seq2_has_gap) and (seq1_has_gap == gap_in_seq1):
                    i += 1
                else:
                    break  # Gap switched sequences or we hit a match

            indel_end = i - 1
            indel_length = indel_end - indel_start + 1

            # We'll count scored positions after determining if it's a homopolymer

            # Process indel using left-right homopolymer extension algorithm
            indel_result = _process_indel_left_right(
                seq1_aligned, seq2_aligned, indel_start, indel_end, 
                adjustment_params, scoring_format
            )
            
            # Add results from indel processing
            edits += indel_result['edits']
            scored_positions += indel_result['scored_positions']
            score_aligned.append(indel_result['score_string'])
            
        else:
            # Substitution - count normally (scoring region already checked above)
            scored_positions += 1
            edits += 1
            score_aligned.append(scoring_format.substitution)
            i += 1

    # Calculate coverage as fraction of sequence used in alignment region
    seq1_coverage = seq1_coverage_positions / seq1_total_length if seq1_total_length > 0 else 0.0
    seq2_coverage = seq2_coverage_positions / seq2_total_length if seq2_total_length > 0 else 0.0

    # Create scoring codes string for visualization
    score_aligned_str = ''.join(score_aligned)

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

