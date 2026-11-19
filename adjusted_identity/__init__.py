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
from enum import Enum
from functools import lru_cache


class ScoringMode(str, Enum):
    """Controls whether terminal gaps are included in identity scoring.

    LOCAL: Score only the overlap region (default, current behavior).
           Terminal gaps are excluded from identity calculation.
    GLOBAL: Score the full alignment including terminal gaps.
            Uses Needleman-Wunsch (NW) alignment instead of semi-global.
    """
    LOCAL = "local"
    GLOBAL = "global"
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
class AlignmentMetrics:
    """Identity metrics for a pre-aligned pair, without visualization strings.

    Returned by score_alignment_metrics(). Each field equals the same-named
    field of the AlignmentResult that score_alignment() returns for the same
    aligned inputs and AdjustmentParams, whatever its scoring_format or
    adjust_gaps arguments.

    Fields:
        identity: Identity score (0.0-1.0) based on adjustment parameters
        mismatches: Number of mismatches/edits counted
        scored_positions: Number of positions used for identity calculation (denominator)
        seq1_coverage: Fraction of seq1 covered by alignment (0.0-1.0)
        seq2_coverage: Fraction of seq2 covered by alignment (0.0-1.0)
    """
    identity: float
    mismatches: int
    scored_positions: int
    seq1_coverage: float
    seq2_coverage: float


@dataclass(frozen=True)
class ScoringFormat:
    """Format codes for alignment scoring visualization."""
    match: str = '|'                    # Exact match (A=A, C=C, G=G, T=T)
    ambiguous_match: str = '='          # Ambiguous nucleotide match (different IUPAC codes with intersection)
    substitution: str = ' '             # Nucleotide substitution
    indel_start: str = ' '              # First position of indel (scored)
    indel_extension: str = '-'          # Indel positions skipped due to normalization
    homopolymer_extension: str = '='    # Homopolymer length difference
    short_hp_edit: str = ' '            # Short HP length diff counted as edit (below hp_normalize_min_length)
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
        end_skip_distance: End-trimming distance in nucleotides (not alignment positions).
                          Scoring begins/ends at the position where both sequences have
                          reached their end_skip_distance-th nucleotide from that end,
                          so the first/last end_skip_distance - 1 nucleotides of each
                          sequence are excluded from scoring.
                          Only activates when sequences have ≥ 2×end_skip_distance nucleotides.
                          Set to 0 to disable end trimming completely.
        max_repeat_motif_length: Maximum length of repeat motifs to detect (e.g., 2 for dinucleotides).
                                Set to 1 to only detect homopolymers, 2 for dinucleotides, etc.
        hp_normalize_min_length: Minimum homopolymer run length at which HP normalization applies.
                                When min(L1, L2) < this value for a homopolymer (motif length 1) run,
                                the length difference is counted as an edit instead of being normalized.
                                Default 1 normalizes all HP runs (matches pre-existing behavior).
    """
    normalize_homopolymers: bool = True      # Ignore homopolymer length differences
    handle_iupac_overlap: bool = True        # Allow different ambiguity codes to match via intersection
    normalize_indels: bool = True            # Count contiguous indels as single events
    end_skip_distance: int = 0               # Nucleotides to skip from each end (0 = disabled by default)
    max_repeat_motif_length: int = 2         # Maximum repeat motif length to detect (default: dinucleotides)
    scoring_mode: ScoringMode = ScoringMode.LOCAL  # LOCAL scores overlap only; GLOBAL includes terminal gaps
    hp_normalize_min_length: int = 1         # Min HP run length for normalization (1 = normalize all)

    def __post_init__(self):
        """Validate parameter combinations and coerce types."""
        # Coerce string to ScoringMode enum
        if isinstance(self.scoring_mode, str) and not isinstance(self.scoring_mode, ScoringMode):
            object.__setattr__(self, 'scoring_mode', ScoringMode(self.scoring_mode))
        if self.normalize_homopolymers and self.max_repeat_motif_length < 1:
            raise ValueError(
                f"Contradictory configuration: normalize_homopolymers=True requires "
                f"max_repeat_motif_length >= 1 to detect homopolymers, got {self.max_repeat_motif_length}"
            )
        if self.hp_normalize_min_length < 1:
            raise ValueError(
                f"hp_normalize_min_length must be >= 1, got {self.hp_normalize_min_length}"
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
        left_motif_length: Length of motif used for left extension (1 = homopolymer, 0 = no extension)
        right_motif_length: Length of motif used for right extension (1 = homopolymer, 0 = no extension)
    """
    left_extension_count: int
    right_extension_count: int
    core_content: str
    is_pure_extension: bool
    left_motif_length: int = 0
    right_motif_length: int = 0


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
        score_result: Scoring result dict with keys 'edits', 'scored_positions',
            'both_pure_extension' and 'short_hp_edit'
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

    Results are memoized: the character alphabet is tiny, and this function is
    called in scoring hot loops, so lookups hit the cache almost always.

    Args:
        nuc1 (str): First nucleotide (single character)
        nuc2 (str): Second nucleotide (single character)
        enable_iupac_intersection (bool): Allow different ambiguity codes to match via intersection

    Returns:
        tuple: (is_match, is_ambiguous) where:
            - is_match: True if nucleotides are equivalent
            - is_ambiguous: True if match involves ambiguity codes (not exact A=A, C=C, G=G, T=T)
    """
    key = (nuc1, nuc2, enable_iupac_intersection)
    cached = _NUC_EQUIV_CACHE.get(key)
    if cached is None:
        cached = _compute_nucleotides_equivalent(nuc1, nuc2, enable_iupac_intersection)
        _NUC_EQUIV_CACHE[key] = cached
    return cached


_NUC_EQUIV_CACHE = {}


def _compute_nucleotides_equivalent(nuc1, nuc2, enable_iupac_intersection):
    """Uncached implementation backing _are_nucleotides_equivalent."""
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


_CIGAR_OP_REVERSED = re.compile(r'([=XIDM])(\d+)')


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

    # Only the trailing operations matter, so read the CIGAR backwards and stop
    # at the first match/mismatch operation instead of tokenizing all of it.
    # In the reversed string each operation reads as its type followed by its
    # (reversed) length.
    gap_length = 0
    gap_type = None  # 'D' for deletions (gaps in query), 'I' for insertions (gaps in target)

    for cigar_op in _CIGAR_OP_REVERSED.finditer(cigar_string[::-1]):
        op = cigar_op.group(1)
        if op == 'D' or op == 'I':
            # First gap operation sets the type
            if gap_type is None:
                gap_type = op
            # Mixed gap types means we've hit the end of contiguous gaps
            elif gap_type != op:
                # Reset since we found mixed types - not a clean suffix gap
                gap_length = 0
                break
            gap_length += int(cigar_op.group(2)[::-1])
        else:
            # Match/mismatch operations end the gap region
            break

    if gap_length > 0:
        gap_in_query = (gap_type == 'D')  # D = deletion = gap in query
        return (gap_length, gap_in_query)

    return None


# =============================================================================
# Column classification
# =============================================================================
#
# Scoring starts by labelling every alignment column with one class byte. The
# scoring region, the variant ranges, repeat context, match counting and the
# score string are then all read from those labels with C-level bytes methods
# (find, count, translate, slicing) instead of a Python loop over the columns.
#
#   M  both sequences have the same standard nucleotide (exact match)
#   N  both have the same non-standard symbol, e.g. R/R (a match, drawn as ambiguous)
#   Q  different symbols that are equivalent, e.g. R/A (a match, drawn as ambiguous)
#   X  both have content and it is not equivalent
#   I  gap in seq1 only
#   J  gap in seq2 only
#   D  gap in both (dual gap)
#
# M, N and Q are "match" columns; X, I, J and D are "variant" columns. Repeat
# context may only be read from M and N columns (where the two sequences hold
# the same symbol).

_CLS_MATCH, _CLS_SAME_AMBIGUOUS, _CLS_EQUIVALENT, _CLS_MISMATCH, _CLS_GAP1, _CLS_GAP2, _CLS_DUAL_GAP = b'MNQXIJD'
_GAP_CLASSES = b'IJD'

# The symbols the table-driven classifier understands (upper or lower case),
# besides the gap character. Each gets a number 1-15; a gap is 0.
_FAST_SYMBOLS = 'ACGTRYSWKMBDHVN'
_UNENCODABLE = 0xFF

_CLASS_TABLES = None


def _build_class_tables():
    """
    Build the lookup tables behind _classify_columns (once, on first use).

    Returns (and stores in _CLASS_TABLES) a tuple of four 256-byte tables:
        [0] symbol-pair byte -> class byte, handle_iupac_overlap=False
        [1] symbol-pair byte -> class byte, handle_iupac_overlap=True
        [2] seq1 character -> symbol number in the high nibble
        [3] seq2 character -> symbol number in the low nibble
    Characters outside the supported alphabet encode as _UNENCODABLE.

    The class tables are filled in by asking _are_nucleotides_equivalent about
    every symbol pair, so that function stays the single definition of
    equivalence.
    """
    global _CLASS_TABLES

    encode_high = bytearray([_UNENCODABLE]) * 256
    encode_low = bytearray([_UNENCODABLE]) * 256
    encode_high[ord('-')] = encode_low[ord('-')] = 0
    for number, symbol in enumerate(_FAST_SYMBOLS, 1):
        for char in (symbol, symbol.lower()):
            encode_high[ord(char)] = number << 4
            encode_low[ord(char)] = number

    class_tables = []
    for handle_iupac in (False, True):
        table = bytearray(256)
        for pair in range(256):
            number1, number2 = pair >> 4, pair & 15
            if number1 == 0 and number2 == 0:
                column_class = _CLS_DUAL_GAP
            elif number1 == 0:
                column_class = _CLS_GAP1
            elif number2 == 0:
                column_class = _CLS_GAP2
            else:
                is_match, is_ambiguous = _are_nucleotides_equivalent(
                    _FAST_SYMBOLS[number1 - 1], _FAST_SYMBOLS[number2 - 1], handle_iupac)
                if not is_match:
                    column_class = _CLS_MISMATCH
                elif not is_ambiguous:
                    column_class = _CLS_MATCH
                elif number1 == number2:
                    column_class = _CLS_SAME_AMBIGUOUS
                else:
                    column_class = _CLS_EQUIVALENT
            table[pair] = column_class
        class_tables.append(bytes(table))

    _CLASS_TABLES = (class_tables[0], class_tables[1], bytes(encode_high), bytes(encode_low))
    return _CLASS_TABLES


def _classify_columns(seq1_aligned, seq2_aligned, handle_iupac):
    """
    Label every alignment column with a class byte (see the table above).

    Args:
        seq1_aligned, seq2_aligned: Aligned sequences with gaps (same length)
        handle_iupac: Whether different ambiguity codes may match via intersection

    Returns:
        bytes: One class byte per column
    """
    if (type(seq1_aligned) is str and type(seq2_aligned) is str
            and seq1_aligned.isascii() and seq2_aligned.isascii()):
        tables = _CLASS_TABLES or _build_class_tables()
        encoded1 = seq1_aligned.encode('ascii').translate(tables[2])
        encoded2 = seq2_aligned.encode('ascii').translate(tables[3])
        if _UNENCODABLE not in encoded1 and _UNENCODABLE not in encoded2:
            # seq1's symbol numbers sit in the high nibbles and seq2's in the
            # low nibbles, so OR-ing the two as big integers gives one byte per
            # column that identifies the symbol pair. A table lookup turns that
            # into the column class.
            paired = int.from_bytes(encoded1, 'big') | int.from_bytes(encoded2, 'big')
            return paired.to_bytes(len(encoded1), 'big').translate(
                tables[1 if handle_iupac else 0])

    return _classify_columns_generic(seq1_aligned, seq2_aligned, handle_iupac)


def _classify_columns_generic(seq1_aligned, seq2_aligned, handle_iupac):
    """
    Column-by-column equivalent of _classify_columns.

    Handles input the table-driven path cannot encode: characters outside the
    IUPAC alphabet, non-ASCII text and str subclasses.
    """
    classes = bytearray()
    append = classes.append

    for char1, char2 in zip(seq1_aligned, seq2_aligned):
        if char1 == '-':
            append(_CLS_DUAL_GAP if char2 == '-' else _CLS_GAP1)
        elif char2 == '-':
            append(_CLS_GAP2)
        elif char1 == char2 and char1 in _STANDARD_NUCS:
            append(_CLS_MATCH)
        else:
            is_match, is_ambiguous = _are_nucleotides_equivalent(char1, char2, handle_iupac)
            if not is_match:
                append(_CLS_MISMATCH)
            elif not is_ambiguous:
                append(_CLS_MATCH)
            elif char1 == char2 or char1.upper() == char2.upper():
                append(_CLS_SAME_AMBIGUOUS)
            else:
                append(_CLS_EQUIVALENT)

    return bytes(classes)


def _class_translation(mapping, default):
    """Build a 256-byte translation table from {class byte: output byte}."""
    table = bytearray([default]) * 256
    for column_class, output in mapping.items():
        table[column_class] = output
    return bytes(table)


_VARIANT = ord('V')

# Two-letter views of the class string, used to walk variant ranges with
# bytes.find: 'M' at match columns and 'V' at variant columns.
_ANY_VARIANT_TABLE = _class_translation(
    {_CLS_MISMATCH: _VARIANT, _CLS_GAP1: _VARIANT, _CLS_GAP2: _VARIANT, _CLS_DUAL_GAP: _VARIANT},
    _CLS_MATCH)
# Same, but dual-gap columns do not count as variant columns. Used to find
# where a variant range with content starts, skipping dual-gap columns.
_CONTENT_VARIANT_TABLE = _class_translation(
    {_CLS_MISMATCH: _VARIANT, _CLS_GAP1: _VARIANT, _CLS_GAP2: _VARIANT},
    _CLS_MATCH)


def _find_scoring_region(seq1_aligned, seq2_aligned, adjustment_params):
    """
    Find the [start, end] region of the alignment where mismatches should be counted.

    Implements MycoBLAST "digital end trimming" by starting/ending the scoring region
    at the position where both sequences have reached their end_skip_distance-th
    nucleotide (not alignment position) from that end - i.e., the first/last
    end_skip_distance - 1 nucleotides of each sequence are excluded - to avoid
    counting sequencing artifacts near read ends.

    Scoring modes:
    - LOCAL + end_skip_distance=0: Score only overlap region (both sequences have content)
    - LOCAL + end_skip_distance>0: Score overlap minus trimmed ends
    - GLOBAL + end_skip_distance=0: Score full alignment (0, len-1)
    - GLOBAL + end_skip_distance>0: Score full alignment minus trimmed ends

    IMPORTANT: This function counts NUCLEOTIDES (non-gap characters), not alignment positions.
    End trimming only activates when both sequences have >= end_skip_distance nucleotides
    available to skip from each end.

    Args:
        seq1_aligned, seq2_aligned: Aligned sequences with gaps (must be same length)
        adjustment_params: AdjustmentParams controlling scoring_mode and end_skip_distance

    Returns:
        tuple: (scoring_start, scoring_end) - inclusive range of alignment positions to score
    """
    alignment_length = len(seq1_aligned)
    is_global = adjustment_params.scoring_mode == ScoringMode.GLOBAL

    # GLOBAL: start from the full alignment bounds
    scoring_start = 0
    scoring_end = alignment_length - 1

    if not is_global:
        # LOCAL: start from the overlap region (both sequences have content)
        for pos in range(alignment_length):
            if seq1_aligned[pos] != '-' and seq2_aligned[pos] != '-':
                scoring_start = pos
                break

        for pos in range(alignment_length - 1, -1, -1):
            if seq1_aligned[pos] != '-' and seq2_aligned[pos] != '-':
                scoring_end = pos
                break

    return _trim_scoring_region(seq1_aligned, seq2_aligned, scoring_start, scoring_end,
                                adjustment_params.end_skip_distance)


def _trim_scoring_region(seq1_aligned, seq2_aligned, scoring_start, scoring_end,
                         end_skip_distance):
    """
    Apply end trimming to a scoring region.

    From each end of [scoring_start, scoring_end], moves inward to the position
    where both sequences have reached their end_skip_distance-th nucleotide.

    Args:
        seq1_aligned, seq2_aligned: Aligned sequences with gaps
        scoring_start, scoring_end: Region before trimming (inclusive)
        end_skip_distance: Nucleotides to skip from each end (0 = no trimming)

    Returns:
        tuple: (scoring_start, scoring_end) after trimming (inclusive)
    """
    if end_skip_distance == 0:
        return scoring_start, scoring_end

    alignment_length = len(seq1_aligned)

    # From the determined bounds, skip end_skip_distance nucleotides inward
    seq1_count = seq2_count = 0
    trimmed_start = scoring_start
    for pos in range(scoring_start, alignment_length):
        if seq1_aligned[pos] != '-':
            seq1_count += 1
        if seq2_aligned[pos] != '-':
            seq2_count += 1
        if seq1_count >= end_skip_distance and seq2_count >= end_skip_distance:
            trimmed_start = pos
            break

    seq1_count = seq2_count = 0
    trimmed_end = scoring_end
    for pos in range(scoring_end, -1, -1):
        if seq1_aligned[pos] != '-':
            seq1_count += 1
        if seq2_aligned[pos] != '-':
            seq2_count += 1
        if seq1_count >= end_skip_distance and seq2_count >= end_skip_distance:
            trimmed_end = pos
            break

    return trimmed_start, trimmed_end


def _left_context(classes, seq1_aligned, position, length):
    """
    Extract up to 'length' nucleotide characters from match positions before 'position'.

    Context is only valid from positions where both sequences hold the same
    symbol (ignoring case), so it is read from seq1.

    Rules:
    - Same symbol in both (non-gap) → use the character
    - Dual-gap (both '-') → skip (not a match, but not a conflict)
    - Any other case (gap or different symbols) → hit another variant range, stop

    Args:
        classes: Column classes of the alignment (from _classify_columns)
        seq1_aligned: First aligned sequence with gaps
        position: Starting position (exclusive, work backward from here)
        length: Maximum number of nucleotide characters to extract

    Returns:
        str: Context string in left-to-right order (e.g., "AAA") containing up
            to 'length' characters - shorter if the sequence start or another
            variant range (non-match position) is reached first - or None if
            no context characters are available at all.

    Examples (for seq1/seq2 as shown):
        "AGG-AC" / "AG-GAC", position 2, length 1 → 'G'  (position 1 is a match)
        "AGT-AC" / "AX-GAC", position 2, length 1 → None  (position 1 is G vs X)
    """
    context_chars = []
    pos = position - 1
    collected = 0

    while collected < length and pos >= 0:
        column_class = classes[pos]
        if column_class == _CLS_MATCH or column_class == _CLS_SAME_AMBIGUOUS:
            context_chars.append(seq1_aligned[pos])
            collected += 1
        elif column_class != _CLS_DUAL_GAP:
            # We've hit a variant range (gap or mismatch)
            break
        pos -= 1

    if not context_chars:
        return None

    return ''.join(reversed(context_chars))


def _right_context(classes, seq1_aligned, position, length):
    """
    Extract up to 'length' nucleotide characters from match positions after 'position'.

    The mirror image of _left_context: same rules, working forward.

    Args:
        classes: Column classes of the alignment (from _classify_columns)
        seq1_aligned: First aligned sequence with gaps
        position: Starting position (exclusive, work forward from here)
        length: Maximum number of nucleotide characters to extract

    Returns:
        str: Context string in left-to-right order (e.g., "TTT") containing up
            to 'length' characters - shorter if the sequence end or another
            variant range (non-match position) is reached first - or None if
            no context characters are available at all.

    Examples (for seq1/seq2 as shown):
        "AGA--TT" / "AGAT-TT", position 4, length 2 → 'TT'  (positions 5-6 are matches)
        "AGA--TC" / "AGAT-TG", position 4, length 2 → 'T'  (position 6 is C vs G)
    """
    context_chars = []
    pos = position + 1
    max_pos = len(classes)
    collected = 0

    while collected < length and pos < max_pos:
        column_class = classes[pos]
        if column_class == _CLS_MATCH or column_class == _CLS_SAME_AMBIGUOUS:
            context_chars.append(seq1_aligned[pos])
            collected += 1
        elif column_class != _CLS_DUAL_GAP:
            # We've hit a variant range (gap or mismatch)
            break
        pos += 1

    if not context_chars:
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


# Caching
#
# Variant ranges are overwhelmingly short (a substitution, a one-base indel)
# and drawn from a tiny alphabet, so the same work recurs constantly and is
# cached: comparing a chunk of an allele with a motif, scoring a pair of
# alleles in their context, and (further down) generating the output for a
# variant range. Each cache is a bounded LRU of this many entries, keyed on
# plain strings, ints and bools; all of them full hold a few megabytes. Long
# variant ranges rarely repeat, so they skip the caches rather than crowd out
# the short ones.
_CACHE_SIZE = 4096
_VARIANT_CACHE_MAX_ALLELES = 8   # combined length of the two alleles
_OUTPUT_CACHE_MAX_SPAN = 16      # columns in the variant range
_SHARED_ANALYSIS_MAX_CORE = 8    # characters of core content

# AlleleAnalysis is immutable and the same few results come up constantly, so
# equal results share one instance instead of being constructed each time
_allele_analysis = lru_cache(maxsize=_CACHE_SIZE)(AlleleAnalysis)
_EMPTY_ALLELE_ANALYSIS = AlleleAnalysis(0, 0, '', True)
_LEFT_EXTENSION_ANALYSIS = AlleleAnalysis(1, 0, '', True, 1, 0)    # one character extending the left context
_RIGHT_EXTENSION_ANALYSIS = AlleleAnalysis(0, 1, '', True, 0, 1)   # one character extending the right context


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
    # Fast path: identical strings always match (equal characters are
    # equivalent regardless of IUPAC handling)
    if chunk == motif:
        return True
    for c1, c2 in zip(chunk, motif):
        is_match, _ = _are_nucleotides_equivalent(c1, c2, handle_iupac)
        if not is_match:
            return False
    return True


# Motifs are a character or two long, so the same comparisons recur constantly
_motif_matches_cached = lru_cache(maxsize=_CACHE_SIZE)(_motif_matches)


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
        return _EMPTY_ALLELE_ANALYSIS  # Empty allele = pure extension

    n = len(allele)

    # Fast path: single-character alleles (the most common case in divergent
    # alignments - single substitution or single-base indel columns). Motifs
    # longer than one character cannot consume a one-character allele, and the
    # degenerate-motif collapse reduces same-character motifs to length 1, so
    # extension testing reduces to comparing against the adjacent context
    # character on each side.
    if n == 1:
        char = allele[0]
        if left_context and _are_nucleotides_equivalent(
                char, left_context[-1], handle_iupac)[0]:
            return _LEFT_EXTENSION_ANALYSIS
        if right_context and _are_nucleotides_equivalent(
                char, right_context[0], handle_iupac)[0]:
            return _RIGHT_EXTENSION_ANALYSIS
        return _allele_analysis(0, 0, allele, False, 0, 0)

    left_consumed = 0
    right_consumed = 0
    left_motif_len_used = 0
    right_motif_len_used = 0
    motif_matches = _motif_matches_cached if type(handle_iupac) is bool else _motif_matches

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
                if chunk == motif or motif_matches(chunk, motif, handle_iupac):
                    consumed += motif_len
                    pos += motif_len
                else:
                    break

            if consumed > 0:
                left_consumed = consumed
                left_motif_len_used = motif_len
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
                if chunk == motif or motif_matches(chunk, motif, handle_iupac):
                    consumed += motif_len
                    pos -= motif_len
                else:
                    break

            if consumed > 0:
                right_consumed = consumed
                right_motif_len_used = motif_len
                break

    # Core content is what remains (string slicing, no list conversion)
    core_start = left_consumed
    core_end = n - right_consumed
    core_content = allele[core_start:core_end] if core_end > core_start else ''

    if len(core_content) > _SHARED_ANALYSIS_MAX_CORE:
        return AlleleAnalysis(
            left_extension_count=left_consumed,
            right_extension_count=right_consumed,
            core_content=core_content,
            is_pure_extension=False,
            left_motif_length=left_motif_len_used,
            right_motif_length=right_motif_len_used,
        )
    return _allele_analysis(left_consumed, right_consumed, core_content,
                            len(core_content) == 0,
                            left_motif_len_used, right_motif_len_used)


# How far to look on each side of a variant range for the end of a homopolymer
# run before falling back to scanning the rest of the sequence
_HP_RUN_WINDOW = 32


def _hp_run_length(seq_aligned, vr_start, vr_end, base, run_chars=None):
    """
    Count the length of the homopolymer run of `base` in `seq_aligned` that
    spans and touches the variant range [vr_start, vr_end].

    Walks left from vr_start-1 and right from vr_end+1, skipping gap columns
    and counting consecutive non-gap occurrences of `base` (case-insensitive).
    Also counts `base` chars within the variant range (skipping gaps).

    Args:
        seq_aligned: Aligned sequence with '-' for gaps
        vr_start, vr_end: Variant range boundaries (inclusive)
        base: The homopolymer base to measure (single character)
        run_chars: Optional shortcut for a plain ASCII str sequence and a
            non-gap base: the upper- and lower-case forms of `base`, which are
            then the only characters that count as the base. Lets the run be
            measured with str methods instead of a walk over the columns.

    Returns:
        Total number of `base` characters in the run, or 0 if `base` is empty.
    """
    if not base:
        return 0

    if run_chars is not None:
        # Occurrences within the variant range
        inside = seq_aligned[vr_start:vr_end + 1]
        count = inside.count(run_chars[0])
        if run_chars[1] != run_chars[0]:
            count += inside.count(run_chars[1])

        # Run to the left: with gaps removed, the trailing base characters
        before = seq_aligned[max(vr_start - _HP_RUN_WINDOW, 0):vr_start].replace('-', '')
        run = len(before) - len(before.rstrip(run_chars))
        if run == len(before) and vr_start > _HP_RUN_WINDOW:
            # The run reaches the edge of the window: measure it in full
            before = seq_aligned[:vr_start].replace('-', '')
            run = len(before) - len(before.rstrip(run_chars))
        count += run

        # Run to the right: with gaps removed, the leading base characters
        after_start = vr_end + 1
        after = seq_aligned[after_start:after_start + _HP_RUN_WINDOW].replace('-', '')
        run = len(after) - len(after.lstrip(run_chars))
        if run == len(after) and after_start + _HP_RUN_WINDOW < len(seq_aligned):
            after = seq_aligned[after_start:].replace('-', '')
            run = len(after) - len(after.lstrip(run_chars))
        return count + run

    base_upper = base.upper()
    n = len(seq_aligned)
    count = 0

    pos = vr_start - 1
    while pos >= 0:
        c = seq_aligned[pos]
        if c == '-':
            pos -= 1
            continue
        if c.upper() == base_upper:
            count += 1
            pos -= 1
        else:
            break

    for pos in range(vr_start, vr_end + 1):
        c = seq_aligned[pos]
        if c == '-':
            continue
        if c.upper() == base_upper:
            count += 1

    pos = vr_end + 1
    while pos < n:
        c = seq_aligned[pos]
        if c == '-':
            pos += 1
            continue
        if c.upper() == base_upper:
            count += 1
            pos += 1
        else:
            break

    return count


def _short_hp_demote(allele1, analysis1, allele2, analysis2, adjustment_params,
                     seq1_aligned, seq2_aligned, vr_start, vr_end, ascii_sequences=False):
    """
    Decide whether a homopolymer extension should be demoted to a counted edit
    because the shorter side's HP run length is below hp_normalize_min_length.

    Returns True only when motif length 1 (true homopolymer) extensions are
    involved and min(L1, L2) < hp_normalize_min_length for any base used.

    ascii_sequences: pass True when both sequences are plain ASCII str, to
    measure the runs with str methods (see _hp_run_length).
    """
    threshold = adjustment_params.hp_normalize_min_length
    if threshold <= 1:
        return False

    bases_to_check = []

    # Left-side HP base
    left_allele = None
    if analysis1.left_motif_length == 1 and analysis1.left_extension_count > 0:
        left_allele = allele1
    elif analysis2.left_motif_length == 1 and analysis2.left_extension_count > 0:
        left_allele = allele2
    if left_allele:
        bases_to_check.append(left_allele[0])

    # Right-side HP base
    right_allele = None
    if analysis1.right_motif_length == 1 and analysis1.right_extension_count > 0:
        right_allele = allele1
    elif analysis2.right_motif_length == 1 and analysis2.right_extension_count > 0:
        right_allele = allele2
    if right_allele:
        b = right_allele[-1]
        if b not in bases_to_check:
            bases_to_check.append(b)

    if not bases_to_check:
        return False

    for base in bases_to_check:
        # In ASCII, the characters that upper-case to the base's upper-case
        # form are exactly its two case forms
        run_chars = base.upper() + base.lower() if ascii_sequences else None
        # Demote when the shorter of the two runs is below the threshold
        if _hp_run_length(seq1_aligned, vr_start, vr_end, base, run_chars) < threshold:
            return True
        if _hp_run_length(seq2_aligned, vr_start, vr_end, base, run_chars) < threshold:
            return True
    return False


def _score_alleles(allele1, analysis1, allele2, analysis2,
                   normalize_homopolymers, normalize_indels, handle_iupac):
    """
    Score a variant range given two allele analyses using Occam's razor.

    Scoring rules (when normalize_homopolymers=True):
    - Both pure extensions → 0 edits (homopolymer equivalent)
    - One pure extension, other has core → score the core as edits
    - Both have core → compare cores

    When normalize_homopolymers=False:
    - Extensions are treated as regular indels

    The short-HP threshold (hp_normalize_min_length) is not applied here: it
    depends on the sequence around the variant range, so the caller applies it
    on top of this result (see _compute_alignment_metrics).

    Args:
        allele1, allele2: The allele strings
        analysis1, analysis2: AlleleAnalysis for each allele
        normalize_homopolymers: Whether repeat extensions are normalized away
        normalize_indels: Whether contiguous indels count as a single edit
        handle_iupac: Whether to use IUPAC intersection matching

    Returns:
        tuple: (edits, scored_positions, both_pure_extension)
    """
    # When homopolymer normalization is disabled, treat extensions as indels
    if not normalize_homopolymers:
        # Total content from both alleles (treating as a regular indel region)
        total_len = max(len(allele1), len(allele2))
        if total_len == 0:
            return (0, 0, False)
        if normalize_indels:
            return (1, 1, False)
        return (total_len, total_len, False)

    # Both pure extensions -> homopolymer equivalent
    if analysis1.is_pure_extension and analysis2.is_pure_extension:
        return (0, 0, True)

    # One pure extension, other has core -> count core as edits
    if analysis1.is_pure_extension:
        core = analysis2.core_content
        if normalize_indels:
            return (1 if core else 0, 1 if core else 0, False)
        return (len(core), len(core), False)

    if analysis2.is_pure_extension:
        core = analysis1.core_content
        if normalize_indels:
            return (1 if core else 0, 1 if core else 0, False)
        return (len(core), len(core), False)

    # Both have core -> compare cores
    core1, core2 = analysis1.core_content, analysis2.core_content

    if core1 == core2:
        # Cores are identical - just extension differences
        return (0, len(core1), False)

    # Cores differ - compute edit count
    min_len = min(len(core1), len(core2))
    max_len = max(len(core1), len(core2))

    # Count substitutions in overlapping region
    substitutions = sum(1 for i in range(min_len)
                        if not _are_nucleotides_equivalent(core1[i], core2[i], handle_iupac)[0])

    # Handle length difference as indel
    if max_len > min_len:
        if normalize_indels:
            indel_edits = 1
            indel_scored = 1
        else:
            indel_edits = max_len - min_len
            indel_scored = max_len - min_len
    else:
        indel_edits = 0
        indel_scored = 0

    return (substitutions + indel_edits, min_len + indel_scored, False)


class _VariantScore:
    """
    How one variant range scores: the analysis of its two alleles and the
    resulting edit count, before any short-HP demotion.

    This depends only on the two alleles, the repeat context on each side and
    the scoring parameters, and the same few combinations recur constantly, so
    instances are cached (see _score_variant) and shared between alignments.
    They are never modified after construction.

    shape is what output generation needs to know about the two analyses:
    (left1, right1, left2, right2, cores_match) - the extension counts of each
    allele and whether the two cores are identical.
    """
    __slots__ = ('analysis1', 'analysis2', 'edits', 'scored_positions',
                 'both_pure_extension', 'shape')


def _new_variant_score(allele1, allele2, left_context, right_context,
                       max_motif_length, handle_iupac,
                       normalize_homopolymers, normalize_indels):
    """
    Analyze and score one variant range.

    Args:
        allele1, allele2: Non-gap content of each sequence in the variant range
        left_context, right_context: Repeat context on each side (or None)
        max_motif_length: Maximum repeat motif length to try
        handle_iupac: Whether to use IUPAC intersection matching
        normalize_homopolymers: Whether repeat extensions are normalized away
        normalize_indels: Whether contiguous indels count as a single edit

    Returns:
        _VariantScore
    """
    analysis1 = _analyze_allele(allele1, left_context, right_context,
                                max_motif_length, handle_iupac)
    analysis2 = _analyze_allele(allele2, left_context, right_context,
                                max_motif_length, handle_iupac)

    variant = _VariantScore()
    variant.analysis1 = analysis1
    variant.analysis2 = analysis2
    (variant.edits, variant.scored_positions,
     variant.both_pure_extension) = _score_alleles(
        allele1, analysis1, allele2, analysis2,
        normalize_homopolymers, normalize_indels, handle_iupac)
    variant.shape = (analysis1.left_extension_count, analysis1.right_extension_count,
                     analysis2.left_extension_count, analysis2.right_extension_count,
                     analysis1.core_content == analysis2.core_content)
    return variant


# Only call the cached form with plain str/None, int and bool arguments
_score_variant = lru_cache(maxsize=_CACHE_SIZE)(_new_variant_score)


# Score codes: an intermediate, format-independent form of the score string.
# Match columns keep their class byte (M, N, Q) and dual gaps keep D; the other
# codes are below. _codes_to_markers turns codes into ScoringFormat markers.
_CODE_MATCH = _CLS_MATCH
_CODE_AMBIGUOUS = _CLS_SAME_AMBIGUOUS
_CODE_DUAL_GAP = _CLS_DUAL_GAP
(_CODE_END_TRIMMED, _CODE_SUBSTITUTION, _CODE_HP_EXTENSION, _CODE_SHORT_HP_EDIT,
 _CODE_INDEL_START, _CODE_INDEL_EXTENSION) = b'TSHPOE'


def _build_variant_codes(pattern, shape, short_hp_demoted,
                         normalize_homopolymers, normalize_indels):
    """
    Generate score codes for the columns of a variant range.

    The visualization reflects how positions were scored:
    - Extension positions show the extension code
    - Core positions show match if cores match, substitution if they differ
    - Gap positions show the extension code if the other sequence is an extension
    - When the short-HP threshold demoted this variant range to a counted edit,
      extension positions use the short-HP-edit code instead

    Args:
        pattern: Class bytes of the variant range's columns (X, I, J or D each)
        shape: (left1, right1, left2, right2, cores_match) from the allele
            analyses - extension counts for each allele and whether the cores
            are identical
        short_hp_demoted: Whether the short-HP threshold demoted this range
        normalize_homopolymers: Whether repeat extensions are normalized away
        normalize_indels: Whether contiguous indels count as a single edit

    Returns:
        bytes: One score code per column
    """
    left1, right1, left2, right2, cores_match = shape

    # Columns holding each allele's characters, in order
    allele1_positions = [pos for pos, column_class in enumerate(pattern)
                         if column_class == _CLS_MISMATCH or column_class == _CLS_GAP2]
    allele2_positions = [pos for pos, column_class in enumerate(pattern)
                         if column_class == _CLS_MISMATCH or column_class == _CLS_GAP1]

    # Build extension position sets for each sequence
    # Note: Extension and core positions form a partition of all content positions,
    # so is_core = not is_ext (no need to build separate core position sets)
    seq1_ext_positions = (set(allele1_positions[:left1]) |
                          set(allele1_positions[-right1:] if right1 > 0 else []))
    seq2_ext_positions = (set(allele2_positions[:left2]) |
                          set(allele2_positions[-right2:] if right2 > 0 else []))

    # Choose extension code based on normalization settings. When the short-HP
    # threshold demoted this variant range to a counted edit, render extension
    # positions as short-HP edits so the visualization matches the scoring.
    if short_hp_demoted:
        ext_code = _CODE_SHORT_HP_EDIT
    elif normalize_homopolymers:
        ext_code = _CODE_HP_EXTENSION
    else:
        ext_code = _CODE_INDEL_EXTENSION

    codes = bytearray()
    seen_core_start = False  # Track if we've seen the first core position (for indel normalization)

    for pos, column_class in enumerate(pattern):
        # Case 1: Dual-gap - both sequences have gap (MSA artifact, not scored)
        if column_class == _CLS_DUAL_GAP:
            codes.append(_CODE_DUAL_GAP)
            continue

        # Case 2: Both have content
        if column_class == _CLS_MISMATCH:
            # When homopolymer normalization is disabled, is_ext is always False (all content is core)
            is_ext1 = normalize_homopolymers and pos in seq1_ext_positions
            is_ext2 = normalize_homopolymers and pos in seq2_ext_positions
            # Extension and core partition content positions, so is_core = not is_ext
            is_core1 = not is_ext1
            is_core2 = not is_ext2

            if is_ext1 and is_ext2:
                # Both extensions - show extension code
                codes.append(ext_code)
            elif is_core1 and is_core2:
                # Both core - show match or mismatch based on core comparison
                codes.append(_CODE_MATCH if cores_match else _CODE_SUBSTITUTION)
            else:
                # One is extension, one is core - show based on whether cores match
                if cores_match:
                    # Cores match - show extension code for seq1's extension, match for seq1's core
                    codes.append(ext_code if is_ext1 else _CODE_MATCH)
                else:
                    # Cores differ - mismatch counted
                    codes.append(_CODE_SUBSTITUTION)
            continue

        # Case 3: seq1 has gap, seq2 has content
        if column_class == _CLS_GAP1:
            is_ext2 = normalize_homopolymers and pos in seq2_ext_positions

            if is_ext2:
                # seq2 is extension - show extension code
                codes.append(ext_code)
            else:
                # seq2 is core - check if cores match (only when HP normalization enabled)
                if normalize_homopolymers and cores_match:
                    # Cores match - seq1 (gap) absorbs extension → show extension code
                    codes.append(ext_code)
                else:
                    # Cores differ or HP normalization disabled - show indel codes
                    if normalize_indels and seen_core_start:
                        codes.append(_CODE_INDEL_EXTENSION)
                    else:
                        codes.append(_CODE_INDEL_START)
                        seen_core_start = True
            continue

        # Case 4: seq1 has content, seq2 has gap
        is_ext1 = normalize_homopolymers and pos in seq1_ext_positions

        if is_ext1:
            # seq1 is extension - show extension code
            codes.append(ext_code)
        else:
            # seq1 is core - check if cores match (only when HP normalization enabled)
            if normalize_homopolymers and cores_match:
                # Cores match - seq1 (core) is the matching core → show match code
                codes.append(_CODE_MATCH)
            else:
                # Cores differ or HP normalization disabled - show indel codes
                if normalize_indels and seen_core_start:
                    codes.append(_CODE_INDEL_EXTENSION)
                else:
                    codes.append(_CODE_INDEL_START)
                    seen_core_start = True

    return bytes(codes)


_variant_codes = lru_cache(maxsize=_CACHE_SIZE)(_build_variant_codes)


def _analyze_alignment(seq1_aligned, seq2_aligned, adjustment_params=None):
    """
    Unified analysis of alignment extracting variant range information and computing metrics.

    Shares its metric computation (_compute_alignment_metrics) with
    score_alignment() and score_alignment_metrics(), ensuring consistent metrics
    regardless of output mode.

    Returns an AlignmentAnalysis object with all computed metrics and detailed
    information on every variant range. A variant range is a maximal contiguous
    run of columns that are not matches, so runs of dual-gap columns are
    reported too (with empty alleles).

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

    scored_ranges = []
    (identity, edits, scored_positions, seq1_coverage, seq2_coverage,
     scoring_start, scoring_end, _) = _compute_alignment_metrics(
        seq1_aligned, seq2_aligned, adjustment_params, scored_ranges,
        include_dual_gap_ranges=True)

    variant_range_infos = []
    for start, end, allele1, allele2, variant, vr_edits, vr_positions, short_hp in scored_ranges:
        _, allele1_positions = _extract_allele(seq1_aligned, start, end - 1)
        _, allele2_positions = _extract_allele(seq2_aligned, start, end - 1)
        variant_range_infos.append(VariantRangeInfo(
            start=start,
            end=end - 1,
            # Bounding match positions (if within the scoring region)
            left_bound_pos=start - 1 if start > scoring_start else -1,
            right_bound_pos=end if end <= scoring_end else -1,
            allele1=allele1,
            allele2=allele2,
            allele1_positions=tuple(allele1_positions),
            allele2_positions=tuple(allele2_positions),
            analysis1=variant.analysis1,
            analysis2=variant.analysis2,
            score_result={
                'edits': vr_edits,
                'scored_positions': vr_positions,
                'both_pure_extension': variant.both_pure_extension and not short_hp,
                'short_hp_edit': short_hp,
            }
        ))

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


def _compute_alignment_metrics(seq1_aligned, seq2_aligned, adjustment_params,
                               scored_ranges=None, include_dual_gap_ranges=False):
    """
    Single source of truth for alignment metrics.

    Used by score_alignment(), score_alignment_metrics() and _analyze_alignment().

    Args:
        seq1_aligned (str): First sequence with gaps ('-') inserted
        seq2_aligned (str): Second sequence with gaps ('-') inserted
        adjustment_params (AdjustmentParams): Parameters controlling analysis (not None)
        scored_ranges (list, optional): If given, a tuple is appended for each
            variant range, in order (needed for output generation):
            (start, end, allele1, allele2, variant, edits, scored_positions, short_hp)
            where columns [start, end) are the range, variant is its
            _VariantScore, and edits/scored_positions/short_hp are its
            contribution after any short-HP demotion.
        include_dual_gap_ranges (bool): Dual-gap columns never affect the
            metrics. By default they are skipped: runs made only of dual-gap
            columns are not visited, and a variant range starts at its first
            column with content. If True, every maximal run of non-match
            columns is visited as a variant range, dual-gap columns included.

    Returns:
        tuple: (identity, mismatches, scored_positions, seq1_coverage, seq2_coverage,
                scoring_start, scoring_end, classes)
            where classes is the column class string from _classify_columns
    """
    # Input validation
    alignment_length = len(seq1_aligned)
    if alignment_length != len(seq2_aligned):
        raise ValueError(f"Aligned sequences must have same length: seq1={len(seq1_aligned)}, seq2={len(seq2_aligned)}")

    handle_iupac = adjustment_params.handle_iupac_overlap
    classes = _classify_columns(seq1_aligned, seq2_aligned, handle_iupac)

    # Find alignment bounds: first and last column where both sequences have
    # content (the whole alignment if there is no such column)
    alignment_end = len(classes.rstrip(_GAP_CLASSES)) - 1
    if alignment_end < 0:
        alignment_start = 0
        alignment_end = alignment_length - 1
    else:
        alignment_start = alignment_length - len(classes.lstrip(_GAP_CLASSES))

    # Calculate coverage counts
    seq1_total_length = alignment_length - seq1_aligned.count('-')
    seq2_total_length = alignment_length - seq2_aligned.count('-')
    span_length = alignment_end + 1 - alignment_start
    if span_length == alignment_length:
        # No overhang at either end: the whole of each sequence is covered
        seq1_coverage_positions = seq1_total_length
        seq2_coverage_positions = seq2_total_length
    else:
        seq1_coverage_positions = span_length - seq1_aligned.count('-', alignment_start, alignment_end + 1)
        seq2_coverage_positions = span_length - seq2_aligned.count('-', alignment_start, alignment_end + 1)

    # Find scoring region (see _find_scoring_region for the rules)
    if adjustment_params.scoring_mode == ScoringMode.GLOBAL:
        scoring_start = 0
        scoring_end = alignment_length - 1
    else:
        scoring_start = alignment_start
        scoring_end = alignment_end
    end_skip_distance = adjustment_params.end_skip_distance
    if end_skip_distance != 0:
        scoring_start, scoring_end = _trim_scoring_region(
            seq1_aligned, seq2_aligned, scoring_start, scoring_end, end_skip_distance)
    region_end = scoring_end + 1

    # Variant ranges are maximal runs of non-match columns within the scoring
    # region. Walk them with bytes.find on two-letter views of the classes:
    # a range starts at the next 'V' of variant_starts and ends at the next 'M'
    # of any_variant.
    any_variant = classes.translate(_ANY_VARIANT_TABLE)
    if include_dual_gap_ranges or _CLS_DUAL_GAP not in classes:
        variant_starts = any_variant
    else:
        # Start each range at its first column with content. Leading dual-gap
        # columns hold no characters, are skipped by the repeat-context and
        # HP-run walks, and are drawn the same inside or outside a range, so
        # leaving them out changes nothing. Runs of only dual-gap columns have
        # two empty alleles and always score 0 edits over 0 positions.
        variant_starts = classes.translate(_CONTENT_VARIANT_TABLE)

    # Every match column in the scoring region is a scored position
    edits = 0
    scored_positions = any_variant.count(b'M', scoring_start, region_end)

    find_start = variant_starts.find
    start = find_start(b'V', scoring_start, region_end)
    if start >= 0:
        find_end = any_variant.find

        # A variant range has a left bounding match (and so left context) only
        # if a match column precedes it within the scoring region
        first_match = find_end(b'M', scoring_start, region_end)
        if first_match < 0:
            first_match = region_end

        max_motif_length = adjustment_params.max_repeat_motif_length
        normalize_homopolymers = adjustment_params.normalize_homopolymers
        normalize_indels = adjustment_params.normalize_indels
        check_short_hp = normalize_homopolymers and adjustment_params.hp_normalize_min_length > 1
        ascii_sequences = (check_short_hp
                           and type(seq1_aligned) is str and type(seq2_aligned) is str
                           and seq1_aligned.isascii() and seq2_aligned.isascii())

        # Context lookups index and slice the alignment directly when
        # max_motif_length is an ordinary positive count (it always is with
        # homopolymer normalization on); otherwise the column walk decides
        direct_context = type(max_motif_length) is int and max_motif_length >= 1
        if direct_context:
            exact_match_run = b'M' * max_motif_length

        # The variant cache is keyed on plain values; anything unusual is
        # scored directly
        cacheable = (type(seq1_aligned) is str and type(seq2_aligned) is str
                     and type(max_motif_length) is int and type(handle_iupac) is bool
                     and type(normalize_homopolymers) is bool and type(normalize_indels) is bool)

        while start >= 0:
            end = find_end(b'M', start, region_end)
            if end < 0:
                end = region_end

            allele1 = seq1_aligned[start:end].replace('-', '')
            allele2 = seq2_aligned[start:end].replace('-', '')

            # Get context: the longest available run of up to max_motif_length
            # characters from matching columns on each side, or None. Left
            # context needs a match column before the range within the scoring
            # region, right context one after it.
            left_context = None
            right_context = None
            allele_length = len(allele1) + len(allele2)
            if allele_length == 0:
                # Only dual gaps: two empty alleles score the same whatever
                # the context
                pass
            elif not direct_context:
                if first_match < start:
                    left_context = _left_context(classes, seq1_aligned, start, max_motif_length)
                if end < region_end:
                    right_context = _right_context(classes, seq1_aligned, end - 1, max_motif_length)
            elif end - start == 1 or (len(allele1) <= 1 and len(allele2) <= 1):
                # A single-character allele is only ever compared with the
                # nearest context character on each side (see _analyze_allele),
                # so one character of context is enough. Keeping it to one also
                # keeps the number of distinct cases small for the cache.
                if first_match < start:
                    column_class = classes[start - 1]
                    if column_class == _CLS_MATCH or column_class == _CLS_SAME_AMBIGUOUS:
                        left_context = seq1_aligned[start - 1]
                    elif column_class == _CLS_DUAL_GAP:
                        left_context = _left_context(classes, seq1_aligned, start, 1)
                if end < region_end:
                    # The range ends at a match column
                    column_class = classes[end]
                    if column_class == _CLS_MATCH or column_class == _CLS_SAME_AMBIGUOUS:
                        right_context = seq1_aligned[end]
            else:
                # The usual case is a run of exact matches on each side, read
                # with one slice; anything else goes to the column walk
                if first_match < start:
                    if (start >= max_motif_length
                            and classes[start - max_motif_length:start] == exact_match_run):
                        left_context = seq1_aligned[start - max_motif_length:start]
                    else:
                        left_context = _left_context(classes, seq1_aligned, start, max_motif_length)
                if end < region_end:
                    if classes[end:end + max_motif_length] == exact_match_run:
                        right_context = seq1_aligned[end:end + max_motif_length]
                    else:
                        right_context = _right_context(classes, seq1_aligned, end - 1, max_motif_length)

            # Analyze the alleles and score the variant range
            if cacheable and allele_length <= _VARIANT_CACHE_MAX_ALLELES:
                variant = _score_variant(
                    allele1, allele2, left_context, right_context,
                    max_motif_length, handle_iupac, normalize_homopolymers, normalize_indels)
            else:
                variant = _new_variant_score(
                    allele1, allele2, left_context, right_context,
                    max_motif_length, handle_iupac, normalize_homopolymers, normalize_indels)
            vr_edits = variant.edits
            vr_positions = variant.scored_positions

            # Short-HP threshold check: if the shorter HP run is below the
            # configured minimum, demote pure-extension variant ranges to a
            # counted edit. Applies to pure-extension cases only (where the
            # length diff would otherwise normalize to 0 edits).
            short_hp = False
            if (check_short_hp
                    and (variant.analysis1.is_pure_extension or variant.analysis2.is_pure_extension)
                    and _short_hp_demote(allele1, variant.analysis1, allele2, variant.analysis2,
                                         adjustment_params, seq1_aligned, seq2_aligned,
                                         start, end - 1, ascii_sequences)):
                short_hp = True
                total_len = max(len(allele1), len(allele2))
                if total_len == 0:
                    vr_edits = vr_positions = 0
                elif normalize_indels:
                    vr_edits = vr_positions = 1
                else:
                    vr_edits = vr_positions = total_len

            edits += vr_edits
            scored_positions += vr_positions

            if scored_ranges is not None:
                scored_ranges.append((start, end, allele1, allele2, variant,
                                      vr_edits, vr_positions, short_hp))

            start = find_start(b'V', end, region_end)

    # Calculate coverage
    seq1_coverage = seq1_coverage_positions / seq1_total_length if seq1_total_length > 0 else 0.0
    seq2_coverage = seq2_coverage_positions / seq2_total_length if seq2_total_length > 0 else 0.0

    # Calculate identity
    identity = 1.0 - (edits / scored_positions) if scored_positions > 0 else 1.0

    return (identity, edits, scored_positions, seq1_coverage, seq2_coverage,
            scoring_start, scoring_end, classes)


# =============================================================================
# Output generation
# =============================================================================

_marker_table_memo = None


def _marker_table(scoring_format):
    """
    Build the translation from score codes to a ScoringFormat's markers.

    Returns:
        A 256-byte table for bytes.translate when every marker is a single
        Latin-1 character, otherwise a dict for str.translate.
    """
    global _marker_table_memo

    # ScoringFormat is immutable, so the table for the format used last time
    # (almost always the default) can be reused
    memo = _marker_table_memo
    if memo is not None and memo[0] is scoring_format:
        return memo[1]

    markers = {
        _CODE_MATCH: scoring_format.match,
        _CODE_AMBIGUOUS: scoring_format.ambiguous_match,
        _CLS_EQUIVALENT: scoring_format.ambiguous_match,
        _CODE_SUBSTITUTION: scoring_format.substitution,
        _CODE_INDEL_START: scoring_format.indel_start,
        _CODE_INDEL_EXTENSION: scoring_format.indel_extension,
        _CODE_HP_EXTENSION: scoring_format.homopolymer_extension,
        _CODE_SHORT_HP_EDIT: scoring_format.short_hp_edit,
        _CODE_END_TRIMMED: scoring_format.end_trimmed,
        _CODE_DUAL_GAP: scoring_format.dual_gap,
    }
    if all(type(marker) is str and len(marker) == 1 and ord(marker) < 256
           for marker in markers.values()):
        table = bytearray(256)
        for code, marker in markers.items():
            table[code] = ord(marker)
        table = bytes(table)
    else:
        table = markers

    if type(scoring_format) is ScoringFormat:
        _marker_table_memo = (scoring_format, table)
    return table


def _codes_to_markers(codes, scoring_format):
    """Turn a bytes string of score codes into the score string for a ScoringFormat."""
    table = _marker_table(scoring_format)
    if type(table) is bytes:
        return codes.translate(table).decode('latin-1')
    return codes.decode('ascii').translate(table)


def _generate_annotated_score_string(classes, scored_ranges, scoring_start, scoring_end,
                                     adjustment_params, scoring_format):
    """
    Generate "annotated" score string for original alignment (adjust_gaps=False).

    This generates a scoring visualization where gap positions are "annotated" with
    markers based on the variant range analysis. The visualization reflects how
    positions were scored using the unified analysis.

    Args:
        classes: Column classes of the alignment (from _compute_alignment_metrics)
        scored_ranges: Variant ranges from _compute_alignment_metrics
        scoring_start, scoring_end: Scoring region boundaries (inclusive)
        adjustment_params: AdjustmentParams for scoring behavior
        scoring_format: ScoringFormat for visualization

    Returns:
        str: Scoring visualization string matching original alignment length
    """
    normalize_homopolymers = adjustment_params.normalize_homopolymers
    normalize_indels = adjustment_params.normalize_indels
    cacheable = type(normalize_homopolymers) is bool and type(normalize_indels) is bool

    # Add end-trimmed codes for positions before scoring region
    pieces = [bytes([_CODE_END_TRIMMED]) * scoring_start]
    append = pieces.append

    # Process the scoring region: variant ranges interleaved with the columns
    # between them. Those in-between columns are matches or dual gaps, whose
    # class bytes already are their score codes.
    seg_start = scoring_start
    for start, end, _, _, variant, _, _, short_hp in scored_ranges:
        if start > seg_start:
            append(classes[seg_start:start])

        if cacheable and end - start <= _OUTPUT_CACHE_MAX_SPAN:
            append(_variant_codes(classes[start:end], variant.shape, short_hp,
                                  normalize_homopolymers, normalize_indels))
        else:
            append(_build_variant_codes(classes[start:end], variant.shape, short_hp,
                                        normalize_homopolymers, normalize_indels))
        seg_start = end

    if seg_start <= scoring_end:
        append(classes[seg_start:scoring_end + 1])

    # Add end-trimmed codes for positions after scoring region
    append(bytes([_CODE_END_TRIMMED]) * (len(classes) - 1 - scoring_end))

    return _codes_to_markers(b''.join(pieces), scoring_format)


def _build_adjusted_variant(allele1, allele2, core1, core2, shape, short_hp_demoted,
                            normalize_homopolymers, normalize_indels, handle_iupac):
    """
    Adjust a single variant range so gaps match the scoring analysis, and
    generate the score codes for the adjusted columns.

    Uses middle-padded alignment:
    - Left extensions: pad shorter with gaps at END (closest to core)
    - Core: left-align, pad shorter with gaps at END
    - Right extensions: pad shorter with gaps at START (closest to core)

    Gap columns are coded according to how they were scored:
    - short-HP demoted range: short-HP edit for all gap columns
    - extension block (left/right of core): extension (not counted)
    - core block: counted indel -> indel start for the first gap column, indel
      extension for the rest when normalize_indels (matching the annotated
      visualization)

    Args:
        allele1, allele2: Non-gap content of each sequence in the variant range
        core1, core2: Core content of each allele from its AlleleAnalysis
        shape: (left1, right1, left2, right2, cores_match) from the allele
            analyses (see _VariantScore)
        short_hp_demoted: Whether the short-HP threshold demoted this range
        normalize_homopolymers: Whether repeat extensions are normalized away
        normalize_indels: Whether contiguous indels count as a single edit
        handle_iupac: Whether to use IUPAC intersection matching

    Returns:
        tuple: (adjusted_seq1_part, adjusted_seq2_part, codes) for this variant range
    """
    # Extract extension and core portions
    # When homopolymer normalization is disabled, treat all content as core
    # (extensions only affect scoring when normalization is enabled, so gap
    # adjustment should not rewrite positions based on unrecognized extensions)
    if normalize_homopolymers:
        left1, right1, left2, right2, _ = shape
    else:
        left1 = right1 = left2 = right2 = 0
        core1 = allele1
        core2 = allele2

    ext1_left = allele1[:left1] if left1 > 0 else ""
    ext1_right = allele1[len(allele1) - right1:] if right1 > 0 else ""

    ext2_left = allele2[:left2] if left2 > 0 else ""
    ext2_right = allele2[len(allele2) - right2:] if right2 > 0 else ""

    # Left extensions (middle-padded: gaps at END, closest to core), then core
    # content (left-aligned, shorter padded with gaps at end), then right
    # extensions (middle-padded: gaps at START, closest to core)
    max_left = max(len(ext1_left), len(ext2_left))
    max_core = max(len(core1), len(core2))
    max_right = max(len(ext1_right), len(ext2_right))
    adjusted1 = (ext1_left.ljust(max_left, '-') + core1.ljust(max_core, '-')
                 + ext1_right.rjust(max_right, '-'))
    adjusted2 = (ext2_left.ljust(max_left, '-') + core2.ljust(max_core, '-')
                 + ext2_right.rjust(max_right, '-'))

    # Where the core block sits within the adjusted columns
    left_ext_len = max(left1, left2)
    core_len = max_core

    codes = bytearray()
    seen_core_indel = False
    for offset, (char1, char2) in enumerate(zip(adjusted1, adjusted2)):
        gap1 = char1 == '-'
        gap2 = char2 == '-'

        if gap1 and gap2:
            codes.append(_CODE_DUAL_GAP)
        elif not gap1 and not gap2:
            # Both have content - check if they match
            is_match, is_ambiguous = _are_nucleotides_equivalent(char1, char2, handle_iupac)
            if is_match:
                codes.append(_CODE_AMBIGUOUS if is_ambiguous else _CODE_MATCH)
            else:
                codes.append(_CODE_SUBSTITUTION)
        elif short_hp_demoted:
            codes.append(_CODE_SHORT_HP_EDIT)
        elif left_ext_len <= offset < left_ext_len + core_len:
            # Core indel - counted as an edit
            if normalize_indels and seen_core_indel:
                codes.append(_CODE_INDEL_EXTENSION)
            else:
                codes.append(_CODE_INDEL_START)
                seen_core_indel = True
        else:
            # Extension block - normalized away, not counted
            codes.append(_CODE_HP_EXTENSION)

    return adjusted1, adjusted2, bytes(codes)


_adjusted_variant = lru_cache(maxsize=_CACHE_SIZE)(_build_adjusted_variant)


# Score codes for columns left as they are in the gap-adjusted alignment, by
# column class. A gap in one sequence there (unaligned flanks / end-trimmed
# regions) reads as a repeat continuation of adjacent context, so it is drawn
# extension-style. Indexed by whether homopolymers are normalized.
_UNADJUSTED_CODE_TABLES = tuple(
    _class_translation({
        _CLS_MATCH: _CODE_MATCH,
        _CLS_SAME_AMBIGUOUS: _CODE_AMBIGUOUS,
        _CLS_EQUIVALENT: _CODE_AMBIGUOUS,
        _CLS_MISMATCH: _CODE_SUBSTITUTION,
        _CLS_GAP1: gap_code,
        _CLS_GAP2: gap_code,
        _CLS_DUAL_GAP: _CODE_DUAL_GAP,
    }, 0)
    for gap_code in (_CODE_INDEL_EXTENSION, _CODE_HP_EXTENSION)
)
_DUAL_GAP_CODE = bytes([_CODE_DUAL_GAP])


def _generate_adjusted_output(seq1_aligned, seq2_aligned, classes, scored_ranges,
                              scoring_start, scoring_end, adjustment_params, scoring_format):
    """
    Adjust alignment so gap positions match the variant range scoring analysis,
    and generate its score string.

    This produces a "canonical" alignment where the scoring string is intuitive
    to interpret position-by-position.

    Algorithm:
    1. For positions before scoring region: emit as-is (end-trimmed)
    2. For positions in scoring region:
       - Match positions: emit both characters as-is
       - Dual-gap positions: dropped
       - Variant range: emit adjusted variant range
    3. For positions after scoring region: emit as-is (end-trimmed)

    In the score string:
    - Extension positions show extension marker (=)
    - Core positions show match (|) or mismatch ( )
    - Match positions show match (|) or ambiguous match (=)
    - Core indel positions (counted as edits) show indel markers
    - Positions inside a short-HP demoted variant range use short_hp_edit

    Args:
        seq1_aligned, seq2_aligned: Original aligned sequences
        classes: Column classes of the alignment (from _compute_alignment_metrics)
        scored_ranges: Variant ranges from _compute_alignment_metrics
        scoring_start, scoring_end: Scoring region boundaries (inclusive)
        adjustment_params: AdjustmentParams for scoring behavior
        scoring_format: ScoringFormat for visualization

    Returns:
        tuple: (adjusted_seq1, adjusted_seq2, score_aligned)
    """
    normalize_homopolymers = adjustment_params.normalize_homopolymers
    normalize_indels = adjustment_params.normalize_indels
    handle_iupac = adjustment_params.handle_iupac_overlap
    cacheable = (type(seq1_aligned) is str and type(seq2_aligned) is str
                 and type(handle_iupac) is bool
                 and type(normalize_homopolymers) is bool and type(normalize_indels) is bool)

    unadjusted_codes = classes.translate(
        _UNADJUSTED_CODE_TABLES[1 if normalize_homopolymers else 0])

    if scoring_start > scoring_end:
        # Nothing is scored, so nothing is adjusted
        return (seq1_aligned[:], seq2_aligned[:],
                _codes_to_markers(unadjusted_codes, scoring_format))

    # Positions before the scoring region: emit as-is
    pieces1 = [seq1_aligned[:scoring_start]]
    pieces2 = [seq2_aligned[:scoring_start]]
    code_pieces = [unadjusted_codes[:scoring_start]]
    has_dual_gaps = _CLS_DUAL_GAP in classes

    def emit_between_ranges(seg_start, seg_end):
        """Emit columns [seg_start, seg_end) of the scoring region that lie
        between variant ranges: matches as-is, dual gaps dropped."""
        if has_dual_gaps and classes.find(_DUAL_GAP_CODE, seg_start, seg_end) >= 0:
            # These columns are matches or dual gaps, so the only gap
            # characters in either sequence are the dual gaps
            pieces1.append(seq1_aligned[seg_start:seg_end].replace('-', ''))
            pieces2.append(seq2_aligned[seg_start:seg_end].replace('-', ''))
            code_pieces.append(unadjusted_codes[seg_start:seg_end].replace(_DUAL_GAP_CODE, b''))
        else:
            pieces1.append(seq1_aligned[seg_start:seg_end])
            pieces2.append(seq2_aligned[seg_start:seg_end])
            code_pieces.append(unadjusted_codes[seg_start:seg_end])

    seg_start = scoring_start
    for start, end, allele1, allele2, variant, _, _, short_hp in scored_ranges:
        if start > seg_start:
            emit_between_ranges(seg_start, start)

        # Emit adjusted variant range
        build = (_adjusted_variant if cacheable and end - start <= _OUTPUT_CACHE_MAX_SPAN
                 else _build_adjusted_variant)
        adjusted1, adjusted2, codes = build(
            allele1, allele2, variant.analysis1.core_content, variant.analysis2.core_content,
            variant.shape, short_hp, normalize_homopolymers, normalize_indels, handle_iupac)
        pieces1.append(adjusted1)
        pieces2.append(adjusted2)
        code_pieces.append(codes)
        seg_start = end

    if seg_start <= scoring_end:
        emit_between_ranges(seg_start, scoring_end + 1)

    # Positions after the scoring region: emit as-is
    pieces1.append(seq1_aligned[scoring_end + 1:])
    pieces2.append(seq2_aligned[scoring_end + 1:])
    code_pieces.append(unadjusted_codes[scoring_end + 1:])

    return (''.join(pieces1), ''.join(pieces2),
            _codes_to_markers(b''.join(code_pieces), scoring_format))


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
    scored_ranges = []
    (identity, edits, scored_positions, seq1_coverage, seq2_coverage,
     scoring_start, scoring_end, classes) = _compute_alignment_metrics(
        seq1_aligned, seq2_aligned, adjustment_params, scored_ranges)

    # Generate output based on adjust_gaps mode
    if adjust_gaps:
        # Adjusted output: gap positions match scoring analysis
        seq1_aligned, seq2_aligned, score_aligned = _generate_adjusted_output(
            seq1_aligned, seq2_aligned, classes, scored_ranges,
            scoring_start, scoring_end, adjustment_params, scoring_format)
    else:
        # Annotated output: original alignment with annotated scoring visualization
        score_aligned = _generate_annotated_score_string(
            classes, scored_ranges, scoring_start, scoring_end,
            adjustment_params, scoring_format)

    return AlignmentResult(
        identity=identity,
        mismatches=edits,
        scored_positions=scored_positions,
        seq1_coverage=seq1_coverage,
        seq2_coverage=seq2_coverage,
        seq1_aligned=seq1_aligned,
        seq2_aligned=seq2_aligned,
        score_aligned=score_aligned
    )


def score_alignment_metrics(seq1_aligned, seq2_aligned, adjustment_params=None):
    """
    Compute identity metrics for a pre-aligned pair without building visualization strings.

    A faster alternative to score_alignment() for callers that only read the
    metrics, e.g. when scoring every pair of rows in a multiple sequence alignment.

    Guarantees:
    - Every field equals the same-named field of score_alignment(seq1_aligned,
      seq2_aligned, adjustment_params) for every AdjustmentParams combination,
      whatever score_alignment's scoring_format or adjust_gaps arguments (the
      metrics never depend on them).
    - Same input validation: raises ValueError if the aligned sequences differ in length.
    - Pure and deterministic: the result depends only on the two strings and
      adjustment_params (a frozen, hashable dataclass), so it can be cached on them.

    Args:
        seq1_aligned (str): First sequence with gaps ('-') inserted
        seq2_aligned (str): Second sequence with gaps ('-') inserted
        adjustment_params (AdjustmentParams, optional): Parameters controlling which adjustments to apply.
                                                       Defaults to DEFAULT_ADJUSTMENT_PARAMS.

    Returns:
        AlignmentMetrics: Dataclass containing identity, mismatches, scored_positions,
        seq1_coverage and seq2_coverage
    """
    if adjustment_params is None:
        adjustment_params = DEFAULT_ADJUSTMENT_PARAMS

    (identity, edits, scored_positions, seq1_coverage, seq2_coverage,
     _, _, _) = _compute_alignment_metrics(seq1_aligned, seq2_aligned, adjustment_params)

    return AlignmentMetrics(
        identity=identity,
        mismatches=edits,
        scored_positions=scored_positions,
        seq1_coverage=seq1_coverage,
        seq2_coverage=seq2_coverage,
    )


_CIGAR_DIGITS = '0123456789'
_CIGAR_DELETE_DIGITS = str.maketrans('', '', _CIGAR_DIGITS)
_CIGAR_OPS_TO_SPACES = str.maketrans('=XDI', '    ')


def _expand_alignment(align_result, query, target):
    """
    Build the gapped query and target strings for an edlib alignment.

    Gives the same two strings as edlib.getNiceAlignment() (which also builds
    a match line this module never reads, one string concatenation per CIGAR
    operation).

    Args:
        align_result: Result of edlib.align(query, target, ..., task="path")
        query, target: The sequences that were aligned

    Returns:
        tuple: (query_aligned, target_aligned)
    """
    cigar = align_result.get('cigar')
    locations = align_result.get('locations')
    if (type(cigar) is str and cigar and locations
            and type(query) is str and type(target) is str):
        operations = cigar.translate(_CIGAR_DELETE_DIGITS)
        lengths = cigar.translate(_CIGAR_OPS_TO_SPACES).split()
        # Well-formed means lengths and operations strictly alternate, with
        # nothing else in the string
        if (len(lengths) == len(operations) and cigar[0] in _CIGAR_DIGITS
                and cigar[-1] in '=XDI' and cigar.isascii()
                and not operations.strip('=XDI')):
            target_pos = locations[0][0]
            if target_pos is None:
                target_pos = 0
            query_pos = 0

            # A gap in one sequence does not interrupt the other, so each
            # aligned string is contiguous slices of its sequence with gap
            # runs in between: the query is cut at deletions, the target at
            # insertions.
            query_parts = []
            target_parts = []
            query_slice_start = query_pos
            target_slice_start = target_pos

            for length, operation in zip(map(int, lengths), operations):
                if operation == '=' or operation == 'X':
                    query_pos += length
                    target_pos += length
                elif operation == 'D':
                    # Target characters opposite a gap in the query
                    query_parts.append(query[query_slice_start:query_pos])
                    query_parts.append('-' * length)
                    query_slice_start = query_pos
                    target_pos += length
                else:
                    # 'I': query characters opposite a gap in the target
                    target_parts.append(target[target_slice_start:target_pos])
                    target_parts.append('-' * length)
                    target_slice_start = target_pos
                    query_pos += length

            query_parts.append(query[query_slice_start:query_pos])
            target_parts.append(target[target_slice_start:target_pos])
            return ''.join(query_parts), ''.join(target_parts)

    # Anything unexpected: let edlib handle (or reject) it
    alignment = edlib.getNiceAlignment(align_result, query, target)
    return alignment['query_aligned'], alignment['target_aligned']


def _align_nw(seq1, seq2):
    """
    Needleman-Wunsch (global) alignment wrapper using edlib.

    Unlike align_edlib_bidirectional, NW alignment guarantees both sequences
    are fully consumed — no bidirectional trimming trick needed.

    Args:
        seq1, seq2: DNA sequences to align

    Returns:
        dict: {'aligned_seq1': str, 'aligned_seq2': str} or None if alignment fails
    """
    if len(seq1) == 0 or len(seq2) == 0:
        return None

    result = edlib.align(seq1, seq2, mode="NW", task="path")

    if result['editDistance'] == -1:
        return None

    aligned_seq1, aligned_seq2 = _expand_alignment(result, seq1, seq2)
    return {
        'aligned_seq1': aligned_seq1,
        'aligned_seq2': aligned_seq2,
    }


def align_edlib_bidirectional(seq1, seq2):
    """
    Multi-stage alignment optimization using CIGAR-based suffix detection.

    Process:
    1. Reverse complement both sequences
    2. Semi-global (HW) alignment with task=path, parse CIGAR for suffix gaps
    3. Trim sequences based on gap info (RC suffix = forward prefix)
    4. Reverse complement back to forward orientation
    5. Semi-global (HW) alignment with task=path for final result
    6. Re-attach trimmed prefix and any seq2 flanks left uncovered by the
       HW alignment, padding the other sequence with gaps

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
    seq2_prefix_trimmed = 0

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

            # If trimming consumed an entire sequence, the sequences share no
            # alignable region (e.g. completely dissimilar inputs where the
            # whole query is a terminal insertion) - treat as alignment failure
            if not rc_seq1 or not rc_seq2:
                return None

    # Step 5: Back to forward orientation. Trimming the RC suffix is equivalent
    # to trimming the forward prefix, so slice the originals instead of paying
    # for two more reverse-complement passes.
    current_seq1 = seq1[seq1_prefix_trimmed:]
    current_seq2 = seq2[seq2_prefix_trimmed:]

    # Step 6: Final forward alignment with task=path
    result = edlib.align(current_seq1, current_seq2, mode="HW", task="path")

    if result['editDistance'] == -1:
        return None  # Sentinel value for failed alignment

    # Step 8: Get nice alignment
    seq1_aligned, seq2_aligned = _expand_alignment(result, current_seq1, current_seq2)

    # Step 9: Re-attach seq2 (target) flanks left uncovered by the HW alignment.
    # HW mode aligns the query to a substring of the target; target regions
    # outside result['locations'] are omitted from the nice alignment entirely,
    # so without this they would vanish from seq2_aligned and inflate
    # seq2_coverage.
    loc_start, loc_end = result['locations'][0]
    if loc_start > 0:
        seq1_aligned = '-' * loc_start + seq1_aligned
        seq2_aligned = current_seq2[:loc_start] + seq2_aligned
    uncovered_tail = len(current_seq2) - 1 - loc_end
    if uncovered_tail > 0:
        seq1_aligned = seq1_aligned + '-' * uncovered_tail
        seq2_aligned = seq2_aligned + current_seq2[loc_end + 1:]

    # Step 10: Re-attach the prefix removed during RC-stage trimming
    # Only one sequence can have its prefix trimmed

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

    # Choose alignment algorithm based on scoring mode
    if adjustment_params.scoring_mode == ScoringMode.GLOBAL:
        align_result = _align_nw(seq1, seq2)
    else:
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

