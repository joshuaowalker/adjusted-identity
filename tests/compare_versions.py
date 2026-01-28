#!/usr/bin/env python3
"""
Utilities for comparing scoring behavior across versions.

This module provides functions to:
1. Load the origin/main scorer (pre-adjust_gaps parameter)
2. Compare scoring results between versions
3. Format differences for debugging

Two main comparison types:
- VersionComparison: Compare v0.1.x and v0.2.x behavior
- HeadComparison: Compare origin/main (no adjust_gaps) with HEAD (has adjust_gaps)
"""

import subprocess
import sys
import tempfile
import importlib.util
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, Dict, Any, Tuple, List

# Import current version
from adjusted_identity import (
    score_alignment as score_head,
    DEFAULT_ADJUSTMENT_PARAMS,
    AdjustmentParams,
)

# Alias for backward compatibility
score_v02x = score_head


@dataclass
class VersionComparison:
    """Result of comparing v0.1.x and v0.2.x scoring."""

    seq1_aligned: str
    seq2_aligned: str
    params: AdjustmentParams

    # v0.1.x results
    v01x_identity: float
    v01x_mismatches: int
    v01x_scored_positions: int
    v01x_score_aligned: str

    # v0.2.x results
    v02x_identity: float
    v02x_mismatches: int
    v02x_scored_positions: int
    v02x_score_aligned: str

    @property
    def dual_gap_count(self) -> int:
        """Count of dual-gap positions (marked as '.' in v0.2.x score_aligned)."""
        return self.v02x_score_aligned.count('.')

    @property
    def v01x_adjusted_scored_positions(self) -> int:
        """v0.1.x scored_positions adjusted by excluding dual-gaps."""
        return self.v01x_scored_positions - self.dual_gap_count

    @property
    def v01x_adjusted_identity(self) -> float:
        """v0.1.x identity recalculated excluding dual-gaps from denominator."""
        if self.v01x_adjusted_scored_positions == 0:
            return 1.0 if self.v01x_mismatches == 0 else 0.0
        return 1.0 - (self.v01x_mismatches / self.v01x_adjusted_scored_positions)

    @property
    def match(self) -> bool:
        """True if both versions produce equivalent results.

        Compares:
        1. Mismatch counts must be equal
        2. Identity values must match when v0.1.x is adjusted for dual-gap exclusion
        """
        if self.v01x_mismatches != self.v02x_mismatches:
            return False
        # Compare identity with v0.1.x adjusted for dual-gap exclusion
        return abs(self.v01x_adjusted_identity - self.v02x_identity) < 1e-9

    @property
    def identity_diff(self) -> float:
        """Difference in identity (v0.2.x - v0.1.x adjusted)."""
        return self.v02x_identity - self.v01x_adjusted_identity

    def format_diff(self) -> str:
        """Format a detailed diff report for debugging."""
        lines = [
            "=" * 70,
            "VERSION COMPARISON DIFFERENCE",
            "=" * 70,
            "",
            "Sequences:",
            f"  seq1: {self.seq1_aligned}",
            f"  seq2: {self.seq2_aligned}",
            "",
            "Parameters:",
            f"  normalize_homopolymers: {self.params.normalize_homopolymers}",
            f"  normalize_indels: {self.params.normalize_indels}",
            f"  handle_iupac_overlap: {self.params.handle_iupac_overlap}",
            f"  end_skip_distance: {self.params.end_skip_distance}",
            "",
            "v0.1.x Results:",
            f"  identity:         {self.v01x_identity:.6f}",
            f"  adjusted_identity:{self.v01x_adjusted_identity:.6f}  (excluding {self.dual_gap_count} dual-gaps)",
            f"  mismatches:       {self.v01x_mismatches}",
            f"  scored_positions: {self.v01x_scored_positions} (adjusted: {self.v01x_adjusted_scored_positions})",
            f"  score_aligned:    {self.v01x_score_aligned}",
            "",
            "v0.2.x Results:",
            f"  identity:         {self.v02x_identity:.6f}",
            f"  mismatches:       {self.v02x_mismatches}",
            f"  scored_positions: {self.v02x_scored_positions}",
            f"  dual_gap_count:   {self.dual_gap_count}",
            f"  score_aligned:    {self.v02x_score_aligned}",
            "",
            "Difference (v0.2.x - v0.1.x adjusted):",
            f"  identity_diff:    {self.identity_diff:+.6f}",
            f"  mismatch_diff:    {self.v02x_mismatches - self.v01x_mismatches:+d}",
            "",
            "Alignment Visualization:",
            f"  seq1:         {self.seq1_aligned}",
            f"  seq2:         {self.seq2_aligned}",
            f"  v01x_score:   {self.v01x_score_aligned}",
            f"  v02x_score:   {self.v02x_score_aligned}",
            "=" * 70,
        ]
        return "\n".join(lines)


@dataclass
class HeadComparison:
    """Result of comparing origin/main (no adjust_gaps) with HEAD (has adjust_gaps)."""

    seq1_aligned: str
    seq2_aligned: str
    params: AdjustmentParams

    # origin/main results (no adjust_gaps parameter)
    origin_identity: float
    origin_mismatches: int
    origin_scored_positions: int
    origin_score_aligned: str
    origin_seq1_aligned: str
    origin_seq2_aligned: str

    # HEAD with adjust_gaps=False
    head_false_identity: float
    head_false_mismatches: int
    head_false_scored_positions: int
    head_false_score_aligned: str
    head_false_seq1_aligned: str
    head_false_seq2_aligned: str

    # HEAD with adjust_gaps=True
    head_true_identity: float
    head_true_mismatches: int
    head_true_scored_positions: int
    head_true_score_aligned: str
    head_true_seq1_aligned: str
    head_true_seq2_aligned: str

    @property
    def adjust_false_identical(self) -> bool:
        """True if HEAD (adjust_gaps=False) is identical to origin/main."""
        return (
            self.origin_identity == self.head_false_identity and
            self.origin_mismatches == self.head_false_mismatches and
            self.origin_scored_positions == self.head_false_scored_positions and
            self.origin_score_aligned == self.head_false_score_aligned and
            self.origin_seq1_aligned == self.head_false_seq1_aligned and
            self.origin_seq2_aligned == self.head_false_seq2_aligned
        )

    @property
    def adjust_true_metrics_match(self) -> bool:
        """True if HEAD (adjust_gaps=True) metrics match origin/main."""
        return (
            self.origin_identity == self.head_true_identity and
            self.origin_mismatches == self.head_true_mismatches and
            self.origin_scored_positions == self.head_true_scored_positions
        )

    @property
    def adjust_true_mismatches_match(self) -> bool:
        """True if HEAD (adjust_gaps=True) mismatch count matches origin/main.

        This is a weaker check than adjust_true_metrics_match - it only verifies
        the mismatch count is the same, allowing for small differences in
        scored_positions (and thus identity) that can occur when gap rewriting
        changes the alignment interpretation slightly.
        """
        return self.origin_mismatches == self.head_true_mismatches

    def format_adjust_false_diff(self) -> str:
        """Format a detailed diff report for adjust_gaps=False differences."""
        lines = [
            "=" * 70,
            "HEAD COMPARISON: adjust_gaps=False DIFFERENCE",
            "=" * 70,
            "",
            "Input Sequences:",
            f"  seq1: {self.seq1_aligned}",
            f"  seq2: {self.seq2_aligned}",
            "",
            "Parameters:",
            f"  normalize_homopolymers: {self.params.normalize_homopolymers}",
            f"  normalize_indels: {self.params.normalize_indels}",
            f"  handle_iupac_overlap: {self.params.handle_iupac_overlap}",
            f"  end_skip_distance: {self.params.end_skip_distance}",
            "",
            "origin/main Results (no adjust_gaps param):",
            f"  identity:         {self.origin_identity:.6f}",
            f"  mismatches:       {self.origin_mismatches}",
            f"  scored_positions: {self.origin_scored_positions}",
            f"  score_aligned:    {self.origin_score_aligned}",
            f"  seq1_aligned:     {self.origin_seq1_aligned}",
            f"  seq2_aligned:     {self.origin_seq2_aligned}",
            "",
            "HEAD Results (adjust_gaps=False):",
            f"  identity:         {self.head_false_identity:.6f}",
            f"  mismatches:       {self.head_false_mismatches}",
            f"  scored_positions: {self.head_false_scored_positions}",
            f"  score_aligned:    {self.head_false_score_aligned}",
            f"  seq1_aligned:     {self.head_false_seq1_aligned}",
            f"  seq2_aligned:     {self.head_false_seq2_aligned}",
            "",
            "Differences:",
        ]

        # Show specific differences
        if self.origin_identity != self.head_false_identity:
            lines.append(f"  identity:         {self.origin_identity:.6f} -> {self.head_false_identity:.6f}")
        if self.origin_mismatches != self.head_false_mismatches:
            lines.append(f"  mismatches:       {self.origin_mismatches} -> {self.head_false_mismatches}")
        if self.origin_scored_positions != self.head_false_scored_positions:
            lines.append(f"  scored_positions: {self.origin_scored_positions} -> {self.head_false_scored_positions}")
        if self.origin_score_aligned != self.head_false_score_aligned:
            lines.append(f"  score_aligned:    '{self.origin_score_aligned}' -> '{self.head_false_score_aligned}'")
        if self.origin_seq1_aligned != self.head_false_seq1_aligned:
            lines.append(f"  seq1_aligned:     '{self.origin_seq1_aligned}' -> '{self.head_false_seq1_aligned}'")
        if self.origin_seq2_aligned != self.head_false_seq2_aligned:
            lines.append(f"  seq2_aligned:     '{self.origin_seq2_aligned}' -> '{self.head_false_seq2_aligned}'")

        lines.append("=" * 70)
        return "\n".join(lines)

    def format_adjust_true_diff(self) -> str:
        """Format a detailed diff report for adjust_gaps=True metric differences."""
        lines = [
            "=" * 70,
            "HEAD COMPARISON: adjust_gaps=True METRICS DIFFERENCE",
            "=" * 70,
            "",
            "Input Sequences:",
            f"  seq1: {self.seq1_aligned}",
            f"  seq2: {self.seq2_aligned}",
            "",
            "Parameters:",
            f"  normalize_homopolymers: {self.params.normalize_homopolymers}",
            f"  normalize_indels: {self.params.normalize_indels}",
            f"  handle_iupac_overlap: {self.params.handle_iupac_overlap}",
            f"  end_skip_distance: {self.params.end_skip_distance}",
            "",
            "origin/main Results:",
            f"  identity:         {self.origin_identity:.6f}",
            f"  mismatches:       {self.origin_mismatches}",
            f"  scored_positions: {self.origin_scored_positions}",
            "",
            "HEAD Results (adjust_gaps=True):",
            f"  identity:         {self.head_true_identity:.6f}",
            f"  mismatches:       {self.head_true_mismatches}",
            f"  scored_positions: {self.head_true_scored_positions}",
            "",
            "Note: aligned strings are expected to differ (that's the point of adjust_gaps).",
            f"  origin score:     {self.origin_score_aligned}",
            f"  HEAD score:       {self.head_true_score_aligned}",
            "",
            "Metric Differences:",
        ]

        # Show specific metric differences
        if self.origin_identity != self.head_true_identity:
            lines.append(f"  identity:         {self.origin_identity:.6f} -> {self.head_true_identity:.6f}")
        if self.origin_mismatches != self.head_true_mismatches:
            lines.append(f"  mismatches:       {self.origin_mismatches} -> {self.head_true_mismatches}")
        if self.origin_scored_positions != self.head_true_scored_positions:
            lines.append(f"  scored_positions: {self.origin_scored_positions} -> {self.head_true_scored_positions}")

        lines.append("=" * 70)
        return "\n".join(lines)


# Cache for loaded origin/main module
_origin_module = None

# Cache for loaded v0.1.x module (alias for backward compatibility)
_v01x_module = None


def load_origin_scorer():
    """
    Load the score_alignment function from git origin/main (without adjust_gaps parameter).

    Returns a function with the same signature as score_alignment (minus adjust_gaps).
    The module is cached after first load.
    """
    global _origin_module

    if _origin_module is not None:
        return _origin_module.score_alignment

    # Get origin/main version of the module
    try:
        origin_code = subprocess.check_output(
            ["git", "show", "origin/main:adjusted_identity/__init__.py"],
            text=True,
            cwd=Path(__file__).parent.parent,
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"Failed to load origin/main version: {e}\n"
            "Make sure you have fetched from origin (git fetch origin)"
        )

    # Write to a temporary file and load as module
    with tempfile.NamedTemporaryFile(
        mode='w',
        suffix='.py',
        delete=False,
        prefix='adjusted_identity_origin_'
    ) as f:
        f.write(origin_code)
        temp_path = f.name

    try:
        spec = importlib.util.spec_from_file_location("adjusted_identity_origin", temp_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        _origin_module = module
        return module.score_alignment
    except Exception as e:
        raise RuntimeError(f"Failed to load origin/main module: {e}")


def load_v01x_scorer():
    """
    Load the v0.1.x score_alignment function from git origin/main.

    Returns a function with the same signature as score_alignment.
    The module is cached after first load.
    """
    global _v01x_module

    if _v01x_module is not None:
        return _v01x_module.score_alignment

    # Get origin/main version of the module
    try:
        origin_code = subprocess.check_output(
            ["git", "show", "origin/main:adjusted_identity/__init__.py"],
            text=True,
            cwd=Path(__file__).parent.parent,
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"Failed to load origin/main version: {e}\n"
            "Make sure you have fetched from origin (git fetch origin)"
        )

    # Write to a temporary file and load as module
    with tempfile.NamedTemporaryFile(
        mode='w',
        suffix='.py',
        delete=False,
        prefix='adjusted_identity_v01x_'
    ) as f:
        f.write(origin_code)
        temp_path = f.name

    try:
        spec = importlib.util.spec_from_file_location("adjusted_identity_v01x", temp_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        _v01x_module = module
        return module.score_alignment
    except Exception as e:
        raise RuntimeError(f"Failed to load v0.1.x module: {e}")


def compare_pair(
    seq1_aligned: str,
    seq2_aligned: str,
    params: Optional[AdjustmentParams] = None,
) -> VersionComparison:
    """
    Compare v0.1.x and v0.2.x scoring for a sequence pair.

    Args:
        seq1_aligned: First aligned sequence (with gaps)
        seq2_aligned: Second aligned sequence (with gaps)
        params: AdjustmentParams to use (default: DEFAULT_ADJUSTMENT_PARAMS)

    Returns:
        VersionComparison with results from both versions
    """
    if params is None:
        params = DEFAULT_ADJUSTMENT_PARAMS

    # Load v0.1.x scorer
    score_v01x = load_v01x_scorer()

    # Run both versions
    result_v01x = score_v01x(seq1_aligned, seq2_aligned, params)
    result_v02x = score_v02x(seq1_aligned, seq2_aligned, params)

    return VersionComparison(
        seq1_aligned=seq1_aligned,
        seq2_aligned=seq2_aligned,
        params=params,
        v01x_identity=result_v01x.identity,
        v01x_mismatches=result_v01x.mismatches,
        v01x_scored_positions=result_v01x.scored_positions,
        v01x_score_aligned=result_v01x.score_aligned,
        v02x_identity=result_v02x.identity,
        v02x_mismatches=result_v02x.mismatches,
        v02x_scored_positions=result_v02x.scored_positions,
        v02x_score_aligned=result_v02x.score_aligned,
    )


def read_fasta(filepath: str) -> list:
    """
    Read FASTA file and return list of (name, sequence) tuples.
    """
    sequences = []
    current_name = None
    current_seq = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name is not None:
                    sequences.append((current_name, ''.join(current_seq)))
                current_name = line[1:]
                current_seq = []
            else:
                current_seq.append(line)

        if current_name is not None:
            sequences.append((current_name, ''.join(current_seq)))

    return sequences


def compare_fasta_pairwise(
    fasta_path: str,
    params: Optional[AdjustmentParams] = None,
    max_pairs: Optional[int] = None,
    stop_on_diff: bool = True,
) -> Tuple[int, int, list]:
    """
    Compare all pairwise combinations in a FASTA file.

    Args:
        fasta_path: Path to aligned FASTA file
        params: AdjustmentParams to use
        max_pairs: Maximum number of pairs to compare (None = all)
        stop_on_diff: If True, stop on first difference and return

    Returns:
        Tuple of (total_pairs, diff_count, list_of_differences)
    """
    from itertools import combinations

    sequences = read_fasta(fasta_path)
    differences = []
    total = 0

    for (name1, seq1), (name2, seq2) in combinations(sequences, 2):
        comparison = compare_pair(seq1, seq2, params)
        total += 1

        if not comparison.match:
            differences.append((name1, name2, comparison))
            if stop_on_diff:
                return total, len(differences), differences

        if max_pairs is not None and total >= max_pairs:
            break

    return total, len(differences), differences


def compare_head_pair(
    seq1_aligned: str,
    seq2_aligned: str,
    params: Optional[AdjustmentParams] = None,
) -> HeadComparison:
    """
    Compare origin/main and HEAD scoring for a sequence pair.

    Args:
        seq1_aligned: First aligned sequence (with gaps)
        seq2_aligned: Second aligned sequence (with gaps)
        params: AdjustmentParams to use (default: DEFAULT_ADJUSTMENT_PARAMS)

    Returns:
        HeadComparison with results from origin/main and HEAD (both adjust_gaps modes)
    """
    if params is None:
        params = DEFAULT_ADJUSTMENT_PARAMS

    # Load origin/main scorer (no adjust_gaps parameter)
    score_origin = load_origin_scorer()

    # Run origin/main (no adjust_gaps param)
    result_origin = score_origin(seq1_aligned, seq2_aligned, params)

    # Run HEAD with adjust_gaps=False (should be identical to origin/main)
    result_head_false = score_head(seq1_aligned, seq2_aligned, params, adjust_gaps=False)

    # Run HEAD with adjust_gaps=True (metrics should match, strings may differ)
    result_head_true = score_head(seq1_aligned, seq2_aligned, params, adjust_gaps=True)

    return HeadComparison(
        seq1_aligned=seq1_aligned,
        seq2_aligned=seq2_aligned,
        params=params,
        # origin/main
        origin_identity=result_origin.identity,
        origin_mismatches=result_origin.mismatches,
        origin_scored_positions=result_origin.scored_positions,
        origin_score_aligned=result_origin.score_aligned,
        origin_seq1_aligned=result_origin.seq1_aligned,
        origin_seq2_aligned=result_origin.seq2_aligned,
        # HEAD adjust_gaps=False
        head_false_identity=result_head_false.identity,
        head_false_mismatches=result_head_false.mismatches,
        head_false_scored_positions=result_head_false.scored_positions,
        head_false_score_aligned=result_head_false.score_aligned,
        head_false_seq1_aligned=result_head_false.seq1_aligned,
        head_false_seq2_aligned=result_head_false.seq2_aligned,
        # HEAD adjust_gaps=True
        head_true_identity=result_head_true.identity,
        head_true_mismatches=result_head_true.mismatches,
        head_true_scored_positions=result_head_true.scored_positions,
        head_true_score_aligned=result_head_true.score_aligned,
        head_true_seq1_aligned=result_head_true.seq1_aligned,
        head_true_seq2_aligned=result_head_true.seq2_aligned,
    )


def compare_head_fasta_pairwise(
    fasta_path: str,
    params: Optional[AdjustmentParams] = None,
    max_pairs: Optional[int] = None,
    stop_on_diff: bool = True,
    check_mode: str = "adjust_false",
) -> Tuple[int, int, List[Tuple[str, str, HeadComparison]]]:
    """
    Compare all pairwise combinations in a FASTA file for HEAD comparison.

    Args:
        fasta_path: Path to aligned FASTA file
        params: AdjustmentParams to use
        max_pairs: Maximum number of pairs to compare (None = all)
        stop_on_diff: If True, stop on first difference and return
        check_mode: "adjust_false" to check adjust_gaps=False is identical,
                   "adjust_true" to check adjust_gaps=True metrics match

    Returns:
        Tuple of (total_pairs, diff_count, list_of_differences)
    """
    from itertools import combinations

    sequences = read_fasta(fasta_path)
    differences = []
    total = 0

    for (name1, seq1), (name2, seq2) in combinations(sequences, 2):
        comparison = compare_head_pair(seq1, seq2, params)
        total += 1

        # Check based on mode
        if check_mode == "adjust_false":
            is_match = comparison.adjust_false_identical
        else:  # adjust_true
            is_match = comparison.adjust_true_metrics_match

        if not is_match:
            differences.append((name1, name2, comparison))
            if stop_on_diff:
                return total, len(differences), differences

        if max_pairs is not None and total >= max_pairs:
            break

    return total, len(differences), differences


if __name__ == "__main__":
    # Quick test
    print("Loading origin/main scorer...")
    scorer = load_origin_scorer()
    print(f"Loaded: {scorer}")

    # Test version comparison
    seq1 = "AAATTTGGG"
    seq2 = "AAATTTGGG"
    comparison = compare_pair(seq1, seq2)
    print(f"\nVersion comparison (identical sequences):")
    print(f"  Match: {comparison.match}")
    print(f"  v0.1.x identity: {comparison.v01x_identity}")
    print(f"  v0.2.x identity: {comparison.v02x_identity}")

    # Test HEAD comparison
    head_comparison = compare_head_pair(seq1, seq2)
    print(f"\nHEAD comparison (identical sequences):")
    print(f"  adjust_gaps=False identical: {head_comparison.adjust_false_identical}")
    print(f"  adjust_gaps=True metrics match: {head_comparison.adjust_true_metrics_match}")
    print(f"  origin identity: {head_comparison.origin_identity}")
    print(f"  HEAD (False) identity: {head_comparison.head_false_identity}")
    print(f"  HEAD (True) identity: {head_comparison.head_true_identity}")
