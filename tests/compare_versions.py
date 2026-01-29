#!/usr/bin/env python3
"""
Utilities for comparing scoring behavior between v0.2.4 baseline and HEAD.

This module provides functions to:
1. Load the v0.2.4 baseline scorer (pre-adjust_gaps parameter)
2. Compare scoring results between v0.2.4 and HEAD
3. Format differences for debugging
"""

import subprocess
import tempfile
import importlib.util
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, Tuple, List

# Import current version
from adjusted_identity import (
    score_alignment as score_head,
    DEFAULT_ADJUSTMENT_PARAMS,
    AdjustmentParams,
)


@dataclass
class HeadComparison:
    """Result of comparing v0.2.4 baseline (no adjust_gaps) with HEAD (has adjust_gaps)."""

    seq1_aligned: str
    seq2_aligned: str
    params: AdjustmentParams

    # v0.2.4 baseline results (no adjust_gaps parameter)
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
        """True if HEAD (adjust_gaps=False) is identical to v0.2.4 baseline."""
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
        """True if HEAD (adjust_gaps=True) metrics match v0.2.4 baseline."""
        return (
            self.origin_identity == self.head_true_identity and
            self.origin_mismatches == self.head_true_mismatches and
            self.origin_scored_positions == self.head_true_scored_positions
        )

    @property
    def adjust_true_mismatches_match(self) -> bool:
        """True if HEAD (adjust_gaps=True) mismatch count matches v0.2.4 baseline.

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
            "v0.2.4 Results (no adjust_gaps param):",
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
            "v0.2.4 Results:",
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


# Cache for loaded v0.2.4 baseline module
_origin_module = None


def load_origin_scorer():
    """
    Load the score_alignment function from git tag v0.2.4 (the baseline for HEAD comparison).

    Returns a function with the same signature as score_alignment (minus adjust_gaps).
    The module is cached after first load.
    """
    global _origin_module

    if _origin_module is not None:
        return _origin_module.score_alignment

    # Get v0.2.4 version of the module
    try:
        origin_code = subprocess.check_output(
            ["git", "show", "v0.2.4:adjusted_identity/__init__.py"],
            text=True,
            cwd=Path(__file__).parent.parent,
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"Failed to load v0.2.4 version: {e}\n"
            "Make sure the v0.2.4 tag exists (git tag -l 'v0.2*')"
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
        raise RuntimeError(f"Failed to load v0.2.4 baseline module: {e}")


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


def compare_head_pair(
    seq1_aligned: str,
    seq2_aligned: str,
    params: Optional[AdjustmentParams] = None,
) -> HeadComparison:
    """
    Compare v0.2.4 baseline and HEAD scoring for a sequence pair.

    Args:
        seq1_aligned: First aligned sequence (with gaps)
        seq2_aligned: Second aligned sequence (with gaps)
        params: AdjustmentParams to use (default: DEFAULT_ADJUSTMENT_PARAMS)

    Returns:
        HeadComparison with results from v0.2.4 baseline and HEAD (both adjust_gaps modes)
    """
    if params is None:
        params = DEFAULT_ADJUSTMENT_PARAMS

    # Load v0.2.4 baseline scorer (no adjust_gaps parameter)
    score_origin = load_origin_scorer()

    # Run v0.2.4 baseline (no adjust_gaps param)
    result_origin = score_origin(seq1_aligned, seq2_aligned, params)

    # Run HEAD with adjust_gaps=False (should be identical to v0.2.4 baseline)
    result_head_false = score_head(seq1_aligned, seq2_aligned, params, adjust_gaps=False)

    # Run HEAD with adjust_gaps=True (metrics should match, strings may differ)
    result_head_true = score_head(seq1_aligned, seq2_aligned, params, adjust_gaps=True)

    return HeadComparison(
        seq1_aligned=seq1_aligned,
        seq2_aligned=seq2_aligned,
        params=params,
        # v0.2.4 baseline
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
    Compare all pairwise combinations in a FASTA file for HEAD comparison against v0.2.4.

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
