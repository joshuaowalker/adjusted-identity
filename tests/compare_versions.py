#!/usr/bin/env python3
"""
Utilities for comparing v0.1.x and v0.2.x scoring behavior.

This module provides functions to:
1. Load the v0.1.x scorer from git origin/main
2. Compare scoring results between versions
3. Format differences for debugging
"""

import subprocess
import sys
import tempfile
import importlib.util
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, Dict, Any, Tuple

# Import current version
from adjusted_identity import (
    score_alignment as score_v02x,
    DEFAULT_ADJUSTMENT_PARAMS,
    AdjustmentParams,
)


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


# Cache for loaded v0.1.x module
_v01x_module = None


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


if __name__ == "__main__":
    # Quick test
    print("Loading v0.1.x scorer...")
    scorer = load_v01x_scorer()
    print(f"Loaded: {scorer}")

    # Test comparison
    seq1 = "AAATTTGGG"
    seq2 = "AAATTTGGG"
    comparison = compare_pair(seq1, seq2)
    print(f"\nTest comparison (identical sequences):")
    print(f"  Match: {comparison.match}")
    print(f"  v0.1.x identity: {comparison.v01x_identity}")
    print(f"  v0.2.x identity: {comparison.v02x_identity}")
