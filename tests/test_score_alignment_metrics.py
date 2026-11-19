"""
Tests for score_alignment_metrics(), the supported metrics-only scoring API.

The contract: for the same aligned inputs and AdjustmentParams, every
AlignmentMetrics field equals the same-named field of score_alignment()
(with either adjust_gaps setting), with the same input validation, and without
building any visualization strings.
"""

import itertools
import random
from dataclasses import FrozenInstanceError, fields
from pathlib import Path

import pytest

import adjusted_identity
from adjusted_identity import (
    AdjustmentParams,
    AlignmentMetrics,
    DEFAULT_ADJUSTMENT_PARAMS,
    RAW_ADJUSTMENT_PARAMS,
    ScoringMode,
    align_and_score,
    score_alignment,
    score_alignment_metrics,
)

METRIC_FIELDS = [f.name for f in fields(AlignmentMetrics)]

MSA_FASTA = (Path(__file__).parent / "data" /
             "ONT10.82-B11-IN25-00187-iNat272809801-c1-RiC500-msa.fasta")


def _param_grid():
    """Every combination of the AdjustmentParams knobs, at a few values each."""
    grid = []
    for hp, iupac, indels, end_skip, motif, mode, hp_min in itertools.product(
            [True, False], [True, False], [True, False], [0, 5], [1, 2, 3],
            [ScoringMode.LOCAL, ScoringMode.GLOBAL], [1, 4]):
        grid.append(AdjustmentParams(
            normalize_homopolymers=hp, handle_iupac_overlap=iupac,
            normalize_indels=indels, end_skip_distance=end_skip,
            max_repeat_motif_length=motif, scoring_mode=mode,
            hp_normalize_min_length=hp_min))
    return grid


PARAM_GRID = _param_grid()


def _param_id(p):
    return (f"hp{int(p.normalize_homopolymers)}-iupac{int(p.handle_iupac_overlap)}"
            f"-indel{int(p.normalize_indels)}-skip{p.end_skip_distance}"
            f"-motif{p.max_repeat_motif_length}-{p.scoring_mode.value}"
            f"-hpmin{p.hp_normalize_min_length}")

# The configuration gaphack 0.7.x scores MSA rows with
GAPHACK_PARAMS = AdjustmentParams(
    end_skip_distance=0, normalize_homopolymers=True, handle_iupac_overlap=True,
    normalize_indels=True, max_repeat_motif_length=1)


def _random_base_sequence(rng, length):
    """Random sequence with homopolymer runs and dinucleotide repeats mixed in."""
    parts = []
    total = 0
    while total < length:
        r = rng.random()
        if r < 0.2:
            chunk = rng.choice("ACGT") * rng.randint(2, 9)
        elif r < 0.3:
            chunk = "".join(rng.choice("ACGT") for _ in range(2)) * rng.randint(2, 5)
        else:
            chunk = rng.choice("ACGT")
        parts.append(chunk)
        total += len(chunk)
    return "".join(parts)[:length]


def _mutate(rng, seq):
    """Apply substitutions, HP length changes, indels, IUPAC codes and truncation."""
    out = list(seq)
    for _ in range(rng.randint(0, 12)):
        if not out:
            break
        i = rng.randrange(len(out))
        op = rng.random()
        if op < 0.3:
            out[i] = rng.choice("ACGT")
        elif op < 0.55:
            out.insert(i, out[i])                      # HP extension
        elif op < 0.7:
            del out[i]
        elif op < 0.8:
            out[i:i] = rng.choice("ACGT") * rng.randint(1, 4)
        else:
            out[i] = rng.choice("RYSWKMBDHVN")
    s = "".join(out)
    if rng.random() < 0.25 and len(s) > 40:
        s = s[rng.randint(0, 15):len(s) - rng.randint(0, 15)]
    return s


def _add_dual_gap_columns(rng, a, b):
    """Insert columns that are gaps in both rows, as an MSA would."""
    a, b = list(a), list(b)
    for _ in range(rng.randint(1, 6)):
        i = rng.randint(0, len(a))
        a.insert(i, "-")
        b.insert(i, "-")
    return "".join(a), "".join(b)


def _synthetic_aligned_pairs(n=40, seed=20261119):
    """Deterministic corpus of aligned pairs exercising every scoring feature."""
    rng = random.Random(seed)
    pairs = []
    for k in range(n):
        base = _random_base_sequence(rng, rng.randint(60, 180))
        s1, s2 = _mutate(rng, base), _mutate(rng, base)
        mode = ScoringMode.LOCAL if k % 2 == 0 else ScoringMode.GLOBAL
        r = align_and_score(s1, s2, AdjustmentParams(scoring_mode=mode))
        if r.mismatches < 0:
            continue
        a, b = r.seq1_aligned, r.seq2_aligned
        if k % 3 == 0:
            a, b = _add_dual_gap_columns(rng, a, b)
        pairs.append((a, b))
    return pairs


SYNTHETIC_PAIRS = _synthetic_aligned_pairs()

EDGE_CASE_PAIRS = [
    ("", ""),
    ("-", "-"),
    ("----", "----"),
    ("A", "A"),
    ("A", "T"),
    ("A", "-"),
    ("-A", "A-"),
    ("AAAA-TT", "AAA--TT"),
    ("TGC-C-TC", "TGCT--TC"),
    ("GATAT--G", "GATATATG"),
    ("ACG--T", "ACG--T"),
    ("---ACGT---", "TTTACGTAAA"),
    ("ACRT", "ACGT"),
    ("ACRT", "ACYT"),
    ("NNNN", "ACGT"),
]


def _metrics_of(result):
    return {name: getattr(result, name) for name in METRIC_FIELDS}


def _assert_matches_score_alignment(pairs, params):
    for a, b in pairs:
        m = score_alignment_metrics(a, b, params)
        assert isinstance(m, AlignmentMetrics)
        got = _metrics_of(m)
        for adjust_gaps in (False, True):
            expected = _metrics_of(score_alignment(a, b, params, adjust_gaps=adjust_gaps))
            assert got == expected, (
                f"\nparams={params}\nadjust_gaps={adjust_gaps}\n"
                f"seq1={a}\nseq2={b}\nmetrics={got}\nscore_alignment={expected}")


class TestScoreAlignmentMetricsBasics:
    """Return type, defaults and validation."""

    def test_returns_alignment_metrics_with_expected_fields(self):
        m = score_alignment_metrics("AAATTTGGG", "AAAATTTGG")
        assert isinstance(m, AlignmentMetrics)
        assert METRIC_FIELDS == ["identity", "mismatches", "scored_positions",
                                 "seq1_coverage", "seq2_coverage"]

    def test_no_visualization_fields(self):
        m = score_alignment_metrics("ACGT", "ACGT")
        for name in ("seq1_aligned", "seq2_aligned", "score_aligned",
                     "variant_ranges", "scoring_start", "scoring_end"):
            assert not hasattr(m, name)

    def test_perfect_match(self):
        m = score_alignment_metrics("ACGTACGT", "ACGTACGT")
        assert m == AlignmentMetrics(identity=1.0, mismatches=0, scored_positions=8,
                                     seq1_coverage=1.0, seq2_coverage=1.0)

    def test_substitution(self):
        m = score_alignment_metrics("ACGTACGT", "ACGAACGT")
        assert m.mismatches == 1
        assert m.scored_positions == 8
        assert m.identity == pytest.approx(7 / 8)

    def test_none_params_means_default(self):
        a, b = "AAAA-TTGC", "AAA--TTCC"
        assert (score_alignment_metrics(a, b) ==
                score_alignment_metrics(a, b, DEFAULT_ADJUSTMENT_PARAMS))

    def test_unequal_lengths_raise_value_error_like_score_alignment(self):
        with pytest.raises(ValueError) as from_metrics:
            score_alignment_metrics("ACGT", "ACG")
        with pytest.raises(ValueError) as from_score:
            score_alignment("ACGT", "ACG")
        assert str(from_metrics.value) == str(from_score.value)

    def test_result_is_frozen(self):
        m = score_alignment_metrics("ACGT", "ACGT")
        with pytest.raises(FrozenInstanceError):
            m.identity = 0.5

    def test_deterministic_and_cacheable(self):
        """Results depend only on the strings and the (hashable) params."""
        a, b = "ACG-TTTAC", "ACGGTT-AC"
        params = AdjustmentParams(max_repeat_motif_length=1, hp_normalize_min_length=3)
        same_params = AdjustmentParams(max_repeat_motif_length=1, hp_normalize_min_length=3)
        cache = {(a, b, params): score_alignment_metrics(a, b, params)}
        assert (a, b, same_params) in cache
        assert cache[(a, b, same_params)] == score_alignment_metrics(a, b, same_params)
        hash(score_alignment_metrics(a, b, params))

    def test_does_not_build_visualization_or_variant_range_detail(self, monkeypatch):
        """The metrics path must not call any output-generation code."""
        def forbidden(*args, **kwargs):
            raise AssertionError("visualization code called from score_alignment_metrics")

        for name in ("_generate_annotated_output", "_generate_adjusted_output",
                     "_generate_annotated_score_string", "_generate_adjusted_score_string",
                     "_generate_variant_score_string", "_adjust_alignment_gaps",
                     "_extract_allele", "VariantRangeInfo", "AlignmentResult"):
            monkeypatch.setattr(adjusted_identity, name, forbidden)

        for a, b in SYNTHETIC_PAIRS[:10] + EDGE_CASE_PAIRS:
            score_alignment_metrics(a, b, DEFAULT_ADJUSTMENT_PARAMS)
            score_alignment_metrics(a, b, RAW_ADJUSTMENT_PARAMS)


class TestScoreAlignmentMetricsEquivalence:
    """Every field equals score_alignment's, for every parameter combination."""

    def test_synthetic_corpus_is_nontrivial(self):
        assert len(SYNTHETIC_PAIRS) >= 35
        default = [score_alignment_metrics(a, b) for a, b in SYNTHETIC_PAIRS]
        assert any(m.mismatches > 0 for m in default)
        assert any(m.identity == 1.0 for m in default)
        assert any("--" in a or "--" in b for a, b in SYNTHETIC_PAIRS)

    @pytest.mark.parametrize("params", PARAM_GRID, ids=_param_id)
    def test_synthetic_corpus_all_param_combinations(self, params):
        _assert_matches_score_alignment(SYNTHETIC_PAIRS, params)

    @pytest.mark.parametrize("params", PARAM_GRID, ids=_param_id)
    def test_edge_cases_all_param_combinations(self, params):
        _assert_matches_score_alignment(EDGE_CASE_PAIRS, params)

    @pytest.mark.skipif(not MSA_FASTA.exists(), reason=f"Local MSA not found: {MSA_FASTA}")
    @pytest.mark.parametrize("params", [
        GAPHACK_PARAMS,
        DEFAULT_ADJUSTMENT_PARAMS,
        RAW_ADJUSTMENT_PARAMS,
        AdjustmentParams(hp_normalize_min_length=6),
        AdjustmentParams(scoring_mode=ScoringMode.GLOBAL, end_skip_distance=20),
    ], ids=["gaphack", "default", "raw", "hp_min_6", "global_end_skip_20"])
    def test_msa_rows_all_pairs(self, params):
        """All pairs among the first rows of a real ONT MSA (local data, skipped in CI)."""
        from .compare_versions import read_fasta
        rows = [seq for _, seq in read_fasta(str(MSA_FASTA))][:25]
        _assert_matches_score_alignment(list(itertools.combinations(rows, 2)), params)
