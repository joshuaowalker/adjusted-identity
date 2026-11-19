"""
Tests for the fast paths inside scoring.

Scoring classifies every alignment column once with C-level bytes operations,
caches the handling of short variant ranges, and uses string methods where the
straightforward code would loop over columns. Each shortcut has a plain
counterpart in the module (a column-by-column classifier, uncached functions, a
character walk) or an outside reference (edlib, the released behaviour). These
tests hold each shortcut to its counterpart, so a fast path cannot drift into
giving different output.
"""

import random
import re
from concurrent.futures import ThreadPoolExecutor
from dataclasses import astuple

import edlib
import pytest

import adjusted_identity
from adjusted_identity import (
    AdjustmentParams,
    AlleleAnalysis,
    DEFAULT_ADJUSTMENT_PARAMS,
    RAW_ADJUSTMENT_PARAMS,
    ScoringFormat,
    ScoringMode,
    _analyze_alignment,
    _are_nucleotides_equivalent,
    _classify_columns,
    _classify_columns_generic,
    _expand_alignment,
    _hp_run_length,
    _left_context,
    _parse_suffix_gap_from_cigar,
    _right_context,
    align_and_score,
    score_alignment,
    score_alignment_metrics,
)

from .test_score_alignment_metrics import EDGE_CASE_PAIRS, PARAM_GRID, SYNTHETIC_PAIRS, _param_id

IUPAC = "ACGTRYSWKMBDHVN"
FAST_ALPHABET = IUPAC + IUPAC.lower() + "-"

# Every eighth combination still covers each value of every parameter several times
PARAM_SAMPLE = PARAM_GRID[::8] + [
    DEFAULT_ADJUSTMENT_PARAMS,
    RAW_ADJUSTMENT_PARAMS,
    AdjustmentParams(end_skip_distance=0, max_repeat_motif_length=1, hp_normalize_min_length=6),
    AdjustmentParams(normalize_homopolymers=False, max_repeat_motif_length=0),
]


def _lowercase_some(rng, seq):
    return "".join(c.lower() if rng.random() < 0.3 else c for c in seq)


def _corpus():
    """Aligned pairs: the metrics corpus plus mixed-case and ambiguity-rich variants."""
    rng = random.Random(20261120)
    pairs = list(SYNTHETIC_PAIRS) + list(EDGE_CASE_PAIRS)
    for a, b in SYNTHETIC_PAIRS[:12]:
        pairs.append((_lowercase_some(rng, a), _lowercase_some(rng, b)))
        pairs.append(("".join(rng.choice("RYN") if c != "-" and rng.random() < 0.08 else c for c in a), b))
    return pairs


CORPUS = _corpus()

# The same corpus with characters the table-driven classifier cannot encode
ODD_CORPUS = [(a.replace("G", "X", 2).replace("T", "é", 1), b.replace("C", "u", 1))
              for a, b in SYNTHETIC_PAIRS[:15]]


def _all_outputs(a, b, params, scoring_format=None):
    """Everything scoring produces for one aligned pair."""
    return (
        astuple(score_alignment(a, b, params, scoring_format)),
        astuple(score_alignment(a, b, params, scoring_format, adjust_gaps=True)),
        astuple(score_alignment_metrics(a, b, params)),
        astuple(_analyze_alignment(a, b, params)),
    )


def _expected_class(char1, char2, handle_iupac):
    """Column class from first principles: gap rules, then nucleotide equivalence."""
    if char1 == "-" and char2 == "-":
        return "D"
    if char1 == "-":
        return "I"
    if char2 == "-":
        return "J"
    is_match, is_ambiguous = _are_nucleotides_equivalent(char1, char2, handle_iupac)
    if not is_match:
        return "X"
    if not is_ambiguous:
        return "M"
    return "N" if char1.upper() == char2.upper() else "Q"


class TestColumnClassifier:
    """_classify_columns labels columns exactly as nucleotide equivalence dictates."""

    @pytest.mark.parametrize("handle_iupac", [True, False])
    def test_every_symbol_pair(self, handle_iupac):
        """All 31 x 31 pairs of supported characters, through the table-driven path."""
        seq1 = "".join(c1 for c1 in FAST_ALPHABET for _ in FAST_ALPHABET)
        seq2 = FAST_ALPHABET * len(FAST_ALPHABET)
        expected = "".join(_expected_class(c1, c2, handle_iupac) for c1, c2 in zip(seq1, seq2))

        assert _classify_columns(seq1, seq2, handle_iupac).decode() == expected
        assert _classify_columns_generic(seq1, seq2, handle_iupac).decode() == expected

    @pytest.mark.parametrize("handle_iupac", [True, False])
    def test_unsupported_characters_fall_back(self, handle_iupac):
        """Characters outside the IUPAC alphabet are classified column by column."""
        alphabet = "ACGTNacgRy-" + "XxUu.*?1 " + "éÉßΣ"
        seq1 = "".join(c1 for c1 in alphabet for _ in alphabet)
        seq2 = alphabet * len(alphabet)
        expected = "".join(_expected_class(c1, c2, handle_iupac) for c1, c2 in zip(seq1, seq2))

        assert _classify_columns(seq1, seq2, handle_iupac).decode() == expected

    def test_specific_classes(self):
        assert _classify_columns("AaR-AX-A", "AAR-NX-C", True) == b"MMNDQNDX"
        # Without IUPAC intersection, two different ambiguity codes do not match
        assert _classify_columns("RN", "YR", True) == b"XQ"
        assert _classify_columns("RN", "YR", False) == b"XX"

    def test_str_subclass_and_plain_str_agree(self):
        class Sequence(str):
            pass

        for a, b in CORPUS[:10]:
            assert (_classify_columns(Sequence(a), Sequence(b), True)
                    == _classify_columns(a, b, True))

    def test_empty(self):
        assert _classify_columns("", "", True) == b""

    def test_random_strings_agree_with_generic(self):
        rng = random.Random(7)
        for _ in range(300):
            n = rng.randint(0, 200)
            a = "".join(rng.choice(FAST_ALPHABET) for _ in range(n))
            b = "".join(rng.choice(FAST_ALPHABET) for _ in range(n))
            for handle_iupac in (True, False):
                assert (_classify_columns(a, b, handle_iupac)
                        == _classify_columns_generic(a, b, handle_iupac))


class TestFastAndGenericScoringAgree:
    """The whole pipeline gives the same output on either classifier."""

    @pytest.mark.parametrize("params", PARAM_SAMPLE, ids=_param_id)
    def test_all_outputs(self, params, monkeypatch):
        fast = [_all_outputs(a, b, params) for a, b in CORPUS]
        monkeypatch.setattr(adjusted_identity, "_classify_columns", _classify_columns_generic)
        generic = [_all_outputs(a, b, params) for a, b in CORPUS]
        assert fast == generic


class TestCachesDoNotChangeResults:
    """Cached and uncached scoring give identical output."""

    CACHED = {
        "_score_variant": "_new_variant_score",
        "_motif_matches_cached": "_motif_matches",
        "_variant_codes": "_build_variant_codes",
        "_adjusted_variant": "_build_adjusted_variant",
    }

    def _uncached(self, monkeypatch):
        for cached_name, plain_name in self.CACHED.items():
            monkeypatch.setattr(adjusted_identity, cached_name,
                                getattr(adjusted_identity, plain_name))
        monkeypatch.setattr(adjusted_identity, "_allele_analysis", AlleleAnalysis)

    @pytest.mark.parametrize("params", PARAM_SAMPLE, ids=_param_id)
    def test_all_outputs(self, params, monkeypatch):
        pairs = CORPUS + ODD_CORPUS
        cold = [_all_outputs(a, b, params) for a, b in pairs]
        warm = [_all_outputs(a, b, params) for a, b in pairs]
        assert cold == warm

        self._uncached(monkeypatch)
        assert [_all_outputs(a, b, params) for a, b in pairs] == cold

    def test_caches_are_bounded(self):
        for name in list(self.CACHED) + ["_allele_analysis"]:
            assert getattr(adjusted_identity, name).cache_info().maxsize is not None

    def test_long_variant_ranges_skip_the_cache(self):
        """A long indel is scored directly rather than stored."""
        a = "ACGTACGT" + "GATTACAGATTACAGATTACA" + "TTGACCA"
        b = "ACGTACGT" + "-" * 21 + "TTGACCA"
        adjusted_identity._score_variant.cache_clear()
        result = score_alignment(a, b)
        assert result.mismatches == 1
        assert adjusted_identity._score_variant.cache_info().currsize == 0

    def test_results_do_not_share_mutable_state(self):
        """Mutating one result's score dict must not leak into later results."""
        a, b = "AAAA-TTGC", "AAA--TTCC"
        first = _analyze_alignment(a, b)
        expected = [dict(vr.score_result) for vr in first.variant_ranges]
        for vr in first.variant_ranges:
            vr.score_result["edits"] = 99
            vr.score_result.clear()
        second = _analyze_alignment(a, b)
        assert [vr.score_result for vr in second.variant_ranges] == expected

    def test_params_that_differ_only_in_type_are_not_confused(self):
        """Truthy non-bool flags and float lengths score like their plain equivalents."""
        class Params:
            normalize_homopolymers = 1
            handle_iupac_overlap = 1
            normalize_indels = 1
            end_skip_distance = 0
            max_repeat_motif_length = 2
            scoring_mode = ScoringMode.LOCAL
            hp_normalize_min_length = 1

        for a, b in CORPUS:
            assert (astuple(score_alignment(a, b, Params()))
                    == astuple(score_alignment(a, b, DEFAULT_ADJUSTMENT_PARAMS)))

    def test_threads(self):
        """Concurrent scoring gives the same results as serial scoring."""
        jobs = [(a, b, p) for p in PARAM_SAMPLE[:6] for a, b in CORPUS]
        serial = [_all_outputs(*job) for job in jobs]
        for name in self.CACHED:
            getattr(adjusted_identity, name).cache_clear()
        with ThreadPoolExecutor(max_workers=8) as pool:
            threaded = list(pool.map(lambda job: _all_outputs(*job), jobs))
        assert threaded == serial


def _strip_dual_gap_columns(a, b):
    kept = [(x, y) for x, y in zip(a, b) if not (x == "-" and y == "-")]
    return "".join(x for x, _ in kept), "".join(y for _, y in kept)


class TestDualGapColumns:
    """Columns that are gaps in both sequences never affect the metrics."""

    DUAL_GAP_PAIRS = [(a, b) for a, b in CORPUS
                      if any(x == "-" and y == "-" for x, y in zip(a, b))]

    def test_corpus_has_dual_gaps(self):
        assert len(self.DUAL_GAP_PAIRS) >= 10

    @pytest.mark.parametrize("params", PARAM_SAMPLE, ids=_param_id)
    def test_metrics_unchanged_by_removing_them(self, params):
        for a, b in self.DUAL_GAP_PAIRS:
            stripped = _strip_dual_gap_columns(a, b)
            assert (score_alignment_metrics(a, b, params)
                    == score_alignment_metrics(*stripped, params)), (a, b)

    def test_analysis_still_reports_dual_gap_ranges(self):
        analysis = _analyze_alignment("AC--GT-A", "AC--GTTA")
        spans = [(vr.start, vr.end, vr.allele1, vr.allele2) for vr in analysis.variant_ranges]
        assert spans == [(2, 3, "", ""), (6, 6, "", "T")]
        assert analysis.variant_ranges[0].score_result == {
            "edits": 0, "scored_positions": 0, "both_pure_extension": True, "short_hp_edit": False}

    def test_dual_gaps_inside_a_variant_range(self):
        """Dual gaps between two differing columns join them into one range."""
        with_gap = score_alignment("ACA--GT", "AC-T-GT")
        without = score_alignment("ACA-GT", "AC-TGT")
        assert (with_gap.mismatches, with_gap.scored_positions) == (
            without.mismatches, without.scored_positions)
        assert [(vr.start, vr.end) for vr in
                _analyze_alignment("ACA--GT", "AC-T-GT").variant_ranges] == [(2, 4)]

    def test_score_strings(self):
        annotated = score_alignment("AC--GT-A", "AC--GTTA")
        assert annotated.score_aligned == "||..||=|"
        adjusted = score_alignment("AC--GT-A", "AC--GTTA", adjust_gaps=True)
        assert (adjusted.seq1_aligned, adjusted.seq2_aligned, adjusted.score_aligned) == (
            "ACGT-A", "ACGTTA", "||||=|")


class TestHomopolymerRunLength:
    """The str-method measurement of a homopolymer run equals the column walk."""

    def test_random_sequences(self):
        rng = random.Random(11)
        for _ in range(3000):
            parts = []
            for _ in range(rng.randint(1, 5)):
                parts.append(rng.choice("AaCT") * rng.choice([1, 2, 5, 31, 32, 33, 34, 70]))
                parts.append("-" * rng.choice([0, 0, 1, 3]))
            seq = "".join(c if rng.random() > 0.05 else "-" for c in "".join(parts))
            start = rng.randrange(len(seq))
            end = min(len(seq) - 1, start + rng.choice([0, 0, 1, 2, 5]))
            for base in "AaC":
                run_chars = base.upper() + base.lower()
                assert (_hp_run_length(seq, start, end, base, run_chars)
                        == _hp_run_length(seq, start, end, base)), (seq, start, end, base)

    def test_run_longer_than_the_window(self):
        seq = "G" + "A" * 100 + "-" + "A" * 100 + "C"
        assert _hp_run_length(seq, 101, 101, "A", "Aa") == 200
        assert _hp_run_length(seq, 101, 101, "A") == 200

    def test_threshold_beyond_the_window(self):
        """A run just short of a large threshold is demoted; just at it, normalized."""
        params = AdjustmentParams(hp_normalize_min_length=60)
        short = score_alignment("C" + "A" * 59 + "-T", "C" + "A" * 60 + "T", params)
        long_ = score_alignment("C" + "A" * 60 + "-T", "C" + "A" * 61 + "T", params)
        assert (short.mismatches, long_.mismatches) == (1, 0)


class TestRepeatContext:
    """Context is read from columns where both sequences hold the same symbol."""

    def _classes(self, a, b):
        return _classify_columns(a, b, True)

    def test_left(self):
        assert _left_context(self._classes("AGG-AC", "AG-GAC"), "AGG-AC", 2, 1) == "G"
        assert _left_context(self._classes("AGT-AC", "AX-GAC"), "AGT-AC", 2, 1) is None
        assert _left_context(self._classes("ACGT", "ACGT"), "ACGT", 4, 2) == "GT"
        assert _left_context(self._classes("ACGT", "ACGT"), "ACGT", 1, 3) == "A"
        assert _left_context(self._classes("ACGT", "ACGT"), "ACGT", 0, 2) is None

    def test_right(self):
        assert _right_context(self._classes("AGA--TT", "AGAT-TT"), "AGA--TT", 4, 2) == "TT"
        assert _right_context(self._classes("AGA--TC", "AGAT-TG"), "AGA--TC", 4, 2) == "T"
        assert _right_context(self._classes("ACGT", "ACGT"), "ACGT", 3, 2) is None

    def test_dual_gaps_are_skipped(self):
        a, b = "AC--G-T", "AC--GTT"
        assert _left_context(self._classes(a, b), a, 5, 3) == "ACG"
        a, b = "T-G--CA", "TTG--CA"
        assert _right_context(self._classes(a, b), a, 1, 3) == "GCA"

    def test_equivalent_but_different_symbols_stop_the_walk(self):
        a, b = "ARG-T", "AAGGT"
        # R/A match for scoring but are not the same symbol, so context stops there
        assert _left_context(self._classes(a, b), a, 3, 3) == "G"
        a, b = "ARG-T", "ARGGT"
        assert _left_context(self._classes(a, b), a, 3, 3) == "ARG"

    def test_zero_length(self):
        assert _left_context(self._classes("ACGT", "ACGT"), "ACGT", 2, 0) is None
        assert _right_context(self._classes("ACGT", "ACGT"), "ACGT", 1, 0) is None


def _random_pair(rng):
    alphabet = rng.choice(["ACGT", "ACGTN", "AC", "ACGTRYacgt"])
    a = "".join(rng.choice(alphabet) for _ in range(rng.randint(1, 80)))
    if rng.random() < 0.3:
        return a, "".join(rng.choice(alphabet) for _ in range(rng.randint(1, 80)))
    b = list(a)
    for _ in range(rng.randint(0, 8)):
        i = rng.randrange(len(b) + 1)
        r = rng.random()
        if r < 0.4 and b:
            b[min(i, len(b) - 1)] = rng.choice(alphabet)
        elif r < 0.7:
            b.insert(i, rng.choice(alphabet) * rng.randint(1, 4))
        elif b:
            del b[min(i, len(b) - 1)]
    b = "".join(b) or "A"
    if rng.random() < 0.3:
        b = "".join(rng.choice(alphabet) for _ in range(rng.randint(1, 9))) + b
    if rng.random() < 0.3:
        a = "".join(rng.choice(alphabet) for _ in range(rng.randint(1, 9))) + a
    return a, b


class TestExpandAlignment:
    """_expand_alignment builds the same strings as edlib.getNiceAlignment."""

    @pytest.mark.parametrize("mode", ["HW", "NW"])
    def test_matches_edlib(self, mode):
        rng = random.Random(5)
        for _ in range(1500):
            query, target = _random_pair(rng)
            result = edlib.align(query, target, mode=mode, task="path")
            nice = edlib.getNiceAlignment(result, query, target)
            assert _expand_alignment(result, query, target) == (
                nice["query_aligned"], nice["target_aligned"]), (query, target, result)

    def test_query_overhanging_target(self):
        result = edlib.align("AAAAGGGG", "AAAA", mode="HW", task="path")
        assert _expand_alignment(result, "AAAAGGGG", "AAAA") == ("AAAAGGGG", "AAAA----")

    def test_unexpected_cigar_is_left_to_edlib(self):
        """Anything not strictly well-formed gets edlib's own handling, errors included."""
        junk = {"cigar": "3=2?", "locations": [(0, 4)]}
        with pytest.raises(Exception) as ours:
            _expand_alignment(junk, "ACGTA", "ACGTA")
        with pytest.raises(Exception) as theirs:
            edlib.getNiceAlignment(junk, "ACGTA", "ACGTA")
        assert str(ours.value) == str(theirs.value)

        for cigar in ("3=junk2=", "=3", "3"):
            loose = {"cigar": cigar, "locations": [(0, 4)]}
            try:
                nice = edlib.getNiceAlignment(loose, "ACGTA", "ACGTA")
            except Exception as error:
                with pytest.raises(Exception) as ours:
                    _expand_alignment(loose, "ACGTA", "ACGTA")
                assert str(ours.value) == str(error)
            else:
                assert _expand_alignment(loose, "ACGTA", "ACGTA") == (
                    nice["query_aligned"], nice["target_aligned"])


def _suffix_gap_reference(cigar_string):
    """The original implementation: tokenize the whole CIGAR, then scan from the end."""
    if not cigar_string:
        return None
    cigar_ops = re.findall(r'(\d+)([=XIDM])', cigar_string)
    if not cigar_ops:
        return None
    gap_length = 0
    gap_type = None
    for length_str, op in reversed(cigar_ops):
        if op in ('D', 'I'):
            if gap_type is None:
                gap_type = op
            elif gap_type != op:
                gap_length = 0
                break
            gap_length += int(length_str)
        else:
            break
    if gap_length > 0:
        return (gap_length, gap_type == 'D')
    return None


class TestSuffixGapParsing:
    """Reading only the tail of the CIGAR gives the same answer as tokenizing all of it."""

    def test_real_cigars(self):
        rng = random.Random(3)
        seen = set()
        for _ in range(1500):
            query, target = _random_pair(rng)
            cigar = edlib.align(query, target, mode="HW", task="path")["cigar"]
            result = _parse_suffix_gap_from_cigar(cigar)
            assert result == _suffix_gap_reference(cigar), cigar
            seen.add(result is None)
        assert seen == {True, False}

    def test_arbitrary_strings(self):
        rng = random.Random(4)
        for _ in range(5000):
            text = "".join(rng.choice("0123456789=XIDM?- ") for _ in range(rng.randint(0, 14)))
            assert _parse_suffix_gap_from_cigar(text) == _suffix_gap_reference(text), text


# Output of adjusted-identity 0.2.10 for inputs the table-driven classifier
# does not handle: (seq1, seq2, params, (mismatches, scored_positions,
# seq1_coverage, seq2_coverage), score_aligned, adjust_gaps=True output)
RELEASED_OUTPUT = [
    ('ACXTAC-GTT', 'ACXTACAGTT', 'default', (1, 10, 1.0, 1.0), '||=||| |||',
     ('ACXTAC-GTT', 'ACXTACAGTT', '||=||| |||')),
    ('ACXTACGT', 'ACNTACGT', 'default', (1, 8, 1.0, 1.0), '|| |||||',
     ('AC-XTACGT', 'ACN-TACGT', '||= |||||')),
    ('ACXTACGT', 'ACNTACGT', 'raw', (1, 8, 1.0, 1.0), '|| |||||',
     ('ACXTACGT', 'ACNTACGT', '|| |||||')),
    ('ACXTACGT', 'ACNTACGT', 'global_no_iupac', (1, 8, 1.0, 1.0), '|| |||||',
     ('AC-XTACGT', 'ACN-TACGT', '||  |||||')),
    ('acgtRYn-AC', 'ACGTYRNAAC', 'default', (0, 7, 1.0, 1.0), '||||====||',
     ('acgtRY-n-AC', 'ACGTY-RNAAC', '|||| ====||')),
    ('acgtRYn-AC', 'ACGTYRNAAC', 'raw', (3, 10, 1.0, 1.0), '||||  = ||',
     ('acgtRYn-AC', 'ACGTYRNAAC', '||||  = ||')),
    ('acgtRYn-AC', 'ACGTYRNAAC', 'global_no_iupac', (2, 9, 1.0, 1.0), '||||  = ||',
     ('acgtRY-n-AC', 'ACGTY-RNAAC', '||||   = ||')),
    ('AC.GT*TT-A', 'ACNGT*TTTA', 'default', (1, 9, 1.0, 1.0), '|| ||=||=|',
     ('AC-.GT*TT-A', 'ACN-GT*TTTA', '||= ||=||=|')),
    ('AC.GT*TT-A', 'ACNGT*TTTA', 'raw', (2, 10, 1.0, 1.0), '|| ||=|| |',
     ('AC.GT*TT-A', 'ACNGT*TTTA', '|| ||=|| |')),
    ('AC.GT*TT-A', 'ACNGT*TTTA', 'global_no_iupac', (2, 10, 1.0, 1.0), '|| ||=|| |',
     ('AC-.GT*TT-A', 'ACN-GT*TTTA', '||  ||=|| |')),
    ('AUGCAUUG-C', 'ATGCATTGGC', 'default', (3, 9, 1.0, 1.0), '| |||  |=|',
     ('AUGCAUUG-C', 'ATGCATTGGC', '| |||  |=|')),
    ('AUGCAUUG-C', 'ATGCATTGGC', 'raw', (4, 10, 1.0, 1.0), '| |||  | |',
     ('AUGCAUUG-C', 'ATGCATTGGC', '| |||  | |')),
    ('ACG\xc9T-AA', 'ACG\xe9TTAA', 'default', (0, 7, 1.0, 1.0), '|||=|=||',
     ('ACG\xc9T-AA', 'ACG\xe9TTAA', '|||=|=||')),
    ('ACG\xc9T-AA', 'ACG\xe9TTAA', 'raw', (1, 8, 1.0, 1.0), '|||=| ||',
     ('ACG\xc9T-AA', 'ACG\xe9TTAA', '|||=| ||')),
    ('xx-ACGT', 'XXXACGT', 'default', (0, 6, 1.0, 1.0), '===||||',
     ('xx-ACGT', 'XXXACGT', '===||||')),
    ('xx-ACGT', 'XXXACGT', 'raw', (1, 7, 1.0, 1.0), '== ||||',
     ('xx-ACGT', 'XXXACGT', '== ||||')),
    ('--AC?T--', 'TTAC?TGG', 'default', (0, 4, 1.0, 0.5), '..||=|..',
     ('--AC?T--', 'TTAC?TGG', '==||=|==')),
    ('--AC?T--', 'TTAC?TGG', 'raw', (0, 4, 1.0, 0.5), '..||=|..',
     ('--AC?T--', 'TTAC?TGG', '--||=|--')),
    ('--AC?T--', 'TTAC?TGG', 'global_no_iupac', (2, 6, 1.0, 0.5), ' -||=| -',
     ('--AC?T--', 'TTAC?TGG', ' -||=| -')),
]

RELEASED_PARAMS = {
    "default": DEFAULT_ADJUSTMENT_PARAMS,
    "raw": RAW_ADJUSTMENT_PARAMS,
    "global_no_iupac": AdjustmentParams(scoring_mode="global", handle_iupac_overlap=False,
                                        hp_normalize_min_length=3),
}


class TestUnusualCharacters:
    """Lower case, non-IUPAC and non-ASCII input scores as it did in 0.2.10."""

    @pytest.mark.parametrize("case", RELEASED_OUTPUT,
                             ids=lambda c: f"{c[0].encode('ascii', 'replace').decode()}-{c[2]}")
    def test_released_output(self, case):
        seq1, seq2, params_name, metrics, score_aligned, adjusted = case
        params = RELEASED_PARAMS[params_name]

        result = score_alignment(seq1, seq2, params)
        assert (result.mismatches, result.scored_positions,
                result.seq1_coverage, result.seq2_coverage) == metrics
        assert result.score_aligned == score_aligned
        assert (result.seq1_aligned, result.seq2_aligned) == (seq1, seq2)

        result = score_alignment(seq1, seq2, params, adjust_gaps=True)
        assert (result.seq1_aligned, result.seq2_aligned, result.score_aligned) == adjusted

        fast_metrics = score_alignment_metrics(seq1, seq2, params)
        assert (fast_metrics.mismatches, fast_metrics.scored_positions,
                fast_metrics.seq1_coverage, fast_metrics.seq2_coverage) == metrics


# More 0.2.10 output, for short alignments chosen because each one tells the
# real code apart from a specific plausible mistake (named in the comment):
# (seq1, seq2, non-default params, (mismatches, scored_positions),
#  score_aligned, adjust_gaps=True output)
BOUNDARY_CASES = [
    # terminal gaps scored in GLOBAL mode; short-HP demotion of a multi-base range
    ('---N-NGCTA', 'TGCNRNTC-A',
     {'max_repeat_motif_length': 1, 'scoring_mode': 'global', 'hp_normalize_min_length': 3},
     (4, 8), '   = = | |', ('---N-NGCTA', 'TGCNRNTC-A', '   = = | |')),
    # gap columns outside the scoring region, drawn in the adjusted alignment
    ('-C-A', 'AT--', {'end_skip_distance': 2, 'max_repeat_motif_length': 3},
     (1, 1), '. ..', ('-C-A', 'AT--', '= .=')),
    # short-HP demotion decided by the run in seq2
    ('-AA-AAAAATAAAT', 'TAAAAAAAATA---', {'scoring_mode': 'global', 'hp_normalize_min_length': 3},
     (2, 11), ' ||=|||||||   ', ('-AA-AAAAATAAAT', 'TAAAAAAAATA---', ' ||=|||||||   ')),
    # short-HP demotion decided by the run in seq1
    ('AA--CCAC', 'AAACCCAC',
     {'handle_iupac_overlap': False, 'end_skip_distance': 2, 'max_repeat_motif_length': 1,
      'hp_normalize_min_length': 3},
     (1, 5), '.|  |||.', ('AA--CCAC', 'AAACCCAC', '||  ||||')),
    # end trimming with end_skip_distance=1
    ('AA--T', 'AAA--',
     {'end_skip_distance': 1, 'max_repeat_motif_length': 1, 'scoring_mode': 'global',
      'hp_normalize_min_length': 2},
     (0, 2), '||=..', ('AA--T', 'AAA--', '||=.=')),
    # homopolymer run counted across upper and lower case
    ('Aa-AA', '--aAA',
     {'handle_iupac_overlap': False, 'normalize_indels': False, 'scoring_mode': 'global',
      'hp_normalize_min_length': 3},
     (0, 2), '===||', ('AaAA', '-aAA', '=|||')),
    # end trimming leaves nothing to score (the annotated string is one longer
    # than the alignment here, as it has been since end trimming was added)
    ('AC-AC', 'ACAAC',
     {'end_skip_distance': 3, 'max_repeat_motif_length': 1, 'scoring_mode': 'global',
      'hp_normalize_min_length': 4},
     (0, 0), '......', ('AC-AC', 'ACAAC', '||=||')),
    # a variant range at the start of the scoring region has no left context,
    # even though matching columns precede it outside the region
    ('ATAAATATAAAA-A', 'ATTAATATAAAA-A',
     {'handle_iupac_overlap': False, 'end_skip_distance': 3, 'max_repeat_motif_length': 3,
      'scoring_mode': 'global'},
     (1, 9), '.. ||||||||...', ('AT-AAATATAAAA-A', 'ATT-AATATAAAA-A', '|| =|||||||||.|')),
    # ... and one at the end of the region has no right context
    ('-aCACACC-', 'aaCACA-CC',
     {'handle_iupac_overlap': False, 'end_skip_distance': 2, 'max_repeat_motif_length': 1,
      'scoring_mode': 'global', 'hp_normalize_min_length': 2},
     (1, 5), '..|||| ..', ('-aCACACC-', 'aaCACA-CC', '=||||| |=')),
    # the same, for a multi-base variant range
    ('TC--GT-RTNT', 'TTR-GT-RTNT',
     {'handle_iupac_overlap': False, 'end_skip_distance': 2, 'max_repeat_motif_length': 1,
      'hp_normalize_min_length': 2},
     (1, 6), '. =.||.=|=.', ('TC-GTRTNT', 'TTRGTRTNT', '| =||=|=|')),
    ('ACaaCaaaaCC', 'ACaaCaaa-AC',
     {'normalize_homopolymers': False, 'normalize_indels': False, 'end_skip_distance': 2,
      'max_repeat_motif_length': 1, 'scoring_mode': 'global', 'hp_normalize_min_length': 3},
     (2, 9), '.|||||||  .', ('ACaaCaaaaCC', 'ACaaCaaaA-C', '||||||||| |')),
    # ambiguous matches inside an adjusted variant range
    ('---GG', 'N-N--',
     {'handle_iupac_overlap': False, 'end_skip_distance': 4, 'max_repeat_motif_length': 3},
     (0, 2), ' .---', ('GG', 'NN', '==')),
    # a scoring region with no match column at all
    ('ACAAA', 'ACCAC', {'end_skip_distance': 3, 'hp_normalize_min_length': 3},
     (1, 1), '.. ..', ('ACAAA', 'ACCAC', '|| | ')),
    # motif comparison without IUPAC intersection
    ('GNC--', 'GNRAT',
     {'normalize_homopolymers': False, 'handle_iupac_overlap': False,
      'max_repeat_motif_length': 1, 'scoring_mode': 'global', 'hp_normalize_min_length': 3},
     (1, 3), '|=  -', ('GNC--', 'GNRAT', '|=  -')),
]


class TestBoundaryCases:
    """Scoring-region edges and short-HP demotion score as they did in 0.2.10."""

    @pytest.mark.parametrize("case", BOUNDARY_CASES, ids=lambda c: f"{c[0]}/{c[1]}")
    def test_released_output(self, case):
        seq1, seq2, params_kwargs, metrics, score_aligned, adjusted = case
        params = AdjustmentParams(**params_kwargs)

        # Warm and cold caches alike
        for _ in range(2):
            result = score_alignment(seq1, seq2, params)
            assert (result.mismatches, result.scored_positions) == metrics
            assert result.score_aligned == score_aligned

            result = score_alignment(seq1, seq2, params, adjust_gaps=True)
            assert (result.seq1_aligned, result.seq2_aligned, result.score_aligned) == adjusted

            fast_metrics = score_alignment_metrics(seq1, seq2, params)
            assert (fast_metrics.mismatches, fast_metrics.scored_positions) == metrics


class TestScoringFormats:
    """Markers are applied the same way whatever characters they are."""

    ASCII = ScoringFormat(match='m', ambiguous_match='a', substitution='s', indel_start='i',
                          indel_extension='e', homopolymer_extension='h', short_hp_edit='p',
                          end_trimmed='t', dual_gap='d')
    # Outside Latin-1, so translated through a different table type
    WIDE = ScoringFormat(match='│', ambiguous_match='≈', substitution='×',
                         indel_start='▶', indel_extension='─',
                         homopolymer_extension='≡', short_hp_edit='◇',
                         end_trimmed='…', dual_gap='░')

    @pytest.mark.parametrize("adjust_gaps", [False, True])
    def test_wide_markers_match_ascii_markers(self, adjust_gaps):
        to_wide = str.maketrans({
            getattr(self.ASCII, name): getattr(self.WIDE, name)
            for name in ("match", "ambiguous_match", "substitution", "indel_start",
                         "indel_extension", "homopolymer_extension", "short_hp_edit",
                         "end_trimmed", "dual_gap")})
        codes_seen = set()
        for params in PARAM_SAMPLE:
            for a, b in CORPUS:
                plain = score_alignment(a, b, params, self.ASCII, adjust_gaps)
                wide = score_alignment(a, b, params, self.WIDE, adjust_gaps)
                assert wide.score_aligned == plain.score_aligned.translate(to_wide)
                assert astuple(wide)[:7] == astuple(plain)[:7]
                codes_seen.update(plain.score_aligned)
        # Every kind of marker appears somewhere in the corpus
        expected = set("masiehptd") if not adjust_gaps else set("masiehpd")
        assert codes_seen == expected

    def test_alternating_formats(self):
        """Switching formats between calls never reuses the other format's markers."""
        a, b = CORPUS[0]
        default = score_alignment(a, b).score_aligned
        custom = score_alignment(a, b, scoring_format=self.ASCII).score_aligned
        for _ in range(3):
            assert score_alignment(a, b).score_aligned == default
            assert score_alignment(a, b, scoring_format=self.ASCII).score_aligned == custom
            assert score_alignment(a, b, scoring_format=ScoringFormat()).score_aligned == default


class TestEndToEnd:
    """align_and_score is deterministic across repeated and reordered calls."""

    def test_repeat_and_reorder(self):
        rng = random.Random(9)
        pairs = [_random_pair(rng) for _ in range(120)]
        first = [astuple(align_and_score(a, b)) for a, b in pairs]
        shuffled = list(range(len(pairs)))
        rng.shuffle(shuffled)
        again = {i: astuple(align_and_score(*pairs[i])) for i in shuffled}
        assert [again[i] for i in range(len(pairs))] == first
