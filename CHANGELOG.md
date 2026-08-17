# Changelog

## Version 0.2.9
- **Performance**: ~1.5x faster `align_and_score()` on real barcode data (validated byte-for-byte identical on 6,580 pinned alignments across 7 parameter configurations and both `adjust_gaps` modes)
  - Coverage and scored-position counting use C-level `str.count` on slices instead of per-character Python loops
  - Annotated score strings emit match segments in bulk, falling back to per-column classification only for segments containing gaps or ambiguity codes
  - `_are_nucleotides_equivalent` results are memoized (tiny alphabet, hot-loop call site)
  - `align_edlib_bidirectional` slices the original sequences instead of reverse-complementing back after RC-stage trimming (removes 2 of 4 RC passes)
  - `_analyze_allele` fast path for single-character alleles (the dominant case in divergent alignments); `_motif_matches` short-circuits on exact string equality
  - Context extraction walks once collecting up to `max_repeat_motif_length` characters instead of retrying at each shorter length (internal helper semantics: `_extract_left_context`/`_extract_right_context` now return the longest available context rather than exactly-N-or-None; callers' resulting contexts are unchanged)

## Version 0.2.8
- **Bug Fix**: Sequence regions beyond the overlap are no longer dropped from the alignment (LOCAL mode)
  - edlib's semi-global (HW) alignment omits target flanks outside the aligned location range; `align_edlib_bidirectional()` silently lost these regions, so `seq2_aligned` was truncated and `seq2_coverage` read 1.0 regardless of how much of seq2 was unaligned
  - Uncovered seq2 flanks are now re-attached with gap padding in seq1, making coverage symmetric under argument swap and consistent with `seq1_coverage`. Identity is unaffected (flanks lie outside the LOCAL scoring region); GLOBAL mode (NW alignment) was never affected
  - Note: workflows that filtered on `min(seq1_coverage, seq2_coverage)` will now correctly reject pairs where seq2 has large unaligned flanks
- **Bug Fix**: Completely dissimilar sequences no longer crash `align_and_score()`
  - When the RC-stage CIGAR reported the entire query as a terminal insertion, trimming emptied one sequence and the forward alignment raised an exception inside `edlib.getNiceAlignment`
  - Such pairs now return the existing alignment-failure sentinel (`identity=0.0`, `mismatches=-1`)
- **Bug Fix**: `adjust_gaps=True` no longer renders counted core indels with the extension marker
  - Every single-gap column was rendered as `=` (homopolymer extension, "normalized away") even when the variant range's core content was counted as an edit, contradicting the returned metrics
  - Core indel columns now show the same indel markers as annotated output (`indel_start` / `indel_extension`); extension and short-HP-demoted rendering is unchanged
- **Docs**: Corrected `end_skip_distance` docstrings to describe actual boundary semantics: scoring starts/ends at the position where both sequences reach their end_skip_distance-th nucleotide, i.e. the first/last `end_skip_distance - 1` nucleotides are excluded (behavior unchanged)

## Version 0.2.7
- **New Feature**: Added `hp_normalize_min_length` to `AdjustmentParams` (default `1`) — minimum homopolymer run length at which HP normalization applies
  - When the shorter side's HP run length (`min(L1, L2)`) is below the threshold, the length difference is counted as an edit instead of being normalized away
  - Applies only to true homopolymers (motif length 1); dinucleotide and longer repeat extensions are unaffected
  - Default `1` preserves existing behavior (all HP runs normalized). Set higher (e.g. `6`) to preserve short-HP signal where per-position HP error rates are comparable to the non-HP substitution baseline — see the companion [hp_error_rate_report](https://github.com/joshuaowalker/speconsense/tree/main/docs/hp_error_rate) for the empirical basis
  - Has no effect when `normalize_homopolymers=False`
- **New Field**: Added `short_hp_edit` to `ScoringFormat` (default `' '`) for the visualization marker at demoted short-HP positions; defaults to a space so existing visualizations match other scored mismatches. Override (e.g. `ScoringFormat(short_hp_edit='x')`) for a debug-visible marker
- **Internal**: `AlleleAnalysis` now records the motif length used for each extension (`left_motif_length`, `right_motif_length`), enabling the threshold check to target homopolymers specifically

## Version 0.2.6
- **New Feature**: Added `ScoringMode` enum (`LOCAL` / `GLOBAL`) to control whether terminal gaps are included in identity scoring
  - `LOCAL` (default): Scores only the overlap region, excluding terminal gaps. Uses semi-global (HW) alignment. This is the existing behavior — all callers are unaffected.
  - `GLOBAL`: Scores the full alignment including terminal gaps. Uses Needleman-Wunsch (NW) alignment. Solves the problem where a 220bp truncated sequence reports 100% identity against a 660bp full-length target because only the overlapping third is scored.
  - Set via `AdjustmentParams(scoring_mode=ScoringMode.GLOBAL)` or `AdjustmentParams(scoring_mode="global")` (string coercion supported)
- **New Function**: `_align_nw()` — Needleman-Wunsch alignment wrapper for GLOBAL mode
- **API Change**: `_find_scoring_region()` now accepts `adjustment_params` instead of `end_skip_distance` (private function, no public API impact)

## Version 0.2.5
- **New Feature**: Added `adjust_gaps` parameter to `score_alignment()` and `align_and_score()`
  - When `adjust_gaps=True`, gap positions are rewritten so the output alignment matches the scoring interpretation
  - Makes visualization directly interpretable position-by-position
  - Identity metrics are identical regardless of `adjust_gaps` setting
  - Defaults to `False` for full backward compatibility
- **Bug Fix**: Gap adjustment no longer applies extension-based rewriting when `normalize_homopolymers=False`
  - Previously, `adjust_gaps=True` could insert gaps to show homopolymer/repeat structure even when normalization was disabled
  - Now correctly treats all variant content as core (no extensions) when normalization is off
- **Architecture**: Unified single-pass analysis ensures both output modes produce identical metrics
- **Internal**: Removed deprecated `_score_alignment_impl` function
- **Internal**: Renamed private API functions for clarity

## Version 0.2.4
- Added CI workflow to run tests automatically on push and pull requests
- Fixed setuptools compatibility by removing legacy license classifier

## Version 0.2.3
- **Published to PyPI**: Package now available via `pip install adjusted-identity`
- Added GitHub Actions workflow for automated PyPI publishing with trusted publishing

## Version 0.2.2

- **Removed**: `score_aligned_seq2` field (added in v0.2.1) has been removed
  - Analysis showed it was redundant: same as `score_aligned` 98% of the time
  - Scoring is symmetric: swap seq1/seq2 arguments to get the alternate perspective
  - This simplifies the API and reduces memory overhead

## Version 0.2.1

- **Bug Fix**: Fixed dual-gap handling so they don't split variant ranges (key regression test added)
- **Bug Fix**: Fixed visualization when one position is extension and other is core with matching cores
- Improved visualization for indel normalization: first core position shows ` `, subsequent show `-`

## Version 0.2.0
- **Major Enhancement**: Implemented variant range algorithm for improved homopolymer and repeat motif detection
- **Key behavioral change**: Alternating indel patterns like `TGC-C-TC` vs `TGCT--TC` now correctly score as identity=1.0
  - The algorithm recognizes that C extends the left C context and T extends the right T context
  - Both alleles are valid repeat extensions → 0 edits (Occam's razor principle)
- **Algorithm improvements**:
  - Variant regions are now bounded by non-gap match positions (respects alignment boundaries)
  - Alleles extracted from variant ranges are analyzed for left/right repeat extensions
  - Split scoring: partial extensions allowed (e.g., "AAG" where "AA" extends context scores AA as 0 edits, G as 1 edit)
  - Opposite direction extensions are valid (allele1 extending left + allele2 extending right = both valid)
- **IUPAC integration**: Motif matching uses `_are_nucleotides_equivalent()` so IUPAC codes can extend context
- **Breaking change**: `end_skip_distance` now defaults to 0 (disabled). Set `end_skip_distance=20` to restore previous behavior.
- Removed 218 lines of dead code from previous indel processing implementation

## Version 0.1.7
- **Feature**: Added multi-sequence alignment (MSA) dual-gap support for homopolymer normalization
- Consensus-based context extraction now handles sequences where both have gaps at the same position (common in MSA outputs from spoa, MUSCLE, MAFFT)
- Dual-gap positions ('-' vs '-') are now correctly treated as matches, not indels
- Homopolymer detection uses consensus from both sequences when extracting context
- Added 17 comprehensive tests for MSA edge cases
- 100% backward compatible - all 133 tests pass
- No API changes - existing code works unchanged

## Version 0.1.6
- **Enhancement**: Added validation for contradictory `AdjustmentParams` configuration
- Now raises `ValueError` when `normalize_homopolymers=True` but `max_repeat_motif_length < 1` (which would silently disable homopolymer normalization)
- Added comprehensive test coverage for parameter validation edge cases
- No API changes - existing valid configurations work unchanged

## Version 0.1.5
- **Enhancement**: Added `ambiguous_match` field to `ScoringFormat` to distinguish between exact nucleotide matches and ambiguous matches
- Modified `_are_nucleotides_equivalent()` to return a tuple indicating match type
- Score patterns now show `|` for exact standard nucleotide matches (A=A, C=C, G=G, T=T) and `=` for any matches involving IUPAC ambiguity codes
- No breaking changes - existing code works unchanged but score visualization is more informative

## Version 0.1.4
- **Bug fix**: Fixed overhang scoring behavior when `end_skip_distance=0`
- Now correctly scores only positions where both sequences have content (no gap vs nucleotide scoring)
- Added comprehensive test suite for overhang region handling edge cases
- No API changes - existing code will work unchanged but may see different results for overhang alignments

## Version 0.1.3
- **Bug fix**: Fixed alignment length mismatch error in `align_edlib_bidirectional()`
- Resolved "Aligned sequences must have same length" errors for certain sequence pairs
- Simplified suffix trimming logic by removing unnecessary sequence trimming/reattachment
- No API changes or performance impact

## Version 0.1.2
- **Breaking**: Removed BioPython dependency - now only requires `edlib`
- Implemented custom `reverse_complement()` function with full IUPAC support
- Reduced package size and installation complexity
- Added comprehensive test coverage for reverse complement functionality
- Maintains 100% API compatibility (no code changes needed)

## Version 0.1.1
- Added repeat motif adjustment support (dinucleotide and longer repeats)
- Implemented intelligent motif length detection with degeneracy handling
- Added `max_repeat_motif_length` parameter to AdjustmentParams
- Enhanced left-right indel processing algorithm for mixed motif lengths
- Added comprehensive test coverage for repeat motif scenarios

## Version 0.1.0
- Initial release
- Complete MycoBLAST-style adjustment implementation (except repeat motifs)
- Comprehensive test suite
- Full documentation and examples
