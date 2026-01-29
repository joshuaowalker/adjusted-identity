# Changelog

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
