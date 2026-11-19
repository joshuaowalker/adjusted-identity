# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Python package that implements MycoBLAST-style sequence identity calculations for DNA sequences, specifically designed for mycological DNA barcoding applications. The package provides functions to calculate sequence identity metrics that account for homopolymer length differences and other technical sequencing artifacts.

## Package Structure

```
adjusted-identity/
├── adjusted_identity/          # Main package directory
│   └── __init__.py            # Main module (moved from root)
├── tests/                     # Test suite
│   ├── __init__.py
│   ├── test_score_alignment.py  # Core scoring tests
│   ├── test_fast_paths.py      # Fast paths vs plain counterparts
│   ├── test_alignment.py       # Alignment algorithm tests
│   └── test_utilities.py       # Utility function tests
├── setup.py                   # Setup configuration
├── pyproject.toml             # Modern Python packaging
├── MANIFEST.in               # Package manifest
├── README.md                 # Package documentation
├── LICENSE                   # BSD 2-clause license
├── .gitignore               # Git ignore rules
└── CLAUDE.md               # This file
```

## Development Commands

### Installation
```bash
# Install from PyPI
pip install adjusted-identity

# Development installation with test dependencies
pip install -e ".[dev]"

# Install from GitHub (latest development version)
pip install git+https://github.com/joshuaowalker/adjusted-identity.git
```

### Testing
```bash
# Run all tests
pytest

# Run tests with coverage
pytest --cov=adjusted_identity --cov-report=html

# Run specific test file
pytest tests/test_score_alignment.py

# Run specific test class
pytest tests/test_score_alignment.py::TestBasicMatching

# Run specific test
pytest tests/test_score_alignment.py::TestBasicMatching::test_perfect_match
```

### Package Building
```bash
# Build source distribution and wheel
python -m build

# Check package metadata
python setup.py check

# Upload to PyPI (uses trusted publishing via GitHub Actions)
twine upload dist/*
```

## Core Architecture

### Main Components

1. **Data Classes & Enums** (`adjusted_identity/__init__.py`):
   - `ScoringMode` - Enum: `LOCAL` (overlap only) or `GLOBAL` (full alignment with terminal gaps)
   - `AlignmentResult` - Contains alignment results and identity metrics
   - `AlignmentMetrics` - The five metric fields of `AlignmentResult`, without the strings
   - `ScoringFormat` - Format codes for alignment visualization
   - `AdjustmentParams` - Parameters controlling sequence adjustments (includes `scoring_mode`)

2. **Core Functions** (`adjusted_identity/__init__.py`):
   - `align_and_score()` - Main entry point for sequence comparison (routes to HW or NW alignment based on scoring mode)
   - `align_edlib_bidirectional()` - Multi-stage semi-global alignment optimization (used by LOCAL mode)
   - `_align_nw()` - Needleman-Wunsch global alignment wrapper (used by GLOBAL mode)
   - `score_alignment()` - Scoring with configurable adjustments
   - `score_alignment_metrics()` - Metrics-only scoring (returns `AlignmentMetrics`, no display strings); public API used by gaphack. Must stay field-for-field equal to `score_alignment()`

3. **Adjustment Features**:
   - Homopolymer length normalization - ignores differences in homopolymer run lengths
   - Short-HP threshold (`hp_normalize_min_length`) - opt in to counting length diffs in short HP runs as real edits (default 1 = normalize all)
   - IUPAC ambiguity code handling - allows different ambiguity codes to match
   - End trimming - skips mismatches in terminal regions (disabled by default, set `end_skip_distance` to enable)
   - Indel normalization - counts contiguous indels as single events
   - LOCAL/GLOBAL scoring modes - controls whether terminal gaps affect identity

### How scoring works internally

Scoring is pure Python built on C-level bytes/str methods, with no per-column Python loop on normal input:

1. `_classify_columns()` labels every column with one class byte: `M` same standard nucleotide, `N` same non-standard symbol, `Q` different but equivalent symbols, `X` mismatch, `I`/`J` gap in seq1/seq2 only, `D` gap in both. The table-driven path handles plain `str` input in the IUPAC alphabet; anything else goes to `_classify_columns_generic()`, a column-by-column equivalent. `_are_nucleotides_equivalent()` remains the single definition of equivalence (the tables are built from it).
2. `_compute_alignment_metrics()` (the single source of truth for metrics) reads the overlap, scoring region, variant ranges, repeat context and match count from the classes. Variant ranges are maximal runs of non-match columns; dual-gap columns never affect metrics and are skipped unless `include_dual_gap_ranges=True` (used by `_analyze_alignment`, which reports them).
3. Each variant range is scored by `_score_variant()` (cached) from its two alleles and the repeat context; the short-HP threshold is applied on top because it depends on the surrounding sequence.
4. Output is assembled as format-independent "score codes" (`_generate_annotated_score_string`, `_generate_adjusted_output`) and turned into `ScoringFormat` markers in one translate (`_codes_to_markers`).

Caches (`_score_variant`, `_variant_codes`, `_adjusted_variant`, `_motif_matches_cached`, `_allele_analysis`) are bounded LRUs keyed on plain values. Results must never depend on cache state; `tests/test_fast_paths.py` enforces this.

**Output must stay identical when optimizing.** Check with the local `perf_regression/check.py` harness and a differential run against the previous release's module (`git show <tag>:adjusted_identity/__init__.py`), not only the test suite.

### Key Constants

- `DEFAULT_ADJUSTMENT_PARAMS` - All adjustments enabled (typical use case)
- `RAW_ADJUSTMENT_PARAMS` - No adjustments (traditional sequence identity)
- `IUPAC_CODES` - Nucleotide ambiguity code definitions

## Test Suite Organization

The test suite is comprehensive and serves as documentation:

1. **test_score_alignment.py**: Main scoring function tests
   - `TestBasicMatching` - Perfect matches and substitutions
   - `TestIndelScoring` - Insertion/deletion handling with/without normalization
   - `TestHomopolymerAdjustment` - Homopolymer length differences
   - `TestIUPACAdjustment` - IUPAC ambiguity code handling
   - `TestEndTrimming` - Terminal region mismatch skipping
   - `TestCombinedAdjustments` - Multiple adjustments together
   - `TestEdgeCases` - Error conditions and edge cases
   - `TestDocumentationExamples` - Real-world usage examples
   - `TestScoringMode` - LOCAL vs GLOBAL scoring mode comparisons
   - `TestShortHomopolymerThreshold` - `hp_normalize_min_length` threshold behavior

2. **test_alignment.py**: Alignment algorithm tests
   - `TestAlignEdlibBidirectional` - Bidirectional alignment
   - `TestAlignAndScore` - End-to-end functionality
   - `TestRealWorldScenarios` - Mycological sequence examples
   - `TestNWAlignment` - Needleman-Wunsch alignment wrapper
   - `TestGlobalModeEndToEnd` - GLOBAL scoring mode end-to-end tests

3. **test_score_alignment_metrics.py**: `score_alignment_metrics()` contract tests
   - Field-for-field equality with `score_alignment()` (both `adjust_gaps` modes) over a seeded synthetic corpus × every AdjustmentParams combination, plus the local MSA when present
   - Asserts the metrics path never calls visualization or variant-range-detail code

4. **test_fast_paths.py**: holds every fast path to a plain counterpart
   - Table-driven column classifier vs the column-by-column one, exhaustively over symbol pairs and through the whole pipeline
   - Cached vs uncached scoring; thread safety; dual-gap columns never change metrics
   - `_hp_run_length` string-method path vs the walk; `_expand_alignment` vs `edlib.getNiceAlignment`; CIGAR tail parsing vs full tokenizing
   - Output pinned from 0.2.10 for unusual characters and scoring-region edge cases (each case was chosen because it exposes a specific plausible mistake)

5. **test_utilities.py**: Utility function tests
   - Data class validation and immutability (including ScoringMode)
   - IUPAC nucleotide equivalence
   - CIGAR string parsing
   - Scoring region identification (LOCAL and GLOBAL modes)

## Dependencies

The package requires:
- `edlib>=1.3.9` - Fast sequence alignment library

Previous versions required BioPython, but v0.1.2+ includes a custom reverse complement
implementation with full IUPAC support, eliminating this heavyweight dependency.

Development dependencies:
- `pytest>=7.0.0` - Testing framework
- `pytest-cov>=4.0.0` - Coverage reporting

## Running and Testing Examples

### Basic Usage Test
```bash
python3 -c "
from adjusted_identity import align_and_score
result = align_and_score('AAATTTGGG', 'AAAATTTGGG')
print(f'Identity: {result.identity:.3f}')
print(f'Coverage: {result.seq1_coverage:.3f}')
"
```

### Package Import Test
```bash
python3 -c "import adjusted_identity; print('Package imports successfully')"
```

### Test Specific Functionality
```bash
# Test homopolymer adjustment
pytest tests/test_score_alignment.py::TestHomopolymerAdjustment -v

# Test IUPAC handling
pytest tests/test_score_alignment.py::TestIUPACAdjustment -v
```

## Development Guidelines

- All tests should pass before committing
- Tests serve as documentation - write clear test names and docstrings
- Add new tests for any new functionality
- Follow existing code style and patterns
- Update README.md for user-facing changes