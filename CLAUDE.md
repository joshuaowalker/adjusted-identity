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
   - `ScoringFormat` - Format codes for alignment visualization
   - `AdjustmentParams` - Parameters controlling sequence adjustments (includes `scoring_mode`)

2. **Core Functions** (`adjusted_identity/__init__.py`):
   - `align_and_score()` - Main entry point for sequence comparison (routes to HW or NW alignment based on scoring mode)
   - `align_edlib_bidirectional()` - Multi-stage semi-global alignment optimization (used by LOCAL mode)
   - `_align_nw()` - Needleman-Wunsch global alignment wrapper (used by GLOBAL mode)
   - `score_alignment()` - Scoring with configurable adjustments

3. **Adjustment Features**:
   - Homopolymer length normalization - ignores differences in homopolymer run lengths
   - IUPAC ambiguity code handling - allows different ambiguity codes to match
   - End trimming - skips mismatches in terminal regions (disabled by default, set `end_skip_distance` to enable)
   - Indel normalization - counts contiguous indels as single events
   - LOCAL/GLOBAL scoring modes - controls whether terminal gaps affect identity

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

2. **test_alignment.py**: Alignment algorithm tests
   - `TestAlignEdlibBidirectional` - Bidirectional alignment
   - `TestAlignAndScore` - End-to-end functionality
   - `TestRealWorldScenarios` - Mycological sequence examples
   - `TestNWAlignment` - Needleman-Wunsch alignment wrapper
   - `TestGlobalModeEndToEnd` - GLOBAL scoring mode end-to-end tests

3. **test_utilities.py**: Utility function tests
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