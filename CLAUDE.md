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

1. **Data Classes** (`adjusted_identity/__init__.py:41-108`):
   - `AlignmentResult` - Contains alignment results and identity metrics
   - `ScoringFormat` - Format codes for alignment visualization
   - `AdjustmentParams` - Parameters controlling sequence adjustments

2. **Core Functions** (`adjusted_identity/__init__.py:653-739`):
   - `align_and_score()` - Main entry point for sequence comparison
   - `align_edlib_bidirectional()` - Multi-stage alignment optimization
   - `score_alignment()` - Scoring with configurable adjustments

3. **Adjustment Features**:
   - Homopolymer length normalization - ignores differences in homopolymer run lengths
   - IUPAC ambiguity code handling - allows different ambiguity codes to match
   - End trimming - skips mismatches in terminal regions (disabled by default, set `end_skip_distance` to enable)
   - Indel normalization - counts contiguous indels as single events

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

2. **test_alignment.py**: Alignment algorithm tests
   - `TestAlignEdlibBidirectional` - Bidirectional alignment
   - `TestAlignAndScore` - End-to-end functionality
   - `TestRealWorldScenarios` - Mycological sequence examples

3. **test_utilities.py**: Utility function tests
   - Data class validation and immutability
   - IUPAC nucleotide equivalence
   - CIGAR string parsing
   - Homopolymer detection

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