# Adjusted Identity Calculator for DNA Sequences

A Python package implementing MycoBLAST-style sequence identity calculations for DNA sequences, specifically designed for mycological DNA barcoding applications. This package provides sophisticated sequence alignment and scoring with various adjustments for homopolymer differences, IUPAC ambiguity codes, and sequencing artifacts.

**Based on the MycoBLAST algorithm developed by Stephen Russell and Mycota Lab.** See the foundational article: ["Why NCBI BLAST identity scores can be misleading for fungi"](https://mycotalab.substack.com/p/why-ncbi-blast-identity-scores-can) which explains the theoretical basis and motivation for these sequence preprocessing adjustments.

## Features

- **Homopolymer Length Normalization**: Ignore differences in homopolymer run lengths (e.g., "AAA" vs "AAAA")
- **Repeat Motif Adjustment**: Handle dinucleotide and longer repeat motifs (e.g., "ATATAT" vs "ATATATAT")
- **IUPAC Ambiguity Code Handling**: Allow different ambiguity codes to match via nucleotide intersection
- **End Trimming**: Skip mismatches in terminal regions to avoid sequencing artifacts (automatic for sequences >40bp)
- **Indel Normalization**: Count contiguous indels as single evolutionary events
- **Comprehensive Alignment**: Multi-stage bidirectional alignment optimization using edlib
- **Flexible Configuration**: Enable/disable individual adjustments as needed

## Installation

### From GitHub

```bash
pip install git+https://github.com/joshuaowalker/adjusted-identity.git
```

### Development Installation

```bash
git clone https://github.com/joshuaowalker/adjusted-identity.git
cd adjusted-identity
pip install -e ".[dev]"
```

## Quick Start

### Option 1: Complete Solution (align + score)
```python
from adjusted_identity import align_and_score

# For raw sequences - handles alignment and scoring
result = align_and_score("ATCGAAAAATGTC", "ATCGAAAATGTC")
print(f"Adjusted identity: {result.identity:.3f}")
print(f"Coverage: seq1={result.seq1_coverage:.3f}, seq2={result.seq2_coverage:.3f}")
```

### Option 2: Core Scoring Only (use with your alignments)
```python
from adjusted_identity import score_alignment

# For pre-aligned sequences from BLAST, BioPython, etc.
aligned_seq1 = "ATCG-AAAT"  # From your alignment tool
aligned_seq2 = "ATCGAAAAT"  # From your alignment tool

result = score_alignment(aligned_seq1, aligned_seq2)
print(f"Adjusted identity: {result.identity:.3f}")
```

### Compare with Traditional Identity
```python
from adjusted_identity import align_and_score, RAW_ADJUSTMENT_PARAMS

# Traditional identity (no adjustments)
raw_result = align_and_score("ATCGAAAAATGTC", "ATCGAAAATGTC", RAW_ADJUSTMENT_PARAMS)
print(f"Traditional identity: {raw_result.identity:.3f}")

# Adjusted identity (with MycoBLAST adjustments)  
adj_result = align_and_score("ATCGAAAAATGTC", "ATCGAAAATGTC")
print(f"Adjusted identity: {adj_result.identity:.3f}")

# Examine the scoring pattern
print(f"Score pattern: {adj_result.score_aligned}")
```

## Use Cases

### Mycological DNA Barcoding

Common scenario: ITS sequences with homopolymer differences due to sequencing artifacts.

```python
from adjusted_identity import align_and_score

# Fungal ITS sequences with different homopolymer lengths
its_seq1 = "TCCGTAGGTGAACCTGCGGAAGGATCATTACCGAGTTTAAA"    # 3 A's at end
its_seq2 = "TCCGTAGGTGAACCTGCGGAAGGATCATTACCGAGTTTTAAAA"   # 4 A's at end

result = align_and_score(its_seq1, its_seq2)
print(f"Species identity: {result.identity:.3f}")  # Should be ~1.0 with homopolymer adjustment
```

### Handling Ambiguous Bases

```python
from adjusted_identity import align_and_score

# Sequences with IUPAC ambiguity codes
barcode1 = "ATCGRGTC"  # R = A or G
barcode2 = "ATCGKGTC"  # K = G or T (both R and K contain G)

result = align_and_score(barcode1, barcode2)
print(f"Identity with IUPAC handling: {result.identity:.3f}")  # Should be 1.0
```

### Repeat Motif Handling

```python
from adjusted_identity import align_and_score, AdjustmentParams

# Dinucleotide repeat differences (AT repeat from Russell article)
seq1 = "CGATAT--C"  # Missing one AT unit
seq2 = "CGATATATC"  # Has extra AT unit

# With repeat motif adjustment (default)
result = align_and_score(seq1, seq2)
print(f"Adjusted identity: {result.identity:.3f}")  # Should be 1.0

# Control max repeat motif length
params = AdjustmentParams(
    max_repeat_motif_length=3  # Detect up to trinucleotide repeats (e.g., CAG)
)
result = align_and_score("CAGCAG---TTC", "CAGCAGCAGTTC", params)
```

### Custom Adjustments

```python
from adjusted_identity import align_and_score, AdjustmentParams

# Custom adjustment parameters
custom_params = AdjustmentParams(
    normalize_homopolymers=True,   # Enable homopolymer adjustment
    handle_iupac_overlap=False,    # Disable IUPAC intersection
    normalize_indels=True,         # Enable indel normalization
    end_skip_distance=10,         # Skip 10bp from each end (instead of default 20)
    max_repeat_motif_length=2     # Detect up to dinucleotide repeats (default)
)

result = align_and_score(seq1, seq2, custom_params)
```

## API Architecture

This package provides a **layered API design** that separates sequence alignment from identity scoring, giving you maximum flexibility:

### Core Layer: `score_alignment()`
The **core implementation** that applies MycoBLAST-style adjustments to any pre-aligned sequences. This allows you to use **any alignment library** (BLAST, BioPython, edlib, etc.) with the adjusted identity algorithm.

**Input**: Just needs gapped sequences of equal length  
**Output**: Adjusted identity metrics with detailed scoring information

### Convenience Layer: `align_and_score()`
A **higher-level function** that combines fast edlib alignment with the scoring algorithm. Provides BLAST-like infix alignment that's fast enough for production use without requiring external dependencies.

**Input**: Raw unaligned sequences  
**Output**: Complete alignment and adjusted identity results

---

## API Reference

### Core Function

#### `score_alignment(seq1_aligned, seq2_aligned, adjustment_params=None, scoring_format=None)`

**The core implementation** - applies MycoBLAST-style adjustments to pre-aligned sequences from any source.

**Use this when:**
- You already have alignments from BLAST, BioPython, or other tools
- You want to integrate adjusted identity into existing pipelines
- You need maximum control over the alignment process

**Parameters:**
- `seq1_aligned`, `seq2_aligned` (str): Pre-aligned sequences with gaps (must be same length)
- `adjustment_params` (AdjustmentParams, optional): Adjustment parameters
- `scoring_format` (ScoringFormat, optional): Scoring visualization format

**Returns:**
- `AlignmentResult`: Scoring results and metrics

**Example with BLAST alignment:**
```python
# After getting BLAST alignment results
from adjusted_identity import score_alignment

blast_seq1 = "ATCG-AAAT"  # From BLAST output
blast_seq2 = "ATCGAAAAT"  # From BLAST output

result = score_alignment(blast_seq1, blast_seq2)
print(f"Adjusted identity: {result.identity:.3f}")
```

### Convenience Function

#### `align_and_score(seq1, seq2, adjustment_params=None, scoring_format=None)`

**High-level convenience function** that handles both alignment and scoring in one step.

**Use this when:**
- You want a simple, fast solution without additional alignment tools
- You need BLAST-like performance for production use
- You're comparing raw sequences end-to-end

**Parameters:**
- `seq1`, `seq2` (str): Raw DNA sequences to compare
- `adjustment_params` (AdjustmentParams, optional): Adjustment parameters
- `scoring_format` (ScoringFormat, optional): Scoring visualization format

**Returns:**
- `AlignmentResult`: Contains identity metrics, alignment, and coverage information

**Example:**
```python
from adjusted_identity import align_and_score

result = align_and_score("ATCGAAAATGTC", "ATCGAAAATGTC")
print(f"Identity: {result.identity:.3f}")
```

### Configuration Classes

#### `AdjustmentParams`

Configure which sequence adjustments to apply:

```python
AdjustmentParams(
    normalize_homopolymers=True,    # Ignore homopolymer length differences
    handle_iupac_overlap=True,      # Allow IUPAC ambiguity intersections
    normalize_indels=True,          # Count contiguous indels as single events
    end_skip_distance=20,          # Skip first/last N nucleotides (not positions) from scoring
    max_repeat_motif_length=2      # Maximum repeat motif length to detect (1=homopolymers only, 2=dinucleotides, etc.)
)
```

#### `ScoringFormat`

Customize alignment visualization characters:

```python
ScoringFormat(
    match='|',                     # Exact match or IUPAC equivalent
    substitution=' ',              # Nucleotide substitution
    indel_start=' ',               # First position of indel
    indel_extension='-',           # Indel positions (normalization)
    homopolymer_extension='=',     # Homopolymer length difference
    end_trimmed='.'               # Position outside scoring region
)
```

#### `AlignmentResult`

Results returned by alignment functions:

```python
AlignmentResult(
    identity=0.95,                 # Identity score (0.0-1.0)
    mismatches=2,                  # Number of mismatches counted
    scored_positions=40,           # Positions used for identity calculation
    seq1_coverage=0.98,            # Fraction of seq1 in alignment
    seq2_coverage=0.97,            # Fraction of seq2 in alignment
    seq1_aligned="ATCG-ATCG",      # Aligned sequence 1 with gaps
    seq2_aligned="ATCGATCG-",      # Aligned sequence 2 with gaps
    score_aligned="||||=|||"      # Scoring visualization
)
```

### Constants

- `DEFAULT_ADJUSTMENT_PARAMS`: All adjustments enabled (recommended)
- `RAW_ADJUSTMENT_PARAMS`: No adjustments (traditional identity)
- `DEFAULT_SCORING_FORMAT`: Default visualization characters

## Integration with Other Alignment Tools

The core `score_alignment()` function works with alignments from any source. Here are examples with popular alignment libraries:

### Using with NCBI BLAST

```python
# NOTE: This example is illustrative and not tested
from Bio.Blast import NCBIWWW, NCBIXML
from adjusted_identity import score_alignment

# Run BLAST search (example)
result_handle = NCBIWWW.qblast("blastn", "nt", query_sequence)
blast_records = NCBIXML.parse(result_handle)

for blast_record in blast_records:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            # Extract aligned sequences from BLAST HSP
            query_aligned = hsp.query
            subject_aligned = hsp.sbjct
            
            # Apply adjusted identity scoring
            adj_result = score_alignment(query_aligned, subject_aligned)
            
            print(f"BLAST identity: {hsp.identities/hsp.align_length:.3f}")
            print(f"Adjusted identity: {adj_result.identity:.3f}")
```

### Using with BioPython PairwiseAligner

```python
from Bio import Align
from adjusted_identity import score_alignment

# Create BioPython aligner
aligner = Align.PairwiseAligner()
aligner.match_score = 2
aligner.mismatch_score = -1

# Example sequences
seq1 = "ATCGATCG"
seq2 = "ATCGTCG"  # Missing one nucleotide

# Perform alignment
alignments = aligner.align(seq1, seq2)
best_alignment = alignments[0]

# Extract aligned sequences using indexing (idiomatic BioPython)
seq1_aligned = str(best_alignment[0])
seq2_aligned = str(best_alignment[1])

# Apply adjusted scoring
result = score_alignment(seq1_aligned, seq2_aligned)
print(f"Adjusted identity: {result.identity:.3f}")
```

### Using with Custom/External Aligners

```python
# NOTE: This example is illustrative template code
from adjusted_identity import score_alignment

def process_alignment_file(alignment_file):
    """Process alignments from any external tool."""
    results = []
    
    # Parse your alignment format (FASTA, SAM, custom, etc.)
    for seq1_aligned, seq2_aligned in parse_alignment_file(alignment_file):
        # Ensure sequences are the same length
        assert len(seq1_aligned) == len(seq2_aligned)
        
        # Apply adjusted identity scoring
        result = score_alignment(seq1_aligned, seq2_aligned)
        results.append(result)
    
    return results
```

## Understanding End Trimming Behavior

The `end_skip_distance` parameter implements "digital end trimming" to skip sequencing artifacts near read ends. **Important**: This parameter counts **nucleotides** (non-gap characters), not alignment positions.

### Automatic Activation

End trimming only activates when sequences are long enough:

```python
from adjusted_identity import align_and_score, AdjustmentParams

# Short sequences (< 2 × end_skip_distance nucleotides): NO trimming applied
short_seq1 = "ATCGATCG"      # 8 nucleotides 
short_seq2 = "ATCGATCG"      # 8 nucleotides
result = align_and_score(short_seq1, short_seq2)  # Uses default end_skip_distance=20
print(f"Scored positions: {result.scored_positions}")  # 8 (full sequence)
print(f"Score pattern: {result.score_aligned}")        # "||||||||" (no trimming dots)

# Long sequences (≥ 2 × end_skip_distance nucleotides): Trimming applied  
long_seq1 = "A" * 25 + "TCGX" + "T" * 25    # 54 nucleotides
long_seq2 = "A" * 25 + "TCGA" + "T" * 25    # 54 nucleotides
result = align_and_score(long_seq1, long_seq2)
print(f"Scored positions: {result.scored_positions}")  # ~14 (middle region only)
print(f"Score pattern: {result.score_aligned}")        # ".......|||| |||||......." (dots show trimmed regions)
```

### Nucleotide vs Position Counting

End trimming counts **actual nucleotides** in each sequence, ignoring gaps:

```python
# This alignment has gaps, but nucleotide counting still works correctly
seq1_aligned = "AAA---TCGATCG---TTT"  # 12 nucleotides (ignoring gaps)
seq2_aligned = "---AAATCGATCGTTT---"  # 12 nucleotides (ignoring gaps)

# With end_skip_distance=5: skips first 5 and last 5 nucleotides from each sequence
# Only the middle "TCGATCG" region (2 nucleotides) would be scored
```

### Customizing End Trimming

```python
# Disable end trimming completely
no_trim_params = AdjustmentParams(end_skip_distance=0)
result = align_and_score(long_seq1, long_seq2, no_trim_params)

# Use shorter trimming distance for smaller sequences
short_trim_params = AdjustmentParams(end_skip_distance=5) 
result = align_and_score(medium_seq1, medium_seq2, short_trim_params)
```

**Rule of thumb**: For sequences shorter than `2 × end_skip_distance` nucleotides, end trimming has no effect and the entire alignment is scored.

## Advanced Usage

### Batch Processing

```python
from adjusted_identity import align_and_score

sequences = [
    ("seq1", "ATCGATCGATCG"),
    ("seq2", "ATCGATCGATCC"),
    ("seq3", "ATCGATCGATCG"),
]

reference = sequences[0][1]
results = []

for name, seq in sequences[1:]:
    result = align_and_score(reference, seq)
    results.append({
        'name': name,
        'identity': result.identity,
        'coverage': min(result.seq1_coverage, result.seq2_coverage)
    })

# Sort by identity
results.sort(key=lambda x: x['identity'], reverse=True)
```

### Understanding Scoring Patterns

The `score_aligned` field shows how each position was scored:

```python
result = align_and_score("AAA-TTT", "AAAATTT")
print(result.score_aligned)  # "|||=|||"
# | = match
# = = homopolymer extension (ignored if adjustment enabled)
```

Scoring symbols:
- `|`: Exact match or IUPAC equivalent
- ` `: Nucleotide substitution (counts as mismatch)
- ` `: First position of indel (counts as mismatch)
- `-`: Indel extension positions (ignored if normalization enabled)
- `=`: Homopolymer extension (ignored if adjustment enabled)
- `.`: End-trimmed position (ignored in scoring)

## Testing

Run the comprehensive test suite:

```bash
# Install development dependencies
pip install -e ".[dev]"

# Run tests
pytest

# Run tests with coverage
pytest --cov=adjusted_identity --cov-report=html
```

The test suite includes:
- Unit tests for all adjustment types
- Edge cases and error conditions
- Real-world mycological scenarios
- Performance tests with long sequences
- Documentation examples

## Background

This package implements the sequence preprocessing approach described in the MycoBLAST algorithm by Stephen Russell and Mycota Lab, adapted for general-purpose DNA sequence comparison. The foundational research is detailed in ["Why NCBI BLAST identity scores can be misleading for fungi"](https://mycotalab.substack.com/p/why-ncbi-blast-identity-scores-can).

The adjustments are particularly valuable for:

- **Fungal taxonomy**: ITS sequences often have homopolymer differences
- **DNA barcoding**: Technical artifacts can obscure phylogenetic signal
- **Sequence quality assessment**: End-trimming handles poor-quality regions
- **Phylogenetic analysis**: IUPAC codes preserve ambiguous but valid matches

**Credit**: This implementation is based on the MycoBLAST algorithm developed by Stephen Russell and Mycota Lab. The theoretical framework and biological motivation are thoroughly explained in their foundational article.

## Contributing

1. Fork the repository
2. Create a feature branch: `git checkout -b feature-name`
3. Make changes and add tests
4. Run the test suite: `pytest`
5. Submit a pull request

## Citation

If you use this package in your research, please cite:

```
Walker, J. (2025). Adjusted Identity Calculator for DNA Sequences. 
GitHub: https://github.com/joshuaowalker/adjusted-identity
```

**Please also cite the foundational work:**

```
Russell, S. (2024). Why NCBI BLAST identity scores can be misleading for fungi. 
Mycota Lab. https://mycotalab.substack.com/p/why-ncbi-blast-identity-scores-can
```

## License

BSD 2-Clause License - see [LICENSE](LICENSE) file for details.

## Changelog

### Version 0.1.1
- Added repeat motif adjustment support (dinucleotide and longer repeats)
- Implemented intelligent motif length detection with degeneracy handling
- Added `max_repeat_motif_length` parameter to AdjustmentParams
- Enhanced left-right indel processing algorithm for mixed motif lengths
- Added comprehensive test coverage for repeat motif scenarios

### Version 0.1.0
- Initial release
- Complete MycoBLAST-style adjustment implementation (except repeat motifs)
- Comprehensive test suite
- Full documentation and examples