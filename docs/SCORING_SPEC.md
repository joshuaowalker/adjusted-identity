# Adjusted Identity Scoring Specification

This document specifies the behavior of the adjusted-identity scoring algorithm (v0.2.x variant range algorithm).

---

## 1. Core Concepts

### 1.1 Variant Range

A **variant range** is a maximal contiguous region of an alignment where at least one position is NOT a match.

**Match definition**: A position is a match if:
- Both sequences have non-gap characters that are equivalent (exact match or IUPAC equivalence)
- OR both sequences have gaps (dual-gap)

Variant ranges are bounded by match positions on the left and right (when available).

### 1.2 Allele Extraction

For each variant range, we extract the **allele** from each sequence:
- The allele is the gap-free content (all non-gap characters) within the variant range
- We also track the source positions for visualization

Example:
```
Alignment:  AAA-TT-GGG
               ^^^^^ variant range (positions 3-7)
Allele:     "TT" from positions [4, 5]
```

### 1.3 Extension Detection

An **extension** is a portion of an allele that matches the context (adjacent nucleotides outside the variant range).

- **Left extension**: Characters at the start of the allele that match a repeat of the left context
- **Right extension**: Characters at the end of the allele that match a repeat of the right context

**Critical rule**: Extensions are only detected from **positions of agreement** (match boundaries). The context used for extension detection comes from matched positions outside the variant range, not from within variant ranges or other non-matching regions.

Extensions are detected by trying motif lengths from largest (`max_repeat_motif_length`, default 2) down to 1.

### 1.4 Core Content

The **core** is the portion of an allele that remains after removing left and right extensions. It represents the "unexplained" content that cannot be attributed to homopolymer/repeat expansion.

### 1.5 Occam's Razor Scoring

We apply Occam's razor: prefer the simplest explanation.

- If both alleles are **pure extensions** (no core): 0 edits (they're equivalent expansions of the same repeat)
- If one allele is pure extension and the other has core: the core represents edits
- If both have cores: compare the cores to determine substitutions and length differences

---

## 2. Variant Range Detection

### 2.1 What Constitutes a Match

A position `i` is a **match** for variant range boundary detection if:

| seq1[i] | seq2[i] | Is Match? | Reason |
|---------|---------|-----------|--------|
| `A` | `A` | Yes | Exact match |
| `A` | `a` | Yes | Case-insensitive match |
| `R` | `A` | Yes (if IUPAC enabled) | R={A,G} contains A |
| `R` | `Y` | Depends | Only if intersection exists and `handle_iupac_overlap=True` |
| `-` | `-` | Yes | Dual-gap treated as match |
| `A` | `-` | No | Gap in one sequence |
| `-` | `A` | No | Gap in one sequence |
| `A` | `T` | No | Mismatch |

### 2.2 Dual-Gap Handling

Dual-gap positions (both sequences have `-`) are treated as matches for boundary detection. This prevents MSA artifacts from artificially splitting variant ranges.

### 2.3 Variant Range Boundaries

Each variant range records:
- `start`: First position of the variant range
- `end`: Last position of the variant range (inclusive)
- `left_bound`: Position of the match immediately before start (-1 if at alignment start)
- `right_bound`: Position of the match immediately after end (-1 if at alignment end)

---

## 3. Context Extraction

### 3.1 MSA Consensus Rules

Context is extracted from positions outside the variant range using consensus rules:

| seq1[pos] | seq2[pos] | Result |
|-----------|-----------|--------|
| Same char | Same char | Use it (consensus) |
| `A` | `a` | Use it (case-insensitive match) |
| Char | `-` | Use the char |
| `-` | Char | Use the char |
| `-` | `-` | Skip position (doesn't count toward context length) |
| `A` | `T` | **Conflict** - return None (no valid context) |

### 3.2 Context Direction

- **Left context**: Extracted backwards from variant start, returned in left-to-right order
- **Right context**: Extracted forwards from variant end, returned in left-to-right order

### 3.3 Insufficient Context

If we cannot collect the requested number of context characters (due to reaching alignment boundary or encountering conflicts), context extraction returns `None` and no extension detection is performed for that side.

---

## 4. Extension Analysis

### 4.1 Motif Length Priority

Extension matching tries motif lengths from largest to smallest:
1. Try `max_repeat_motif_length` (default: 2)
2. If no match, try `max_repeat_motif_length - 1`
3. Continue down to 1
4. Use first motif length that matches

### 4.2 Degenerate Motif Handling

If a longer motif consists of all identical characters (e.g., "AA"), it is collapsed to a homopolymer (length 1). This prevents double-counting of homopolymer extensions.

Example: Context "AA" with motif_len=2 → treated as homopolymer "A" with motif_len=1

### 4.3 Complete Motif Requirement

Only **complete** motif matches are consumed as extensions. Partial motifs remain as core content.

Example:
```
Left context: "AT" (dinucleotide)
Allele: "ATATG"
Result: "ATAT" consumed as 2 complete motifs, "G" is core
```

### 4.4 IUPAC in Extensions

When `handle_iupac_overlap=True`, extension matching uses IUPAC equivalence:
- `R` in allele can match `A` in context (since R={A,G})
- Different ambiguity codes can match if they have intersection

---

## 5. Scoring Rules

### 5.1 Pure Extension Cases

| Allele 1 | Allele 2 | Edits | Scored Positions |
|----------|----------|-------|------------------|
| Pure extension | Pure extension | 0 | 0 |
| Pure extension | Has core | len(core2) or 1 if normalized | len(core2) or 1 |
| Has core | Pure extension | len(core1) or 1 if normalized | len(core1) or 1 |
| Empty | Empty | 0 | 0 |

### 5.2 Core Comparison

When both alleles have core content:

1. If cores are identical: 0 edits, scored_positions = len(core)
2. If cores differ:
   - Count substitutions (positions where chars don't match, using IUPAC if enabled)
   - Calculate length difference
   - If `normalize_indels`: length difference counts as 1 edit
   - If not: length difference counts as abs(len1 - len2) edits

### 5.3 Indel Normalization

When `normalize_indels=True`:
- Contiguous indel regions count as 1 edit regardless of length
- This applies to core content that represents insertions/deletions

---

## 6. Parameter Interactions

### 6.1 normalize_homopolymers

| Value | Effect |
|-------|--------|
| `True` (default) | Extensions detected and ignored in scoring |
| `False` | All allele content treated as indels (no extension detection) |

### 6.2 normalize_indels

| Value | Effect |
|-------|--------|
| `True` (default) | Contiguous indels = 1 edit |
| `False` | Each position in indel = 1 edit |

### 6.3 handle_iupac_overlap

| Value | Effect |
|-------|--------|
| `True` (default) | Different IUPAC codes with intersection are equivalent |
| `False` | IUPAC codes must match exactly (or one must be standard nucleotide) |

### 6.4 Combined Effects

| HP Norm | Indel Norm | Behavior |
|---------|------------|----------|
| True | True | Extensions ignored, remaining core normalized |
| True | False | Extensions ignored, core positions counted individually |
| False | True | All content as indels, normalized to 1 per region |
| False | False | All content as indels, counted per position |

---

## 7. Visualization

### 7.1 Score String Markers

| Marker | Meaning | Scored? |
|--------|---------|---------|
| `\|` | Match (exact or IUPAC equivalent) | Yes |
| `=` | Homopolymer/repeat extension OR ambiguous match (IUPAC intersection) | See note |
| `-` | Indel extension (when normalize_homopolymers=False) | No |
| ` ` (space) | Mismatch or indel start | Yes |
| `.` | End-trimmed position OR dual-gap (MSA artifact) | No |

**Note on `=` marker**: The `=` character serves dual purposes:
- Homopolymer/repeat extensions: Not scored (extensions don't count as mismatches)
- Ambiguous IUPAC matches: Scored as a match position

### 7.1.1 Dual-Gap Handling

**Dual-gap positions** (where both sequences have `-`) represent MSA alignment artifacts and are:
- Marked with `.` in the visualization (same as end-trimmed positions)
- NOT counted in `scored_positions`
- NOT counted as mismatches

This is a change from v0.1.x which counted dual-gaps in `scored_positions`.

### 7.1.2 Overhang Handling

**Overhang positions** occur when one sequence has a gap at the alignment boundary while the other has nucleotide content. Examples:

```
seq1: -ATT    seq1: ATT-
seq2: XATT    seq2: ATTX
score: .|||   score: |||.
```

Overhangs are:
- Marked with `.` (same as end-trimmed and dual-gap positions)
- NOT counted in `scored_positions`
- NOT counted as mismatches

**Rationale**: Overhangs represent regions where only one sequence has content, typically occurring at alignment boundaries. Since there's no aligned content from the other sequence to compare against, these positions cannot meaningfully contribute to identity scoring.

### 7.2 Visualization String

The algorithm produces a single score string (`score_aligned`) that visualizes the scoring for each position. To get the scoring from the opposite sequence's perspective, call `score_alignment()` with the arguments swapped.

### 7.3 Asymmetric Visualization

When one position is an extension in one sequence and core in the other:
- Extension position shows extension marker (`=`)
- Core position shows match marker (`|`) if cores match

---

## 8. Examples

The following examples were harvested from empirical analysis of real-world MSA data
(ONT fungal barcoding sequences). Each represents a distinct pattern.

### 8.1 Dual-Gap Not Scored

```
seq1: AA--TT
seq2: AA--TT
score: ||..||
```

- Dual-gaps at positions 2-3 marked with `.`
- `scored_positions = 4` (dual-gaps excluded)
- `identity = 1.0`, `mismatches = 0`

**Rationale**: Dual-gaps are MSA artifacts with no biological meaning.

### 8.2 Opposite Direction Extensions

```
seq1: TGC-C-TC
seq2: TGCT--TC
score: |||==.||
```

- Variant range at positions 3-5
- allele1 = "C", allele2 = "T"
- C extends left context (C), T extends right context (T)
- Both are pure extensions → 0 edits

**Rationale**: Occam's razor - both placements are valid expansions.

### 8.3 Same Direction Extensions

```
seq1: AAG-G-CC
seq2: AAGG--CC
score: |||==.||
```

- Both G's extend the left G context
- `mismatches = 0`

### 8.4 Substitution Not Extension

```
seq1: GG-T-AA
seq2: GG-C-AA
score: ||. .||
```

- Variant range at position 3
- T doesn't extend G or A → core = "T"
- C doesn't extend G or A → core = "C"
- T ≠ C → 1 substitution

**Rationale**: Extensions only from positions of agreement.

### 8.5 One Extension, One Core

```
seq1: AAA-A-TTT
seq2: AAA--XTTT
score: |||.= |||
```

- Position 3: dual-gap (`.`)
- Position 4: allele1 = "A" extends left A context → extension (`=`)
- Position 5: allele2 = "X" doesn't extend A or T → core, mismatch (` `)
- Core content X counts as 1 mismatch

### 8.6 Context Skips Dual-Gaps

```
seq1: TT--T-GG
seq2: TT---TGG
score: ||..==||
```

- Positions 2-3: dual-gaps (`..`)
- Context extraction skips dual-gaps to find T context at position 1
- Positions 4-5: Both T's in the variant range extend T context → both marked as extensions (`==`)
- Result: 0 mismatches

### 8.7 Gap-Free Substitution as Extensions

```
seq1: ATTCA
seq2: ATCCA
score: ||=||
```

- Single mismatch at position 2: T vs C
- T extends left context (T homopolymer)
- C extends right context (C homopolymer)
- Both are pure extensions → 0 edits

**Rationale**: Standard alignment algorithms (e.g., edlib) don't implement homopolymer-aware scoring. If they did, they might produce `ATT-CA` vs `AT-CCA` instead, showing explicit homopolymer extensions. Since we can't modify the aligner, adjusted-identity compensates post-alignment by recognizing that both nucleotides are valid homopolymer extensions of their respective contexts.

This is a key feature of the variant range algorithm: it handles cases where the aligner produced a gap-free alignment but the "substitution" is actually equivalent homopolymer expansions that could have been aligned with gaps.

### 8.8 IUPAC Ambiguous Match

```
seq1: AAARTTT
seq2: AAAGTTT
score: |||=|||
```

- Position 3: R vs G
- R = {A, G}, G is in that set → ambiguous match
- Marked with `=`, scored as match (0 mismatches)
- `scored_positions = 7`, `identity = 1.0`

**Rationale**: IUPAC ambiguity codes represent uncertainty in the base call. When the other sequence has a nucleotide within the ambiguous set, it's a valid match.

---

## 9. Differences from v0.1.x

### 9.1 Dual-Gap Handling

| Scenario | v0.1.x Behavior | v0.2.x Behavior | Rationale |
|----------|-----------------|-----------------|-----------|
| Dual-gap positions (`-` in both sequences) | Counted in `scored_positions`, marked as `\|` | NOT counted in `scored_positions`, marked as `.` | Dual-gaps are MSA artifacts with no biological meaning; excluding them produces more accurate identity scores |

### 9.2 Variant Range Algorithm

| Scenario | v0.1.x Behavior | v0.2.x Behavior | Rationale |
|----------|-----------------|-----------------|-----------|
| Alternating gaps (e.g., `TGC-C-TC` vs `TGCT--TC`) | Processed as separate indel events, may count mismatches | Identified as variant range, both alleles are pure extensions → 0 mismatches | Occam's razor: both placements are valid expansions of the context |
| Gap-free substitution with context match (e.g., `ATTCA` vs `ATCCA`) | Counted as 1 substitution | Recognized as opposite-direction homopolymer extensions → 0 mismatches | Compensates for aligners that don't implement homopolymer-aware scoring; the "substitution" represents equivalent expansions |
| Insertion not extending context (e.g., `-/G` between C and A) | Sometimes incorrectly treated as extension | Correctly identified as mismatch (G doesn't extend C or A) | Insertions must actually match adjacent context to be treated as extensions |
| Complex indel patterns | Position-by-position analysis | Holistic variant range analysis with allele extraction | More accurate for MSA data where gap placement is arbitrary |

---

## Appendix A: Glossary

- **Allele**: Gap-free content extracted from a sequence within a variant range
- **Core**: Portion of allele not explained by extensions
- **Extension**: Portion of allele matching adjacent context (repeat expansion)
- **Variant range**: Contiguous non-match region in alignment
- **Pure extension**: Allele with no core (entirely explained by extensions)
