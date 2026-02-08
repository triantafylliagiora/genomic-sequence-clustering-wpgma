# Genomic Clustering Tool

A Python-based bioinformatics pipeline for phylogenetic analysis of genomic sequences using sequence alignment and hierarchical clustering.

## Overview

This tool performs genomic sequence analysis through four main stages:
1. **FASTA Parsing** - Validates and parses genomic sequence files
2. **Sequence Alignment** - Performs pairwise alignment using the Needleman-Wunsch algorithm
3. **Distance Calculation** - Computes evolutionary distances using Jukes-Cantor correction
4. **Hierarchical Clustering** - Builds phylogenetic trees using WPGMA (Weighted Pair Group Method with Arithmetic Mean)

## Features

-  Robust FASTA file parsing with validation
-  Global sequence alignment (Needleman-Wunsch)
-  Evolutionary distance matrix computation
-  WPGMA clustering for phylogenetic tree construction
-  Handles sequences of varying lengths
-  Input validation and error handling

## Installation

### Prerequisites
- Python 3.7 or higher
- pip package manager

### Setup

1. Clone the repository:
```bash
git clone https://github.com/triantafylliagiora/genomic-sequence-clustering-wpgma.git
cd genomic_clustering
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### Quick Start

Run the complete pipeline with the demo script:

```bash
python demo.py example_sequences.fasta
```

### Individual Module Usage

#### 1. Parse FASTA File
```python
from fasta_parser import ParseSeqFile

sequences = ParseSeqFile("sequences.fasta")
# Returns: [(label1, sequence1), (label2, sequence2), ...]
```

#### 2. Align Sequences
```python
from sequence_alignment import AlignByDP

alignments = AlignByDP(sequences)
# Returns: {(0, 1): (aligned_seq1, aligned_seq2), ...}
```

#### 3. Compute Distance Matrix
```python
from distance_matrix import ComputeDistMatrix

distances = ComputeDistMatrix(alignments)
# Returns: 2D matrix of evolutionary distances
```

#### 4. Perform Clustering
```python
from clustering import Cluster

labels = [seq[0] for seq in sequences]
tree = Cluster(distances, labels)
# Returns: Newick-like tree string
```

## Input Format

FASTA files should follow this format:
```
>Label1 ATCGATCG
>Label2 GCTAGCTA
>Label3 TTAACCGG
```

**Requirements:**
- Each line starts with `>`
- Format: `>Label<whitespace>Sequence`
- Valid nucleotides: A, C, G, T
- Unique labels
- No duplicate sequences

## Algorithms

### Needleman-Wunsch Algorithm
- **Match score:** +5
- **Mismatch penalty:** -2
- **Gap penalty:** -6

### Jukes-Cantor Distance
Corrects for multiple substitutions at the same site:
```
d = -0.75 * ln(1 - (4/3) * p)
```
where `p` is the proportion of differing sites.

### WPGMA Clustering
Builds hierarchical trees by iteratively merging the closest clusters using weighted averaging.

## Project Structure

```
genomic_clustering/
├── fasta_parser.py         # FASTA file parsing and validation
├── sequence_alignment.py   # Needleman-Wunsch alignment
├── distance_matrix.py      # Evolutionary distance computation
├── clustering.py           # WPGMA clustering algorithm
├── demo.py                 # Full pipeline demonstration
├── example_sequences.fasta # Sample input file
├── requirements.txt        # Python dependencies
├── README.md              # This file
└── LICENSE                # License information
```

## Error Handling

The tool validates:
- File existence and format
- Sequence composition (only ACGT)
- Label uniqueness
- Non-empty sequences and labels
- Proper FASTA formatting


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Authors

Triantafyllia Giora
Jesus Esteban

## Acknowledgments

- Needleman-Wunsch algorithm for sequence alignment
- Jukes-Cantor model for evolutionary distance
- WPGMA method for hierarchical clustering
