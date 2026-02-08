"""
Genomic Clustering Demo
Demonstrates the complete pipeline from FASTA parsing to phylogenetic tree generation.
"""

import sys
from fasta_parser import ParseSeqFile
from sequence_alignment import AlignByDP
from distance_matrix import ComputeDistMatrix
from clustering import Cluster


def print_separator():
    """Print a visual separator."""
    print("\n" + "=" * 70 + "\n")


def main():
    """Run the complete genomic clustering pipeline."""
    
    # Check command-line arguments
    if len(sys.argv) < 2:
        print("Usage: python demo.py <fasta_file>")
        print("\nExample:")
        print("  python demo.py example_sequences.fasta")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    
    try:
        # Step 1: Parse FASTA file
        print_separator()
        print("STEP 1: Parsing FASTA file")
        print(f"File: {fasta_file}")
        print_separator()
        
        sequences = ParseSeqFile(fasta_file)
        print(f"Successfully parsed {len(sequences)} sequences:\n")
        
        for idx, (label, seq) in enumerate(sequences):
            print(f"  {idx}. {label}: {seq[:50]}{'...' if len(seq) > 50 else ''} (length: {len(seq)})")
        
        # Step 2: Perform pairwise alignments
        print_separator()
        print("STEP 2: Performing pairwise sequence alignments")
        print("Algorithm: Needleman-Wunsch (match=+5, mismatch=-2, gap=-6)")
        print_separator()
        
        alignments = AlignByDP(sequences)
        print(f"Completed {len(alignments)} pairwise alignments\n")
        
        # Show first alignment as example
        if alignments:
            first_key = list(alignments.keys())[0]
            i, j = first_key
            aligned1, aligned2 = alignments[first_key]
            print(f"Example alignment ({sequences[i][0]} vs {sequences[j][0]}):")
            print(f"  Seq1: {aligned1[:60]}{'...' if len(aligned1) > 60 else ''}")
            print(f"  Seq2: {aligned2[:60]}{'...' if len(aligned2) > 60 else ''}")
        
        # Step 3: Compute distance matrix
        print_separator()
        print("STEP 3: Computing evolutionary distance matrix")
        print("Method: Jukes-Cantor correction")
        print_separator()
        
        distance_matrix = ComputeDistMatrix(alignments)
        print("Distance matrix computed successfully\n")
        
        # Display distance matrix
        labels = [seq[0] for seq in sequences]
        print("Evolutionary distances:")
        print(f"\n{'':15s}", end="")
        for label in labels:
            print(f"{label[:12]:>12s}", end=" ")
        print()
        
        for i, label in enumerate(labels):
            print(f"{label[:15]:15s}", end="")
            for j in range(len(labels)):
                print(f"{distance_matrix[i][j]:12.4f}", end=" ")
            print()
        
        # Step 4: Perform hierarchical clustering
        print_separator()
        print("STEP 4: Building phylogenetic tree")
        print("Method: WPGMA (Weighted Pair Group Method with Arithmetic Mean)")
        print_separator()
        
        tree = Cluster(distance_matrix, labels)
        
        print("Phylogenetic tree (Newick-like format):\n")
        print(f"  {tree}\n")
        
        print_separator()
        print("✓ Pipeline completed successfully!")
        print_separator()
        
    except FileNotFoundError:
        print(f"Error: File '{fasta_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
