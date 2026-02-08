import math
import numpy as np

def remove_gaps(seq1, seq2):
    """Remove the gap positions from a pair of aligned sequences.
    It returns a list of (char1, char2) pairs that represent the nucleotide pairs."""
    pairs = []
    for a,b in zip(seq1, seq2):
        if a != '-' and b!='-':
            pairs.append((a, b))
    return pairs

def count_differences(pairs):
    """Count how many positions differ.
    Returns:
    diff: number of positions where a != b
    total: total number of comparable positions."""
    diff = 0
    for a,b in pairs:
        if a!=b:
            diff+=1
    total = len(pairs)
    return diff, total


def compute_p_distance(diff, total):
    """Calculate p-distance (diff / total) between two sequences.
    For total = 0, the distance is defined as 0."""
    if total == 0:
        return 0.0
    return diff / total

def corrected_distance(p):
    """Apply Jukes-Cantor correction formula (-0.75 * ln(1 - (4/3) * p)).
    Returns 30 for p >= 0.75 as the formula becomes undefined."""
    if p >= 0.75:
        return 30.0
    return -0.75 * math.log(1 - (4/3) * p)

def compute_pair_distance(seq1, seq2):
    """Calculate the evolutionary distance value for one aligned pair."""
    pairs = remove_gaps(seq1, seq2)
    diff, total = count_differences(pairs)
    p = compute_p_distance(diff, total)
    return corrected_distance(p)

def ComputeDistMatrix(aligned_dict):
    """Build the evolutionary distance matrix from aligned sequence pairs.
    Input: dict(tuple(int, int) -> tuple(string, string))
    Output: list(list(float))."""
    if not isinstance(aligned_dict, dict) or not all(
            isinstance(k, tuple) and len(k) == 2 and isinstance(k[0], int) and isinstance(k[1], int)
            and isinstance(v, tuple) and len(v) == 2 and all(isinstance(s, str) for s in v)
            for k, v in aligned_dict.items()
    ):
        raise Exception("malformed input")

    indices = []
    for (i, j) in aligned_dict.keys():
        indices.append(i)
        indices.append(j)
    n = max(indices) + 1

    matrix = np.zeros((n, n))

    for (i, j), (seq1, seq2) in aligned_dict.items():
        d = compute_pair_distance(seq1, seq2)
        matrix[i, j] = d
        matrix[j, i] = d  

    return matrix.tolist()