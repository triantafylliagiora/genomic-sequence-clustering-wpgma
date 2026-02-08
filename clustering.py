import numpy as np
from typing import Tuple, Union


def merge_labels(labels: list[str], row: int, col: int) -> list[str]:
    """Create a new cluster label by merging two existing labels.
    The function takes the labels at positions row and col, wraps them
    into one cluster, removes the original labels and inserts the new label at the start of the list."""
    new_label = f"({labels[row]},{labels[col]})"
    new_labels = [labels[k] for k in range(len(labels)) if k not in (row, col)]
    new_labels.insert(0, new_label)
    return new_labels

def reduction_matrix(matrix: np.ndarray, labels: list[str]) -> Tuple[np.ndarray, list[str]]:
    """This function identifies the pair of clusters with the smallest
    non‑zero distance, merges them into a new cluster and constructs
    a new distance matrix according to the WPGMA rule.
    Special cases:
    If all off‑diagonal entries are zero, the sequences are identical.
    If no valid pair is found, the input matrix is malformed."""
    lowest = float("inf")
    row, col = None, None
    nonzero_entries = False

    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if i != j:
                if matrix[i,j] == 0:
                    continue
                nonzero_entries = True
                if matrix[i,j] < lowest:
                    lowest = matrix[i,j]
                    row, col = i, j
    if not nonzero_entries:
        raise Exception("Sequences are identical. They should be further examined.")
    if row is None or col is None:
        raise Exception("Malformed input")

    b = []
    for i in range(len(matrix)):
        if i not in (row, col):
            b.append((matrix[row,i] + matrix[col,i]) / 2)

    new_matrix = np.delete(matrix, [row, col], axis=0)
    new_matrix = np.delete(new_matrix, [row, col], axis=1)

    size = len(new_matrix)+1
    WPGMA_matrix = np.zeros((size, size))
    WPGMA_matrix[1:,1:] = new_matrix
    WPGMA_matrix[1:,0] = b
    WPGMA_matrix[0,1:] = b
    WPGMA_matrix[0,0] = 0

    new_labels = merge_labels(labels, row, col)

    return WPGMA_matrix, new_labels


def Cluster(matrix: Union[list[list[float]], np.ndarray], labels: list[str]) -> str:
    """Perform WPGMA clustering until one cluster remains.
    Input : list(list(float)), list(string)
    Output : string (representing the tree)"""
    if not isinstance(matrix, (list, np.ndarray)) or not all(isinstance(l, str) for l in labels):
        raise Exception("malformed input")

    matrix = np.array(matrix, dtype=float)
    while len(labels) > 1:
        matrix, labels = reduction_matrix(matrix, labels)
    return labels[0]