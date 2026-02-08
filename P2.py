def sigma(a, b):
    """Define match/mismatch scores after comparing two characters."""
    match, mismatch = 5, -2
    return match if a == b else mismatch


def get_needleman_wunsch(seq1, seq2):
    """Build the Needleman–Wunsch DP matrix for two sequences."""
    match, mismatch, gap = 5, -2, -6
    n, m = len(seq1), len(seq2)
    S = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        S[i][0] = S[i - 1][0] + gap
    for j in range(1, m + 1):
        S[0][j] = S[0][j - 1] + gap

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = S[i - 1][j - 1] + sigma(seq1[i - 1], seq2[j - 1])
            up   = S[i - 1][j] + gap
            left = S[i][j - 1] + gap
            S[i][j] = max(diag, up, left)

    return S


def get_traceback(seq1, seq2):
    """Trace back through the matrix to build the aligned sequences.
    It returns a pair of aligned strings."""
    match, mismatch, gap = 5, -2, -6
    S = get_needleman_wunsch(seq1, seq2)
    aligned_seq1, aligned_seq2 = "", ""
    i, j = len(seq1), len(seq2)

    while i > 0 or j > 0:
        if i > 0 and j > 0 and S[i][j] == S[i-1][j-1] + sigma(seq1[i-1], seq2[j-1]):
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1

        elif i > 0 and S[i][j] == S[i-1][j] + gap:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1

        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1

    return (aligned_seq1, aligned_seq2)


def AlignByDP(list_of_sequences):
    """Align pairwise all pairs of sequences using Needleman–Wunsch.
    Input: list(tuple(str, str)).
    Returns: dict(tuple(int, int) -> tuple(string, string))"""
    if not isinstance(list_of_sequences, list) or not all(isinstance(p, tuple) and len(p) == 2 and isinstance(p[1], str) for p in list_of_sequences):
        raise Exception("malformed input")

    output_dict = {}
    n = len(list_of_sequences)

    for i in range(n):
        for j in range(i + 1, n):
            seq_i = list_of_sequences[i][1]
            seq_j = list_of_sequences[j][1]

            alignment = get_traceback(seq_i, seq_j)

            output_dict[(i, j)] = alignment

    return output_dict