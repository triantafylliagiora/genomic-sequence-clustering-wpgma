from typing import Tuple, Set


def check_line_empty_or_whitespace(line: str) -> bool:
    """Returns True if the line is empty or has whitespace."""
    return not line or not line.strip()


def check_line_starts_with_arrow(line: str) -> None:
    """Make sure the line begins with '>' otherwise raise an error."""
    stripped = line.lstrip()
    if not stripped.startswith('>'):
        raise Exception("malformed input")


def extract_label_and_sequence_parts(content: str) -> Tuple[str, str]:
    """Split the line (after '>') into a label and the sequence.
    The label must come first followed by at least one whitespace and the sequence.
    Raise an error if the format is incorrect."""
    parts = content.split(None, 1)
    if len(parts) < 2:
        raise Exception("malformed input")
    return parts[0], parts[1]


def validate_label_not_empty(label: str) -> None:
    """Raise an exception if the label is empty."""
    if not label:
        raise Exception("malformed input")


def validate_sequence_not_empty(sequence: str) -> None:
    """Raise an exception if the sequence is empty."""
    if not sequence:
        raise Exception("malformed input")

def normalize_sequence_whitespace(sequence: str) -> str:
    """Remove all spaces and tabs from the sequence and convert it to uppercase."""
    sequence = sequence.replace(' ', '').replace('\t', '')
    return sequence.upper()


def validate_sequence_characters(sequence: str) -> None:
    """Check that every character in the sequence is one of A, C, G, T.
    Raise an error if any invalid character is found."""
    for char in sequence:
        if char not in 'ACGT':
            raise Exception("malformed input")


def check_duplicates_label_sequence(label: str, sequence: str, seen_labels: Set[str], seen_sequences: Set[str]) -> None:
    """Ensure that the label unique or raise an error if it has already been seen."""
    if label in seen_labels:
        raise Exception("malformed input")
    else:
        seen_labels.add(label)


def ParseSeqFile(fasta: str) -> list[Tuple[str, str]]:
    """Parse a fasta file where each line must contain:
        >label sequence
    The function validates the formatting, the uniqueness and the allowed characters.
    It returns list(tuple(str, str))."""
    try:
        with open(fasta, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        raise Exception("malformed input")

    parsed_list = []
    seen_labels = set()
    seen_sequences = set()

    for line in lines:
        line = line.rstrip('\n')

        if check_line_empty_or_whitespace(line):
            continue

        check_line_starts_with_arrow(line)
        content = line[1:].lstrip()

        label, seq_part = extract_label_and_sequence_parts(content)
        validate_label_not_empty(label)

        sequence = normalize_sequence_whitespace(seq_part)
        validate_sequence_not_empty(sequence)
        validate_sequence_characters(sequence)
        check_duplicates_label_sequence(label, sequence, seen_labels, seen_sequences)

        parsed_list.append((label, sequence))

    return parsed_list
