from typing import Union
# for annotation of classes

def is_nucleic_acid(seq: str) -> bool:
    """
    Checks whether a given sequence contains only valid DNA or RNA nucleotides.

    Arguments:
        seq: str - nucleotide sequence

    Returns:
        bool - True if the sequence is a valid DNA or RNA sequence
    """
    seq_set = set(seq.upper())
    dna = {"A", "T", "G", "C"}
    rna = {"A", "U", "G", "C"}
    return seq_set.issubset(dna) or seq_set.issubset(rna)


def transcribe(seq: str) -> Union[str, None]:
    """
    Transcribes a DNA sequence into an RNA sequence.

    Arguments:
        seq: str - DNA sequence

    Returns:
        str - transcribed RNA sequence
        None - if the input is invalid or already RNA
    """
    if not is_nucleic_acid(seq):
        print(f"Error: invalid sequence '{seq}'")
        return None
    if "U" in seq.upper():
        print(f"Error: transcription is possible only for DNA, input: '{seq}'")
        return None
    table = str.maketrans({"T": "U", "t": "u"})
    return seq.translate(table)


def reverse(seq: str) -> Union[str, None]:
    """
    Returns the reversed nucleotide sequence.

    Arguments:
        seq: str - nucleotide sequence

    Returns:
        str - reversed sequence
        None - if the input is not a valid nucleic acid
    """
    if not is_nucleic_acid(seq):
        print(f"Error: invalid sequence '{seq}'")
        return None
    return seq[::-1]


def complement(seq: str) -> Union[str, None]:
    """
    Returns the complementary sequence of a DNA or RNA sequence.

    Arguments:
        seq: str - nucleotide sequence

    Returns:
        str - complementary sequence
        None - if the input is not a valid nucleic acid
    """
    if not is_nucleic_acid(seq):
        print(f"Error: invalid sequence '{seq}'")
        return None

    seq_upper = seq.upper()
    is_dna = 'U' not in seq_upper

    if is_dna:
        # DNA complement
        table = str.maketrans({
            'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
            'a': 't', 't': 'a', 'g': 'c', 'c': 'g'
        })
    else:
        # RNA complement
        table = str.maketrans({
            'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G',
            'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g'
        })

    return seq.translate(table)


def reverse_complement(seq: str) -> Union[str, None]:
    """
    Returns the reverse complement of a DNA or RNA sequence.

    Arguments:
        seq: str - nucleotide sequence

    Returns:
        str - reverse complement sequence
        None - if the input is not a valid nucleic acid
    """
    if not is_nucleic_acid(seq):
        print(f"Error: invalid sequence '{seq}'")
        return None

    comp = complement(seq)
    if comp is None:
        return None
    return reverse(comp)


def dna_rna_tools(*args: str) -> Union[str, list[str], None]:
    """
    Executes basic DNA/RNA utility functions such as transcription, reverse, and complement.

    Arguments:
        *args: variable number of arguments.
            The last argument must be the operation name (str),
            all preceding arguments are nucleotide sequences (str).

    Supported operations:
        - "is_nucleic_acid": check if the sequence is valid DNA/RNA
        - "transcribe": convert DNA to RNA
        - "reverse": reverse sequence
        - "complement": get complement sequence
        - "reverse_complement": get reverse complement sequence

    Returns:
        str or list[str] - result(s) of the chosen operation.
        None - if the operation is unknown or input is invalid.
    """
    if len(args) < 2:
        print("Ошибка: нужно передать хотя бы одну последовательность и операцию")
        return None

    *seqs, operation = args

    funcs = {
        "is_nucleic_acid": is_nucleic_acid,
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }

    if operation not in funcs:
        print(f"Ошибка: неизвестная операция '{operation}'")
        return None

    results = [funcs[operation](s) for s in seqs]
    # If all results are None (invalid inputs), return None
    if all(r is None for r in results):
        return None

    return results[0] if len(results) == 1 else results
