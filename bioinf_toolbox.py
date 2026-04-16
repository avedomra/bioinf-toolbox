import os
from abc import ABC, abstractmethod
from typing import Tuple, Union
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


class BiologicalSequence(ABC):
    """Abstract base class for all biological sequences."""

    def __init__(self, sequence: str):
        self._sequence = sequence

    def __len__(self) -> int:
        return len(self._sequence)

    def __getitem__(self, index):
        return self._sequence[index]

    def __str__(self) -> str:
        return f"{self.__class__.__name__}({self._sequence})"

    def __repr__(self) -> str:
        return self.__str__()

    @abstractmethod
    def check_alphabet(self) -> bool:
        """Check if sequence contains only valid characters for this type."""
        pass


class NucleicAcidSequence(BiologicalSequence):
    """
    Class for nucleic acid sequences (DNA/RNA).
    Not intended to be instantiated directly.
    Subclasses must define ALPHABET and COMPLEMENT_MAP class attributes.
    """

    ALPHABET: set = NotImplemented
    COMPLEMENT_MAP: dict = NotImplemented

    def __init__(self, sequence: str):
        if type(self) is NucleicAcidSequence:
            raise NotImplementedError(
                "NucleicAcidSequence is not intended to be instantiated directly. "
                "Use DNASequence or RNASequence instead."
            )
        super().__init__(sequence)

    def check_alphabet(self) -> bool:
        """Check if all characters belong to the valid nucleotide alphabet."""
        return set(self._sequence.upper()).issubset(self.ALPHABET)

    def complement(self) -> "NucleicAcidSequence":
        """Return the complement sequence."""
        comp = self._sequence.translate(str.maketrans(self.COMPLEMENT_MAP))
        return self.__class__(comp)

    def reverse(self) -> "NucleicAcidSequence":
        """Return the reversed sequence."""
        return self.__class__(self._sequence[::-1])

    def reverse_complement(self) -> "NucleicAcidSequence":
        """Return the reverse complement sequence."""
        return self.complement().reverse()


class DNASequence(NucleicAcidSequence):
    """Class for DNA sequences."""

    ALPHABET = {"A", "T", "G", "C"}
    COMPLEMENT_MAP = {
        "A": "T", "T": "A", "G": "C", "C": "G",
        "a": "t", "t": "a", "g": "c", "c": "g",
    }

    def transcribe(self) -> "RNASequence":
        """Transcribe DNA into RNA (T -> U)."""
        rna_seq = self._sequence.translate(str.maketrans({"T": "U", "t": "u"}))
        return RNASequence(rna_seq)


class RNASequence(NucleicAcidSequence):
    """Class for RNA sequences."""

    ALPHABET = {"A", "U", "G", "C"}
    COMPLEMENT_MAP = {
        "A": "U", "U": "A", "G": "C", "C": "G",
        "a": "u", "u": "a", "g": "c", "c": "g",
    }


class AminoAcidSequence(BiologicalSequence):
    """Class for amino acid (protein) sequences."""

    ALPHABET = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy")

    def check_alphabet(self) -> bool:
        """Check if all characters are valid amino acid one-letter codes."""
        return set(self._sequence).issubset(self.ALPHABET)

    def amino_acid_percentage(self) -> dict[str, float]:
        """
        Calculate the percentage of each amino acid in the sequence.

        Returns:
            dict[str, float] - amino acid one-letter codes as keys,
            rounded percentages as values (only present amino acids are included).
        """
        seq_upper = self._sequence.upper()
        total = len(seq_upper)
        result = {}
        for aa in set(seq_upper):
            result[aa] = round(seq_upper.count(aa) / total * 100, 2)
        return result




def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[float, int, Tuple[Union[float, int], Union[float, int]]] = (0, 100),
    length_bounds: Union[float, int, Tuple[Union[float, int], Union[float, int]]] = (0, 2**32),
    quality_threshold: float = 0,
) -> None:
    """
    Filters FASTQ sequences by GC-content, length, and quality using Biopython.

    Arguments:
        input_fastq: path to the input FASTQ file
        output_fastq: name of the output FASTQ file
        gc_bounds: GC-content interval in %, default (0, 100)
        length_bounds: sequence length interval, default (0, 2**32)
        quality_threshold: threshold value of mean phred quality, default 0

    Returns:
        None; writes a new FASTQ file inside the 'filtered' folder
        containing only sequences that pass all filters.
    """
    gc_lower, gc_upper = _parse_bounds(gc_bounds, 0.0, 100.0)
    len_lower, len_upper = _parse_bounds(length_bounds, 0.0, float(2**32))

    os.makedirs("filtered", exist_ok=True)
    output_path = os.path.join("filtered", output_fastq)

    if os.path.exists(output_path):
        raise FileExistsError(f"File '{output_path}' already exists — choose another name")

    with open(output_path, "w") as out_handle:
        for record in SeqIO.parse(input_fastq, "fastq"):
            if not _passes_filters(record, gc_lower, gc_upper,
                                   len_lower, len_upper, quality_threshold):
                continue
            SeqIO.write(record, out_handle, "fastq")


def _parse_bounds(
    bounds: Union[float, int, Tuple[Union[float, int], Union[float, int]]],
    default_lower: float = 0.0,
    default_upper: float = 100.0,
) -> Tuple[float, float]:
    """Parse interval bounds from a single number or a tuple of two numbers."""
    if isinstance(bounds, (int, float)):
        return default_lower, float(bounds)
    elif isinstance(bounds, tuple) and len(bounds) == 2:
        return float(bounds[0]), float(bounds[1])
    return default_lower, default_upper


def _passes_filters(
    record,
    gc_lower: float,
    gc_upper: float,
    len_lower: float,
    len_upper: float,
    quality_threshold: float,
) -> bool:
    """Check whether a SeqRecord passes all three filters."""
    seq_len = len(record.seq)
    if not (len_lower <= seq_len <= len_upper):
        return False

    gc_percent = gc_fraction(record.seq) * 100
    if not (gc_lower <= gc_percent <= gc_upper):
        return False

    qualities = record.letter_annotations["phred_quality"]
    mean_quality = sum(qualities) / len(qualities)
    if mean_quality < quality_threshold:
        return False

    return True