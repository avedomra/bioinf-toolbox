import os
import logging
import argparse
from abc import ABC, abstractmethod
from typing import Tuple, Union
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


# Logging setup

def setup_logger(log_file: str = "filter_fastq.log") -> logging.Logger:
    """Configure and return a logger that writes to *log_file*"""
    logger = logging.getLogger("bioinf_toolbox")
    logger.setLevel(logging.DEBUG)

    if not logger.handlers:
        fh = logging.FileHandler(log_file)
        fh.setLevel(logging.DEBUG)
        formatter = logging.Formatter(
            "%(asctime)s  %(levelname)-8s  %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


logger = setup_logger()


# Abstract base and sequence classes (unchanged)

class BiologicalSequence(ABC):
    """Abstract base class for all biological sequences"""

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
        """Check if sequence contains only valid characters for this type"""
        pass


class NucleicAcidSequence(BiologicalSequence):
    """
    Class for nucleic acid sequences (DNA/RNA)
    Not intended to be instantiated directly
    Subclasses must define ALPHABET and COMPLEMENT_MAP class attributes
    """

    ALPHABET: set = NotImplemented
    COMPLEMENT_MAP: dict = NotImplemented

    def __init__(self, sequence: str):
        if type(self) is NucleicAcidSequence:
            raise NotImplementedError(
                "NucleicAcidSequence is not intended to be instantiated directly "
                "Use DNASequence or RNASequence instead"
            )
        super().__init__(sequence)

    def check_alphabet(self) -> bool:
        """Check if all characters belong to the valid nucleotide alphabet"""
        return set(self._sequence.upper()).issubset(self.ALPHABET)

    def complement(self) -> "NucleicAcidSequence":
        """Return the complement sequence"""
        comp = self._sequence.translate(str.maketrans(self.COMPLEMENT_MAP))
        return self.__class__(comp)

    def reverse(self) -> "NucleicAcidSequence":
        """Return the reversed sequence"""
        return self.__class__(self._sequence[::-1])

    def reverse_complement(self) -> "NucleicAcidSequence":
        """Return the reverse complement sequence."""
        return self.complement().reverse()


class DNASequence(NucleicAcidSequence):
    """Class for DNA sequences"""

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
    """Class for RNA sequences"""

    ALPHABET = {"A", "U", "G", "C"}
    COMPLEMENT_MAP = {
        "A": "U", "U": "A", "G": "C", "C": "G",
        "a": "u", "u": "a", "g": "c", "c": "g",
    }


class AminoAcidSequence(BiologicalSequence):
    """Class for amino acid (protein) sequences"""

    ALPHABET = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy")

    def check_alphabet(self) -> bool:
        """Check if all characters are valid amino acid one-letter codes"""
        return set(self._sequence).issubset(self.ALPHABET)

    def amino_acid_percentage(self) -> dict[str, float]:
        """
        Calculate the percentage of each amino acid in the sequence

        Returns:
            dict[str, float] - amino acid one-letter codes as keys,
            rounded percentages as values (only present amino acids are included)
        """
        seq_upper = self._sequence.upper()
        total = len(seq_upper)
        result = {}
        for aa in set(seq_upper):
            result[aa] = round(seq_upper.count(aa) / total * 100, 2)
        return result


# FastQ filtering

def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[float, int, Tuple[Union[float, int], Union[float, int]]] = (0, 100),
    length_bounds: Union[float, int, Tuple[Union[float, int], Union[float, int]]] = (0, 2**32),
    quality_threshold: float = 0,
) -> None:
    """
    Filters FASTQ sequences by GC-content, length, and quality using Biopython

    Arguments:
        input_fastq: path to the input FASTQ file
        output_fastq: name of the output FASTQ file
        gc_bounds: GC-content interval in %, default (0, 100)
        length_bounds: sequence length interval, default (0, 2**32)
        quality_threshold: threshold value of mean phred quality, default 0

    Returns:
        None; writes a new FASTQ file inside the 'filtered' folder
        containing only sequences that pass all filters
    """
    logger.info(
        "Starting filter_fastq: input='%s', output='%s', "
        "gc_bounds=%s, length_bounds=%s, quality_threshold=%s",
        input_fastq, output_fastq, gc_bounds, length_bounds, quality_threshold,
    )

    gc_lower, gc_upper = _parse_bounds(gc_bounds, 0.0, 100.0)
    len_lower, len_upper = _parse_bounds(length_bounds, 0.0, float(2**32))

    os.makedirs("filtered", exist_ok=True)
    output_path = os.path.join("filtered", output_fastq)

    if os.path.exists(output_path):
        msg = f"File '{output_path}' already exists - choose another name"
        logger.error(msg)
        raise FileExistsError(msg)

    passed = 0
    total = 0
    with open(output_path, "w") as out_handle:
        for record in SeqIO.parse(input_fastq, "fastq"):
            total += 1
            if not _passes_filters(record, gc_lower, gc_upper,
                                   len_lower, len_upper, quality_threshold):
                continue
            SeqIO.write(record, out_handle, "fastq")
            passed += 1

    logger.info(
        "Filtering complete: %d / %d reads passed. Output written to '%s'.",
        passed, total, output_path,
    )


def _parse_bounds(
    bounds: Union[float, int, Tuple[Union[float, int], Union[float, int]]],
    default_lower: float = 0.0,
    default_upper: float = 100.0,
) -> Tuple[float, float]:
    """Parse interval bounds from a single number or a tuple of two numbers"""
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
    """Check whether a SeqRecord passes all three filters"""
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


# CLI entry point

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bioinf_toolbox",
        description="Filter FASTQ reads by GC content, length, and quality",
    )
    parser.add_argument(
        "input_fastq",
        help="Path to the input FASTQ file",
    )
    parser.add_argument(
        "output_fastq",
        help="Name of the output FASTQ file (saved inside the 'filtered/' folder)",
    )
    parser.add_argument(
        "--gc-lower",
        type=float,
        default=0.0,
        metavar="FLOAT",
        help="Lower bound for GC content in %% (default: 0)",
    )
    parser.add_argument(
        "--gc-upper",
        type=float,
        default=100.0,
        metavar="FLOAT",
        help="Upper bound for GC content in %% (default: 100)",
    )
    parser.add_argument(
        "--len-lower",
        type=int,
        default=0,
        metavar="INT",
        help="Minimum read length (default: 0)",
    )
    parser.add_argument(
        "--len-upper",
        type=int,
        default=2**32,
        metavar="INT",
        help="Maximum read length (default: 2^32)",
    )
    parser.add_argument(
        "--quality",
        type=float,
        default=0.0,
        metavar="FLOAT",
        help="Minimum mean Phred quality score (default: 0)",
    )
    parser.add_argument(
        "--log-file",
        type=str,
        default="filter_fastq.log",
        metavar="PATH",
        help="Path to the log file (default: filter_fastq.log)",
    )
    return parser


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()

    # Re-configure logger if a custom log file was requested
    global logger
    logger = setup_logger(args.log_file)

    filter_fastq(
        input_fastq=args.input_fastq,
        output_fastq=args.output_fastq,
        gc_bounds=(args.gc_lower, args.gc_upper),
        length_bounds=(args.len_lower, args.len_upper),
        quality_threshold=args.quality,
    )


if __name__ == "__main__":
    main()
