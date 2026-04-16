# Bioinf Toolbox :beetle:

A lightweight set of tools for basic DNA, RNA, protein, and FASTQ sequence analysis.  
Includes OOP-based utilities for sequence validation, transcription, complement generation, and FASTQ filtering.

---

## Installation

Clone the repository and navigate to the project folder:

```bash
git clone https://github.com/avedomra/bioinf-toolbox
cd bioinf-toolbox
```

Create the conda environment and activate it:

```bash
conda env create -f environment.yml
conda activate bioinf_toolbox
```

---

## Main script `bioinf_toolbox.py`

The primary script that integrates biological sequence classes and FASTQ filtering capabilities.

**Key classes:**

- `BiologicalSequence` - abstract base class defining the interface for all biological sequences
- `NucleicAcidSequence` - base class for nucleic acids with complement, reverse, and reverse complement methods
- `DNASequence` - DNA sequence class with transcription support
- `RNASequence` - RNA sequence class
- `AminoAcidSequence` - protein sequence class

**Key functions:**

- `filter_fastq()` - FASTQ filtering function powered by Biopython

---

## Classes

### `BiologicalSequence`

Abstract base class. Defines the shared interface for all biological sequences:

- `len(seq)` - sequence length
- `seq[i]`, `seq[i:j]` - indexing and slicing
- `str(seq)` - human-readable representation
- `check_alphabet()` - validates the sequence alphabet (abstract)

### `NucleicAcidSequence`

Base class for DNA and RNA sequences. Not intended to be instantiated directly - use `DNASequence` or `RNASequence` instead.

Available methods:

- `check_alphabet()` - check if all characters belong to the valid nucleotide alphabet
- `complement()` - return the complement sequence
- `reverse()` - return the reversed sequence
- `reverse_complement()` - return the reverse complement sequence

### `DNASequence`

Inherits from `NucleicAcidSequence`. Valid alphabet: `A, T, G, C`.

Additional method:

- `transcribe()` - transcribe DNA into an `RNASequence` object (T → U)

### `RNASequence`

Inherits from `NucleicAcidSequence`. Valid alphabet: `A, U, G, C`.

### `AminoAcidSequence`

Implements the `BiologicalSequence` interface for protein sequences. Valid alphabet: standard 20 amino acid one-letter codes.

Available methods:

- `check_alphabet()` - check if all characters are valid amino acid one-letter codes
- `amino_acid_percentage()` - return a dictionary with the percentage of each amino acid present in the sequence

---

## Usage

### Biological Sequences

```python
from bioinf_toolbox import DNASequence, RNASequence, AminoAcidSequence

# DNA
dna = DNASequence("ATGC")
print(dna.complement())         # DNASequence(TACG)
print(dna.reverse())            # DNASequence(CGTA)
print(dna.reverse_complement()) # DNASequence(GCAT)
print(dna.transcribe())         # RNASequence(AUGC)
print(dna.check_alphabet())     # True
print(len(dna))                 # 4
print(dna[1:3])                 # TG

# RNA
rna = RNASequence("AUGC")
print(rna.complement())         # RNASequence(UACG)

# Protein
prot = AminoAcidSequence("AACDEELLLL")
print(prot.amino_acid_percentage())
# {'A': 20.0, 'C': 10.0, 'D': 10.0, 'E': 20.0, 'L': 40.0}
```

### FASTQ Filtering - Python API

Filters FASTQ sequences by GC-content, length, and quality using Biopython.  
Output is written inside the `filtered` folder.

Available parameters:

- `gc_bounds` - required GC-content interval (in %), default `(0, 100)`
- `length_bounds` - required sequence length interval, default `(0, 2**32)`
- `quality_threshold` - minimum required mean quality (Phred33), default `0`

```python
from bioinf_toolbox import filter_fastq

filter_fastq("input.fastq", "output.fastq",
             gc_bounds=(40, 60),
             length_bounds=(256, 4294967296),
             quality_threshold=30)

# Single value means upper bound only
filter_fastq("input.fastq", "output.fastq",
             gc_bounds=50,
             length_bounds=2**32,
             quality_threshold=25)
```

### FASTQ Filtering - Command Line Interface

`filter_fastq` can also be run directly from the command line:

```bash
python bioinf_toolbox.py input.fastq output.fastq \
    --gc-lower 40 --gc-upper 60 \
    --len-lower 100 --len-upper 500 \
    --quality 30 \
    --log-file my_run.log
```

All arguments:

| Argument | Type | Default | Description |
|---|---|---|---|
| `input_fastq` | positional | - | Path to the input FASTQ file |
| `output_fastq` | positional | - | Name of the output file (saved in `filtered/`) |
| `--gc-lower` | float | `0` | Lower bound for GC content (%) |
| `--gc-upper` | float | `100` | Upper bound for GC content (%) |
| `--len-lower` | int | `0` | Minimum read length |
| `--len-upper` | int | `2^32` | Maximum read length |
| `--quality` | float | `0` | Minimum mean Phred quality score |
| `--log-file` | path | `filter_fastq.log` | Path to the log file |

### Logging

All runs log to a file (default: `filter_fastq.log`). Each run records:

- `INFO` - start of filtering with all parameters
- `INFO` - number of reads passed and path to the output file
- `ERROR` - if the output file already exists

Example log output:

```
2026-04-16 12:00:00  INFO      Starting filter_fastq: input='input.fastq', output='out.fastq', ...
2026-04-16 12:00:01  INFO      Filtering complete: 142 / 200 reads passed. Output written to 'filtered/out.fastq'.
```

---

## Tests

Tests are located in `tests/test_filter_fastq.py` and use `pytest`.

To run:

```bash
pytest tests/test_filter_fastq.py -v
```

The test suite covers:

- `_parse_bounds` - single value and tuple inputs
- `_passes_filters` - reads passing all filters, and failing by GC, length, or quality
- `filter_fastq` - output file creation, correct read selection, `FileExistsError` on duplicate output
- logging - log file creation, error logging on duplicate output

---

## Author & Contacts

Masha Domracheva

:frog: email: <m.domracheva2000@yandex.ru>

:snake: GitHub: [avedomra](https://gist.github.com/avedomra)
