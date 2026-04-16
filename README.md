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

- `BiologicalSequence` — abstract base class defining the interface for all biological sequences
- `NucleicAcidSequence` — base class for nucleic acids with complement, reverse, and reverse complement methods
- `DNASequence` — DNA sequence class with transcription support
- `RNASequence` — RNA sequence class
- `AminoAcidSequence` — protein sequence class

**Key functions:**

- `filter_fastq()` — FASTQ filtering function powered by Biopython

---

## Classes

### `BiologicalSequence`

Abstract base class. Defines the shared interface for all biological sequences:

- `len(seq)` — sequence length
- `seq[i]`, `seq[i:j]` — indexing and slicing
- `str(seq)` — human-readable representation
- `check_alphabet()` — validates the sequence alphabet (abstract)

### `NucleicAcidSequence`

Base class for DNA and RNA sequences. Not intended to be instantiated directly — use `DNASequence` or `RNASequence` instead.

Available methods:

- `check_alphabet()` — check if all characters belong to the valid nucleotide alphabet
- `complement()` — return the complement sequence
- `reverse()` — return the reversed sequence
- `reverse_complement()` — return the reverse complement sequence

### `DNASequence`

Inherits from `NucleicAcidSequence`. Valid alphabet: `A, T, G, C`.

Additional method:

- `transcribe()` — transcribe DNA into an `RNASequence` object (T → U)

### `RNASequence`

Inherits from `NucleicAcidSequence`. Valid alphabet: `A, U, G, C`.

### `AminoAcidSequence`

Implements the `BiologicalSequence` interface for protein sequences. Valid alphabet: standard 20 amino acid one-letter codes.

Available methods:

- `check_alphabet()` — check if all characters are valid amino acid one-letter codes
- `amino_acid_percentage()` — return a dictionary with the percentage of each amino acid present in the sequence

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

### FASTQ Filtering

Filters FASTQ sequences by GC-content, length, and quality using Biopython.  
Output is written inside the `filtered` folder.

Available parameters:

- `gc_bounds` — required GC-content interval (in %), default `(0, 100)`
- `length_bounds` — required sequence length interval, default `(0, 2**32)`
- `quality_threshold` — minimum required mean quality (Phred33), default `0`

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

---

## Author & Contacts

Masha Domracheva

:frog: email: <m.domracheva2000@yandex.ru>

:snake: GitHub: [avedomra](https://gist.github.com/avedomra)