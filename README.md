# Bioinf Toolbox :beetle:

A lightweight set of tools for basic DNA and RNA sequence analysis.  
Includes utilities for sequence validation, transcription, complement generation, and FASTQ filtering.

---

## Installation

Clone the repository and navigate to the project folder:

```bash
git clone https://github.com/avedomra/bioinf-toolbox
cd bioinf-toolbox 
```

## Main program `bioinf_toolbox.py`

The primary interface that integrates both DNA/RNA tools and FASTQ filtering capabilities.

Key functions:

- `dna_rna_tools()` - unified interface for DNA/RNA sequence operations
- `filter_fastq()` - FASTQ filtering function

## Modules

### `dna_rna_tools.py`

Core utilities for nucleic acid sequence analysis and transformation.

Available functions:

- `is_nucleic_acid`: check if the sequence is valid DNA/RNA
- `transcribe`: convert DNA to RNA
- `reverse`: reverse sequence
- `complement`: get complement sequence
- `reverse_complement`: get reverse complement sequence

### `fastq_tools.py`

Utilities for filtering and analyzing FASTQ sequencing data.

Available functions:

- `calculate_gc_content()`: computes GC-content percentage
- `calculate_mean_quality()`: calculates average Phred33 quality scores
- `parse_bounds()` - handles flexible boundary input formats
- `check_sequence_validity()` - validates sequences against multiple criteria

## Usage

### DNA/RNA Tools

Executes basic DNA/RNA utility functions such as transcription, reverse, and complement.

```python
from bioinf_toolbox import dna_rna_tools

print(dna_rna_tools("ATGC", "reverse"))
print(dna_rna_tools("ATGC", "complement"))
print(dna_rna_tools("ATGC", "reverse_complement"))
print(dna_rna_tools("ATGC", "transcribe"))
```

### FASTQ Filtering

Filters FASTQ sequences by GC-content, length, and quality.

Available parameters:

- `gc_bounds`: required GC-content interval (in %)
- `length_bounds`: required sequence length interval
- `quality_threshold`: minimum required mean quality (in Phred33 encoded)

```python
from bioinf_toolbox import filter_fastq

filtered_1 = filter_fastq(your_fastq, 
             gc_bounds=(40, 60), 
             length_bounds=(256, 4294967296), 
             quality_threshold=30)

filtered_2 = filter_fastq(your_fastq, 
             gc_bounds=(50), 
             length_bounds=(2**32), 
             quality_threshold=25)
```

## Advanced Usage

You can also directly import modules for more detailed management.

```python
# Direct module imports
from dna_rna_tools import transcribe
from fastq_tools import calculate_gc_content

# Use individual functions
rna_seq = transcribe("ATGC")
gc_percent = calculate_gc_content("AGATACACA")
```

## Author & Contacts

Maria Domracheva

:frog: email: <m.domracheva2000@yandex.ru>
:snake: GitHub: [avedomra](https://gist.github.com/avedomra)

