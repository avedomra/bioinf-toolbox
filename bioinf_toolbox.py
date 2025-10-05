from typing import Dict, Tuple, Union
# for annotation of classes
from fastq_tools import parse_bounds, check_sequence_validity
# for filter_fastq function 
from dna_rna_tools import is_nucleic_acid, transcribe, reverse, complement, reverse_complement
# for dna_rna_tools function 


def filter_fastq(seqs: Dict[str, Tuple[str, str]],
                 gc_bounds: Union[float, int, Tuple[Union[float, int], Union[float, int]]] = (0, 100),
                 length_bounds: Union[float, int, Tuple[Union[float, int], Union[float, int]]] = (0, 2**32),
                 quality_threshold: float = 0) -> Dict[str, Tuple[str, str]]:
    """
    Filters FASTQ sequences by GC-content, length, and quality.
    
    Arguments:
        seqs: dictionary of FASTQ sequences, where the key is the sequence ID, 
        the value is a tuple [nucleotide sequence, quality string]
        gc_bounds: GC-content interval in %, default (0, 100)
        length_bounds: sequence length interval, default (0, 2**32)
        quality_threshold: threshold value of mean quality, default 0
        
    Returns:
        Dict[str, Tuple[str, str]]: filtered dictionary of sequences that satisfy all requirements
    """
    
    parsed_gc_bounds = parse_bounds(gc_bounds, 0.0, 100.0)
    parsed_length_bounds = parse_bounds(length_bounds, 0.0, 2**32)
    
    filtered_seqs = {}
    
    for seq_id, (sequence, quality) in seqs.items():
        if check_sequence_validity(sequence, quality, parsed_gc_bounds, 
                                parsed_length_bounds, quality_threshold):
            filtered_seqs[seq_id] = (sequence, quality)
    
    return filtered_seqs

def dna_rna_tools(*args: str) -> Union[str, list[str], None]:
    """
    Executes basic DNA/RNA utility functions such as transcription, reverse, and complement.

    Arguments:
        *args: variable number of arguments
            The last argument must be the procedure name (str),
            all preceding arguments are nucleotide sequences (str).

    Supported procedures:
        - "is_nucleic_acid": check if the sequence is valid DNA/RNA
        - "transcribe": convert DNA to RNA
        - "reverse": reverse sequence
        - "complement": get complement sequence
        - "reverse_complement": get reverse complement sequence

    Returns:
        The result of the chosen operation.
        Returns a single value if one sequence was passed,
        or a list of results if multiple sequences were provided.
    """

    if len(args) < 2:
        print("Error: Please enter at least one sequence and procedure")
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
        print(f"Error: uncnown operation '{operation}'")
        return None

    results = [funcs[operation](s) for s in seqs]
    # If all results are None (invalid inputs), return None
    if all(r is None for r in results):
        return None

    return results[0] if len(results) == 1 else results
