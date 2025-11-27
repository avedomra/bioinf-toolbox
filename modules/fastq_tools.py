from typing import Tuple, Union
# for annotation of classes in functions

def calculate_gc_content(sequence: str) -> float:
    """
    Calsulates GC-content of sequence in %. 
    
    Arguments:
        sequence: str - nucleotide sequence
        
    Returns:
        float - GC-content in % (0-100)
    """
    sequence = sequence.upper()
    if not sequence: 
        return 0.0
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100


def calculate_mean_quality(quality_string: str) -> float:
    """
    Calculates mean quality of read in phred33 quality score.
    
    Arguments:
        quality_string: str - string of quality scores in phred33
        
    Returns:
        float - mean quality score in phred33
    """
    if not quality_string:
        return 0.0
    quality_scores = [ord(char) - 33 for char in quality_string] # ord - returns numeric code of ASCII character 
    return sum(quality_scores) / len(quality_scores)


def parse_bounds(bounds: Union[float, int, Tuple[Union[float, int], Union[float, int]]], # set all possible input formats
                default_lower: float = 0.0, 
                default_upper: float = 100.0) -> Tuple[float, float]: 
    """
    Parses interval bounds from various input formats.
    
    Arguments:
        bounds: one number (int or float) or tuple of two numbers (int or float)
        default_lower: default lower boundary
        default_upper: default upper boundary
        
    Returns:
        Tuple[float, float] - tuple of two numbers [lower, upper]
    """
    if isinstance(bounds, (int, float)):
        return default_lower, float(bounds)
    elif isinstance(bounds, tuple) and len(bounds) == 2:
        return float(bounds[0]), float(bounds[1])
    else:
        return default_lower, default_upper


def check_sequence_validity(sequence: str, quality: str, 
                          gc_bounds: Tuple[float, float],
                          length_bounds: Tuple[float, float],
                          quality_threshold: float) -> bool:
    """
    Checks that the sequence satisfies all requirements.
    
    Arguments:
        sequence: str - nucleotide sequence 
        quality: string of quality scores in phred33
        gc_bounds: Tuple[float, float] - GC-content bounds
        length_bounds: Tuple[float, float] - length bounds
        quality_threshold: float - threshold value of mean quality
        
    Returns:
        bool - True if the sequence passes all filters
    """
    # check length
    seq_length = len(sequence)
    lower_length, upper_length = length_bounds
    if not (lower_length <= seq_length <= upper_length):
        return False
    
    # check GC-content
    gc_content = calculate_gc_content(sequence)
    lower_gc, upper_gc = gc_bounds
    if not (lower_gc <= gc_content <= upper_gc):
        return False
    
    # check quality
    mean_quality = calculate_mean_quality(quality)
    if mean_quality < quality_threshold:
        return False
    
    return True
