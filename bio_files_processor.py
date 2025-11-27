import os

def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str): 
    
    """
    Converts a multi-line FASTA format to one-line FASTA format.

    Arguments:
        input_fasta: path to the input FASTA file
        output_fasta: output file name (if not provided, will use the same name with '_oneline' suffix)
    """
    
    if not output_fasta:
        name, ext = os.path.splitext(input_fasta)
        output_fasta = f"{name}_oneline{ext if ext else '.fasta'}" 

    with open(input_fasta) as infile, open(output_fasta, "w") as outfile: 
        header = None
        seq = []

        for line in infile:
            line = line.strip()
            if not line:
                continue 

            if line.startswith(">"):
                if header:   
                    outfile.write(f"{header}\n{''.join(seq)}\n") 
                header = line
                seq = []
            else:
                seq.append(line) 

        if header:
            outfile.write(f"{header}\n{''.join(seq)}\n")

def parse_blast_output(input_file: str, output_file: str):
    
    """
    Parses BLAST results and extracts names of the best matches for each query.

    Arguments:
        input_file (str): path to the input file
        output_file (str): path to the output file
    """
    
    proteins = []
    
    with open(input_file, 'r') as infile:
        content = infile.read()
    
    queries = content.split('Query #')

    for query in queries[1:]:
        lines = query.split('\n') 
        alignment_section = False 
        best_match_found = False 

        for line in lines:
            line = line.strip()

            if "Sequences producing significant alignments:" in line:
                alignment_section = True
                continue

            if alignment_section and line and not best_match_found:
                if line.startswith('-') or not line:
                    continue

            parts = [part for part in line.split('  ') if part.strip()]
            if parts:
                if '[' in line:
                    prot_name = line.split('[')[0].strip()
                else: 
                    prot_name = line.split('  ')[0].strip() if '  ' in line else line
                if prot_name.endswith('...'):
                    prot_name = prot_name[:-3]
                    proteins.append(prot_name)
                    best_match_found = True
    
    unique_proteins = sorted(set(proteins))

    with open(output_file, 'w') as outfile:
        for protein in unique_proteins:
            outfile.write(protein + '\n')
