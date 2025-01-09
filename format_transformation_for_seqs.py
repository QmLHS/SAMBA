def convert_sequence_format(input_text):
    """Converts seqs from NCBI format to FASTA-qiime2 format.
    
    Input format:
    ID \n Taxonomy \n Sequence \n
    Output format:
    >ID \n Sequence in uppercase
    """
    lines = input_text.strip().split('\n')
    output = []
    
    for i in range(0, len(lines), 3):
        # last uncomplete entry
        if i + 2 >= len(lines):
            break
            
        seq_id = lines[i].strip()
        sequence = lines[i + 2].strip().upper()
        
        output.append(f">{seq_id}")
        output.append(sequence)
    
    return '\n'.join(output)

def main():
    output_file = "ncbi_sequences_in.fasta"
    
    try:
        with open('ncbi_sequences.txt', 'r') as f:
            input_text = f.read()
        
        converted_text = convert_sequence_format(input_text)
        
        with open(output_file, 'w') as f:
            f.write(converted_text)
        print(f"Conversion successful. Output written to {output_file}")
    except FileNotFoundError:
        print(f"Error: Could not find input file ncbi_sequences.txt")
    except Exception as e:
        print(f"Error occurred: {str(e)}")

if __name__ == "__main__":
    main()
