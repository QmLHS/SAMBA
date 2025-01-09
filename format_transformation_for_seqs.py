def convert_sequence_format(input_text):
    """Converts seqs from NCBI format to FASTA-qiime2 format.
    
    Input format:
    ID \n Taxonomy \n Sequence \n
    Output format:
    >ID \n Fixed-width sequence with dashes for padding
    """
    lines = input_text.strip().split('\n')
    output = []
    
    for i in range(0, len(lines), 3):
        # last uncomplete entry
        if i + 2 >= len(lines):
            break
            
        seq_id = lines[i].strip()
        # taxo = lines[i + 1].strip()
        sequence = lines[i + 2].strip()
        
        output.append(f">{seq_id}")
        
        block_size = 150
        sequence_blocks = []
        for j in range(0, len(sequence), block_size):
            block = sequence[j:j + block_size]
            padded_block = block.ljust(block_size, '-')
            sequence_blocks.append(padded_block)
        output.extend(sequence_blocks)
        output.append('')
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
