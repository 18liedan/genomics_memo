import sys

def extract_masked_ranges(fasta_path):
    with open(fasta_path, 'r') as f:
        chrom = None
        pos = 0
        in_masked_block = False
        start_pos = 0
        
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith(">"):
                # If we were in a masked block at the end of the previous chromosome
                if in_masked_block:
                    print(f"{chrom}\t{start_pos}\t{pos}")
                
                # Reset for new chromosome
                chrom = line[1:].split()[0]
                pos = 0
                in_masked_block = False
            else:
                for char in line:
                    # Check if character is lowercase (soft-masked)
                    if char.islower():
                        if not in_masked_block:
                            start_pos = pos
                            in_masked_block = True
                    else:
                        if in_masked_block:
                            print(f"{chrom}\t{start_pos}\t{pos}")
                            in_masked_block = False
                    pos += 1
        
        # Catch the last block if the file ends on a masked region
        if in_masked_block:
            print(f"{chrom}\t{start_pos}\t{pos}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        extract_masked_ranges(sys.argv[1])
