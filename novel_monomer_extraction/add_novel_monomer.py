def create_new_fasta_with_monomer(input_fasta, output_fasta, new_sequence):
    # Read existing content
    with open(input_fasta, 'r') as f:
        content = f.read()
    
    # Create new monomer entry
    new_entry = f">Z_12_DXZ1*_doubled/0000_0000/R\n{new_sequence}\n"
    
    # Write all content to new file
    with open(output_fasta, 'w') as f:
        f.write(content)
        f.write(new_entry)

# Create the new monomer sequence
new_monomer = 'ACCGTCTGGTTTTTATATGAAGTTCTTTCCTTCACTACCACAGGCCTCAAAGCGGTCCAAATCTCCACTTGCAGATTCTACAAAAAGAGTGATTCCAATCTGCTCTATCAATAGGATTGTTCAACTCCATGGGTTGAATGCCATCCTCACAAATTAGTTTCTGAGAATGC'

# Create new file
create_new_fasta_with_monomer("../test_data/DXZ1_star_monomers.fa", "../test_data/DXZ1_star_monomers_new.fa", new_monomer)