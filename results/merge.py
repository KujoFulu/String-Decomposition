import os
import glob

input_folder = '100_new_monomers'
# input_folder = 'newalgo_new_monomer' 
output_file = 'merged_new.tsv'

# Get all tsv file path under this folder
tsv_files = glob.glob(os.path.join(input_folder, '*.tsv'))

# Initialize a marker
header_saved = False

with open(output_file, 'w', encoding='utf-8') as outfile:
    for filename in tsv_files:
        with open(filename, 'r', encoding='utf-8') as infile:
            # Read all rows
            lines = infile.readlines()
            # if not empty
            if lines:
                if not header_saved:
                    # write head
                    outfile.write(lines[0])
                    header_saved = True
                # write data
                outfile.writelines(lines[1:])
