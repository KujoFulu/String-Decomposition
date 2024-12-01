import pandas as pd
data = pd.read_csv("final_decomposition_fc89af8.tsv", sep='\t', header=None)

# Assign proper column names based on the format provided
column_names = [
    "read-name", "best-monomer", "start-pos", "end-pos", "identity",
    "second-best-monomer", "second-best-monomer-identity", 
    "homo-best-monomer", "homo-identity", 
    "homo-second-best-monomer", "homo-second-best-monomer-identity", "reliability"
]

data.columns = column_names

# Convert 'identity' column to numeric to allow filtering
data['identity'] = pd.to_numeric(data['identity'], errors='coerce')

# Filter rows where 'identity' is less than 90
filtered_data = data[data['identity'] < 90]

output_file_path = 'filtered_rows_identity_less_than_90.csv'
filtered_data.to_csv(output_file_path, index=False)

# Count occurrences of each monomer in the 'best-monomer' column in the filtered data
monomer_counts = filtered_data['best-monomer'].value_counts()

# Display the results
monomer_counts

#B_1_DXZ1*_doubled/94_279/R'        11
#L_11_DXZ1*_doubled/1809_1977/R'    10
#D_3_DXZ1*_doubled/451_620/R'       10
#E_4_DXZ1*_doubled/621_788/R'        9
#F_5_DXZ1*_doubled/789_959/R'        9
#G_6_DXZ1*_doubled/960_1129/R'       5
#A_0_DXZ1*_doubled/1978_2147/R'      5
#H_7_DXZ1*_doubled/1130_1300/R'      5
#I_8_DXZ1*_doubled/1301_1470/R'      5
#J_9_DXZ1*_doubled/1471_1638/R'      4
#K_10_DXZ1*_doubled/1639_1808/R'     4
#C_2_DXZ1*_doubled/280_450/R'        4

# Filter rows where 'identity' is less than 80
filtered_data = data[data['identity'] < 80]

output_file_path = 'filtered_rows_identity_less_than_80.csv'
filtered_data.to_csv(output_file_path, index=False)

# Count occurrences of each monomer in the 'best-monomer' column in the filtered data
monomer_counts = filtered_data['best-monomer'].value_counts()

# Display the results
monomer_counts

#F_5_DXZ1*_doubled/789_959/R'      3
#J_9_DXZ1*_doubled/1471_1638/R'    1
#H_7_DXZ1*_doubled/1130_1300/R'    1
#E_4_DXZ1*_doubled/621_788/R'      1
#D_3_DXZ1*_doubled/451_620/R'      1
#C_2_DXZ1*_doubled/280_450/R'      1
#B_1_DXZ1*_doubled/94_279/R'       1

k = "ACCGTCTGGTTTTTATATGAAGTTCTTTCCTTCACTACCACAGGCCTCAAAGCGGTCCAAATCTCCACTTGCAGATTCTACAAAAAGAGTGTTTGCAAACTGCTCTATCAAAGGAATGTTCAACTCTGGGAGTTGAATGCAATCATCACAGAGCAGTTTCTGAGAATGCT"
f = "TCCGTTTGCCTTTTATATGAAGTTCCTTCCTGTACTACCGTAGGCCTCAAAGCAGTCCAAATCTCCATTTGCAGATTCTACAAAAAGAGTGATTCCAATCTGCTCTATCAATAGGATTGTTCAACTCCATGAGTTGAATGCCATCCTCACAAAGTCGTTTCTGAGAATGCT"

#Build the novel monomer mentioned in the paper
print(k[0:70]+f[len(f)-101:])
#ACCGTCTGGTTTTTATATGAAGTTCTTTCCTTCACTACCACAGGCCTCAAAGCGGTCCAAATCTCCACTTGCAGATTCTACAAAAAGAGTGATTCCAATCTGCTCTATCAATAGGATTGTTCAACTCCATGAGTTGAATGCCATCCTCACAAAGTCGTTTCTGAGAATGCT

#Since there are 3 F regions with low ideneity, I saved them in lowIdentity_F5_regions_reads.fasta
# Load the filtered rows with identity < 80
filtered_data_80 = pd.read_csv('filtered_rows_identity_less_than_80.csv')

# Filter rows for F_5 monomer
f5_regions = filtered_data_80[filtered_data_80['best-monomer'] == "F_5_DXZ1*_doubled/789_959/R'"]

# Read the sequence from the second line of the read.fa file
with open('read.fa', 'r') as fasta_file:
    fasta_content = fasta_file.readlines()
sequence = fasta_content[1].strip()

# Extract the regions and save to a new FASTA file
f5_reads = []
for _, row in f5_regions.iterrows():
    start_pos = int(row['start-pos'])
    end_pos = int(row['end-pos'])
    f5_reads.append(f">{row['read-name']}\n{sequence[start_pos:end_pos]}")

with open('lowIdentity_F5_regions_reads.fasta', 'w') as f:
    f.write("\n".join(f5_reads))