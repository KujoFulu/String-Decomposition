import pandas as pd
import glob
import os

# data = pd.read_csv('merged_new.tsv', sep='\t')
data = pd.read_csv('merged_original.tsv', sep='\t')

# Turn'Match_Percentage' to float
data['Match_Percentage'] = pd.to_numeric(data['Match_Percentage'], errors='coerce')

# Calculate the average of 'Match_Percentage' in 'Read_Name' 
average_match_percentage = data.groupby('Read_Name')['Match_Percentage'].mean().reset_index()

# print the result
print(average_match_percentage)

# Save result to tsv file
# average_match_percentage.to_csv('ave_match_percentage_new.tsv', sep='\t', index=False)
average_match_percentage.to_csv('ave_match_percentage_original.tsv', sep='\t', index=False)
