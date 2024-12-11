import pandas as pd
import glob
import os

# data = pd.read_csv('merged_new.tsv', sep='\t')
data = pd.read_csv('merged_original.tsv', sep='\t')

# 将'Match_Percentage'列转换为浮点数类型（如果不是的话）
data['Match_Percentage'] = pd.to_numeric(data['Match_Percentage'], errors='coerce')

# 计算每个'Read_Name'的'Match_Percentage'均值
average_match_percentage = data.groupby('Read_Name')['Match_Percentage'].mean().reset_index()

# 输出结果
print(average_match_percentage)

# 如果你想将结果保存到一个新的TSV文件中
# average_match_percentage.to_csv('ave_match_percentage_new.tsv', sep='\t', index=False)
average_match_percentage.to_csv('ave_match_percentage_original.tsv', sep='\t', index=False)
