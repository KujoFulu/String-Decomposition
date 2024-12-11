import os
import glob

input_folder = '100_new_monomers'  # 请替换为你的实际路径
# input_folder = 'newalgo_new_monomer' 
output_file = 'merged_new.tsv'

# 获取文件夹下所有的TSV文件列表
tsv_files = glob.glob(os.path.join(input_folder, '*.tsv'))

# 初始化一个标记，判断是否已经写入表头
header_saved = False

with open(output_file, 'w', encoding='utf-8') as outfile:
    for filename in tsv_files:
        with open(filename, 'r', encoding='utf-8') as infile:
            # 读取文件的所有行
            lines = infile.readlines()
            # 如果文件不为空
            if lines:
                if not header_saved:
                    # 写入表头
                    outfile.write(lines[0])
                    header_saved = True
                # 写入数据（跳过表头）
                outfile.writelines(lines[1:])
