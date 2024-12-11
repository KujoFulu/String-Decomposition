import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# helper function to preprocess the data
# for each file:
# 1. take only the first column
# 3. if it is "*_rc", remove the "_rc" part and reverse the order of the list
def preprocess_data(file):
    data = pd.read_csv(file, sep='\t', header=None)
    first_column = data[0].tolist()
    return first_column[1:]

# helper function to process all files in a list
# monoseq_list = [[A, B, C,...], [A_rc, B_rc, C_rc,...], ...]
def process_all_to_lists(file_list):
    monoseq_list = []
    for file in file_list:
        monoseq_list.append(preprocess_data(file))
    return monoseq_list

# draw a heatmap of the monomer sequences
# 1. x-axis: 12 monomers, y-axis: 12 monomers
# 2. each cell represents the count of monomer x following monomer y
def draw_heatmap_original(monoseq_list):
    # define the 24 monomers (A to L; A_rc to L_rc)
    monomers = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L',
                'A_rc', 'B_rc', 'C_rc', 'D_rc', 'E_rc', 'F_rc', 'G_rc', 'H_rc', 'I_rc', 'J_rc', 'K_rc', 'L_rc']

    # create a dictionary to store the counts of each pair of monomers
    monomer_dict = {}
    for monoseq in monoseq_list:
        for i in range(len(monoseq) - 1):
            pair = monoseq[i] + '+' + monoseq[i + 1]
            if pair not in monomer_dict:
                monomer_dict[pair] = 1
            else:
                monomer_dict[pair] += 1

    # create a matrix to store the counts
    count_matrix = np.zeros((24, 24), dtype=int)

    # fill in the matrix
    for pair, count in monomer_dict.items():
        monomer_x, monomer_y = pair.split('+')
        if monomer_x in monomers and monomer_y in monomers:
            x = monomers.index(monomer_x)
            y = monomers.index(monomer_y)
            count_matrix[x, y] = count

    # plot the heatmap
    plt.figure(figsize=(12, 10))
    plt.imshow(count_matrix, cmap='viridis', interpolation='nearest')
    plt.colorbar(label='Count')
    # adjust ticks so that labels are centered in each cell
    plt.xticks(ticks=np.arange(24), labels=monomers, rotation=45, ha='right')
    plt.yticks(ticks=np.arange(24), labels=monomers)
    # draw gridlines to align with cells
    plt.gca().set_xticks(np.arange(-0.5, 24, 1), minor=True)
    plt.gca().set_yticks(np.arange(-0.5, 24, 1), minor=True)
    plt.grid(which='minor', color='white', linestyle='-', linewidth=0.5)
    plt.tick_params(which='minor', size=0)  # hide minor ticks
    # annotate each cell with its count
    for i in range(24):
        for j in range(24):
            count = count_matrix[i, j]
            if count > 0:  # only annotate cells with non-zero counts
                color = 'black' if count > 1000 else 'white'
                plt.text(j, i, str(count), ha='center', va='center', color=color, fontsize=8)
    # set labels and title
    plt.xlabel('Monomer X')
    plt.ylabel('Monomer Y')
    plt.title('Count Heatmap of Monomer Sequences (x follows y)')
    plt.show()

# draw a heatmap of the monomer sequences
# 1. x-axis: 13 monomers, y-axis: 13 monomers
# 2. each cell represents the count of monomer x following monomer y
def draw_heatmap_new(monoseq_list):
    # define the 26 monomers (A to L, Z; A_rc to L_rc, Z_rc)
    monomers = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'Z',
                'A_rc', 'B_rc', 'C_rc', 'D_rc', 'E_rc', 'F_rc', 'G_rc', 'H_rc', 'I_rc', 'J_rc', 'K_rc', 'L_rc', 'Z_rc']

    # create a dictionary to store the counts of each pair of monomers
    monomer_dict = {}
    for monoseq in monoseq_list:
        for i in range(len(monoseq) - 1):
            pair = monoseq[i] + '+' + monoseq[i + 1]
            if pair not in monomer_dict:
                monomer_dict[pair] = 1
            else:
                monomer_dict[pair] += 1

    # create a matrix to store the counts
    count_matrix = np.zeros((26, 26), dtype=int)

    # fill in the matrix
    for pair, count in monomer_dict.items():
        monomer_x, monomer_y = pair.split('+')
        if monomer_x in monomers and monomer_y in monomers:
            x = monomers.index(monomer_x)
            y = monomers.index(monomer_y)
            count_matrix[x, y] = count

    # plot the heatmap
    plt.figure(figsize=(13, 10))
    plt.imshow(count_matrix, cmap='viridis', interpolation='nearest')
    plt.colorbar(label='Count')
    # adjust ticks so that labels are centered in each cell
    plt.xticks(ticks=np.arange(26), labels=monomers, rotation=45, ha='right')
    plt.yticks(ticks=np.arange(26), labels=monomers)
    # draw gridlines to align with cells
    plt.gca().set_xticks(np.arange(-0.5, 26, 1), minor=True)
    plt.gca().set_yticks(np.arange(-0.5, 26, 1), minor=True)
    plt.grid(which='minor', color='white', linestyle='-', linewidth=0.5)
    plt.tick_params(which='minor', size=0)  # hide minor ticks
    # annotate each cell with its count
    for i in range(26):
        for j in range(26):
            count = count_matrix[i, j]
            if count > 0:  # only annotate cells with non-zero counts
                color = 'black' if count > 1000 else 'white'
                plt.text(j, i, str(count), ha='center', va='center', color=color, fontsize=8)
    # set labels and title
    plt.xlabel('Monomer X')
    plt.ylabel('Monomer Y')
    plt.title('Count Heatmap of Monomer Sequences (x follows y)')
    plt.show()


def main():
    # extract all the file paths in the alignment_results folder
    file_list = []
    for file in os.listdir('./results/100_original_monomers'):
        if file.endswith('.tsv'):
            file_list.append('./results/100_original_monomers/' + file)

    # process the data
    monoseq_list = process_all_to_lists(file_list)
    # draw the heatmap
    draw_heatmap_original(monoseq_list)
    # draw_heatmap_new(monoseq_list)

main()