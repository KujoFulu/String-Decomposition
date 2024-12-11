import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# # this function draw the violin plot for each monomer's percent identity across one certain read
# # input file path is a string
# def single_read_violin_plot(file):
#     # read in file
#     data = pd.read_csv(file, sep='\t')

#     # extract monomer and match percentage
#     violin_data = data[['Monomer', 'Match_Percentage']]
#     violin_data['Monomer'] = violin_data['Monomer'].str.extract(r'([A-Z])')

#     # generate a unique color for each monomer
#     monomers = sorted(violin_data['Monomer'].unique())
#     palette = sns.color_palette("husl", len(monomers))

#     # draw violin plot
#     plt.figure(figsize=(12, 6))
#     sns.violinplot(data=violin_data, x='Monomer', y='Match_Percentage', order=sorted(violin_data['Monomer'].unique()), palette=dict(zip(monomers, palette)))
#     plt.title("Violin Plot of Percent Identity Across 12 Monomers")
#     plt.xlabel("Monomers")
#     plt.ylabel("Percent Identity (%)")
#     plt.ylim(40, 105)
#     plt.xticks(rotation=45)
#     plt.tight_layout()
#     plt.show()

# # # test
# # file_name = './alignment_results/chrX-v0_3324_aligned_1890_R_3_23654_91_alignments.tsv'
# # single_read_violin_plot(file_name)

# this function draw the violin plot for each monomer's percent identity across all reads
# input file path is a list of file paths
def combined_violin_plot(file_list):
    # combine data from all files
    combined_data = pd.DataFrame()
    for file in file_list:
        data = pd.read_csv(file, sep='\t')
        combined_data = pd.concat([combined_data, data], ignore_index=True)

    # extract monomer and match percentage
    violin_data = combined_data[['Monomer', 'Match_Percentage']]
    violin_data['Monomer'] = violin_data['Monomer'].str.extract(r'([A-Z])')
    print(violin_data)

    # generate a unique color for each monomer
    monomers = sorted(violin_data['Monomer'].unique())
    palette = sns.color_palette("husl", len(monomers))

    # draw violin plot
    plt.figure(figsize=(12, 6))
    sns.violinplot(data=violin_data, x='Monomer', y='Match_Percentage', order=sorted(violin_data['Monomer'].unique()), palette=dict(zip(monomers, palette)))

    # # draw low identity points as outliers
    # # create a df to extract all lines with identity < 50
    # low_identity = violin_data[violin_data['Match_Percentage'] < 50]

    # # check Z
    # all_Z = violin_data[violin_data['Monomer'] == 'Z']
    # ave = all_Z['Match_Percentage'].mean() # average Z match percentage
    # low_Z = all_Z[all_Z['Match_Percentage'] < 60]
    # print(len(all_Z), len(low_Z), len(ave))

    # # Overlay scatter plot
    # sns.stripplot(
    #     data=low_identity,
    #     x='Monomer',
    #     y='Match_Percentage',
    #     order=monomers,
    #     color='black',  # Scatter point color
    #     size=4,  # Scatter point size
    #     alpha=0.6,  # Transparency for better visibility
    #     jitter=True  # Add jitter to avoid overlap
    # )

    # plot info
    plt.title("Violin Plot of Match Percentage Across 12 Monomers")
    plt.xlabel("Monomers")
    plt.ylabel("Match Percentage (%)")
    plt.ylim(0, 105)
    plt.axhline(y=100, color='red', linestyle='--', linewidth=1.5, label="100% Threshold")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

def main():
    # extract all the file paths in the alignment_results folder
    file_list = []
    for file in os.listdir('./results/100_original_monomers'):
        if file.endswith('.tsv'):
            file_list.append('./results/100_original_monomers/' + file)

    # plot
    combined_violin_plot(file_list)

main()