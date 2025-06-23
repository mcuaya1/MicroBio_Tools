import argparse
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
import os
import scipy.stats as stats

from qiime2 import Metadata
from qiime2 import Artifact

#To clean up asv labels
def asv_label_formatter(asv_list):
    for i in range(len(asv_list)):
        if 'Other' in asv_list[i] or '__':
            if 'g__' in asv_list[i]:
                asv_list[i] = asv_list[i].split(';')[5]
                
            elif 'f__' in asv_list[i]:
                asv_list[i] = asv_list[i].split(';')[4]
                
            elif 'o__' in asv_list[i]:
                asv_list[i] = asv_list[i].split(';')[3]
                
            elif 'c__' in asv_list[i]:
                asv_list[i] = asv_list[i].split(';')[2]
                
            elif 'p__' in asv_list[i]:
                asv_list[i] = asv_list[i].split(';')[1]
            else:
                asv_list[i] = asv_list[i].split(';')[0]
        else:
            asv_list[i]=asv_list[i].split(';')[-1]

def correlation_analysis(map_file,
                         corr_col_0,
                         corr_col_1,
                         samples_ids,
                         plot_title,
                         top_tax_file,
                         output_dir) -> None:

    # Extract top n asvs from csv file
    top_n = pd.read_excel(f"{top_tax_file}")
    top_n = top_n.rename(columns={'Unnamed: 0': 'Treatments'})
    top_n = top_n.set_index('Treatments')
    new_lables = top_n.columns.to_list()
    
    asv_label_formatter(new_lables)
    top_n.columns = new_lables


    # Extract all correlation columns into a list
    corr_cols = corr_col_1.split(",")
    
    # Extract corrleation samples
    corr_map = map_file.get_column(f"{corr_col_0}").drop_missing_values().to_dataframe()
    print(corr_map)
    sample_names = corr_map.index.to_list()
    # Merge all corrleation values with correlation samples
    for i, corr_col in enumerate(corr_cols):
        sample_vals = map_file.get_column(f"{corr_col}").to_dataframe().fillna(0)
    
        corr_map = corr_map.merge(sample_vals,
                             left_index=True,
                             right_index=True)
        corr_map[f'{corr_col}'] = pd.to_numeric(corr_map[f"{corr_col}"],
                             errors='coerce')
    corr_map = corr_map.set_index(f'{corr_col_0}')

    # Remove _Exp tag from sample names
    sample_names = corr_map.index.to_list()
    for i in range(len(sample_names)):
        sample_names[i] = sample_names[i].split('_Exp')[0]

    corr_map.index = sample_names
    corr_map.index.name = 'Samples'
    # Group samples by name and get the mean
    corr_map = corr_map.groupby(corr_map.index).mean()
    samples_ids = samples_ids.split(',')
    corr_map = corr_map[corr_map.index.isin(samples_ids)]
    merged_map = pd.merge(top_n,
                          corr_map,
                          left_index=True,
                          right_index=True,
                          how='left')

    print(corr_map)
    for col in corr_cols:
        fig, ax = plt.subplots(figsize=(15, 10))
        for i in range(len(new_lables)):
            ax.scatter(merged_map[new_lables[i]],
                       merged_map[f"{col}"],
                       label=f"{new_lables[i]}")

        plt.ylabel(f"{col.replace('_',' ')}", fontsize='15')
        plt.xlabel(f'ASV Abundance', fontsize='15')

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        handles, labels = plt.gca().get_legend_handles_labels()
        ax.legend(handles,
                  labels,
                  bbox_to_anchor=(1, 1),
                  frameon=False,
                  title="ASV",
                  alignment='left')
        fig.tight_layout()
        plt.grid(True)
        fig.savefig(f"{output}{col}_corrleation_graph.png")


def validate_data(asv_table):
    if '.qza' in asv_table:
        asv_table = Artifact.load(asv_table)
        return asv_table

    return None


if __name__ == '__main__':
    pd.options.mode.chained_assignment = None
    parser = argparse.ArgumentParser(add_help=False,
                                     prog="correlation-analysis.py",
                                     description="Program to generate custom correlation analysis plots")

    parser.add_argument('-i',
                        "--input-file",
                        required=True,
                        help="Imported feature table",
                        type=str)
    parser.add_argument('-m',
                        "--map-file",
                        required=True,
                        help="Correlation map file",
                        type=str)
    parser.add_argument('-s',
                        "--samples",
                        required=True,
                        help="Samples, please seperate samples with ',' in between",
                        type=str)
    parser.add_argument('-c0',
                        "--correlation-column-0",
                        required=True,
                        help="Colum to corrlate against",
                        type=str)
    parser.add_argument('-c1',
                        "--correlation-column-1",
                        required=True,
                        help="Colum to corrlate against",
                        type=str)
    parser.add_argument('-p',
                        "--plot-title",
                        help="Tilte for plot",
                        type=str)
    parser.add_argument('-d',
                        "--output-dir",
                        required=True,
                        help="Output directory location",
                        type=str)
    parser.add_argument('-t',
                        "--taxa-file",
                        required=True,
                        help="Top N taxa file location",
                        type=str)
    parser.add_argument('-h',
                        '--help',
                        action='help',
                        default=argparse.SUPPRESS,
                        help='Display commands possible with this program.')
    args = parser.parse_args()

    data_file = args.input_file
    map_file = args.map_file
    corr_col_0 = args.correlation_column_0
    corr_col_1 = args.correlation_column_1
    plot_title = args.plot_title
    taxa_file = args.taxa_file
    samples = args.samples
    output = os.path.join(args.output_dir, "correlation-output/")

    if ((asv_table := validate_data(data_file))) and ((map_file := Metadata.load(map_file))):
        if not os.path.exists(output):
            os.mkdir(output)
        correlation_analysis(map_file,
                             corr_col_0,
                             corr_col_1,
                             samples,
                             plot_title,
                             taxa_file,
                             output)
    else:
        print('Invalid data type or map file')
        exit(1)
