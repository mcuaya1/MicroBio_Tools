import argparse
import datetime as dt
import pandas as pd
import os
from skbio import OrdinationResults
from qiime2.plugins import feature_table
from qiime2.plugins import diversity
from qiime2 import Metadata
from qiime2 import Artifact
import matplotlib.pyplot as plt
import numpy as np

#Stats output and sig_test(adonis test)

def beta_diversity(asv_table, map_file, data_column, treatments, plot_tilte, output):
    treatments = tuple(treatments[0].split(','))
    # Filter asv table to include only samples from specified group
    asv_table_filtered = feature_table.methods.filter_samples(table=asv_table,
                                                              metadata=map_file,
                                                              where=f"[{data_column}] IN {treatments}")
    asv_table_filtered = asv_table_filtered.filtered_table

    # Preform Braycurtis metric
    beta_results = diversity.pipelines.beta(
            table=asv_table_filtered,
            metric='braycurtis')

    beta_diversity_table = beta_results.distance_matrix
#    # https://forum.qiime2.org/t/load-distancematrix-artifact-to-dataframe/11660
#    # Turn table into a distance matrix
#
    pcoa_results = diversity.methods.pcoa(distance_matrix=beta_diversity_table)
    pcoa_results = pcoa_results.pcoa
    pcoa_results = pcoa_results.view(OrdinationResults)
    print(pcoa_results)
    eigen_values = pcoa_results.eigvals
    total_eigen_values = eigen_values.sum()
#
# Provides dataframe
# https://medium.com/@conniezhou678/applied-machine-learning-part-12-principal-coordinate-analysis-pcoa-in-python-5acc2a3afe2d
# https://www.tutorialspoint.com/numpy/numpy_matplotlib.htm
    pcoa_results = pcoa_results.samples
#    pcoa_results.columns = pcoa_results.index.to_list()
    fig, ax = plt.subplots(figsize=(15, 10))
    cmap = plt.get_cmap('tab20')
    colors = [cmap(i) for i in np.linspace(0, 1, len(treatments))]
    mapping = {}
    for i in range(len(colors)):
        mapping[treatments[i]] = colors[i]
    column = map_file.get_column(data_column)
    for row in pcoa_results.itertuples():
        # print(row.Index, row[1], row[2])
        label = column.get_value(row.Index)
        ax.scatter(row[1], row[2], color=mapping[label], label=label)
    plt.ylabel(f'Axis 2 [{(eigen_values[1]/total_eigen_values):.2%}]', fontsize='15') 
    plt.xlabel(f'Axis 1 [{(eigen_values[0]/total_eigen_values):.2%}]', fontsize='15') 

# https://stackoverflow.com/questions/13588920/stop-matplotlib-repeating-labels-in-legend
    handles, labels = plt.gca().get_legend_handles_labels()
    uniques = dict(zip(labels, handles))
    ax.legend(uniques.values(), uniques.keys(), bbox_to_anchor=(1, 1), frameon=False, title="Treatments", alignment='left')
    plt.title(f'{plot_tilte}', fontsize='20') 
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.grid(True)
    plt.tight_layout()
    fig.savefig(f"{output}beta_diversity.png")

#Further resources can be found at the following links below:
#https://develop.qiime2.org/en/latest/intro.html
#https://docs.qiime2.org/2024.5/plugins/


def validate_data(asv_table) -> None:
    # Check if data is a qza type
    if '.qza' in asv_table:
        asv_table = Artifact.load(asv_table)
        return asv_table

    return None


if __name__ == '__main__':
    parser = argparse.ArgumentParser(add_help=False, prog="betsa-diversity-genator.py", description="Program to generate custom beta diversity scatter plots")
    parser.add_argument('-i',"--input-file", required=True, help="Imported feature table",type=str)
    parser.add_argument('-m',"--map-file", required=True, help="Map file for data",type=str)
    parser.add_argument('-c',"--column", required=True, help="Colmun to parse for data",type=str)
    parser.add_argument('-p', "--plot-title", help="Tilte for plot",type=str)
    parser.add_argument('-l', "--listing", nargs='+', type=str, help="Set a preferred listing for x axis (Default is nothing)")
    parser.add_argument('-d', "--output-dir", required=True, help="Output directory location",type=str)
    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Display commands possible with this program.')
    args = parser.parse_args()
    data_file = args.input_file
    map_file = args.map_file
    data_column = args.column
    plot_tilte = args.plot_title
    treatments = args.listing
    output = os.path.join(args.output_dir, "beta-diversity/")
    # output=args.output_dir
    # Load in ASV table and map file
    if ((asv_table := validate_data(data_file)) != None) and ((map_file := Metadata.load(map_file)) != None):
        if not os.path.exists(output):
            os.mkdir(output)
        beta_diversity(asv_table, map_file, data_column, treatments, plot_tilte, output)
    else:
        print('Invalid data type or map file')
        exit(1)
