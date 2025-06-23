import argparse
from datetime import datetime
import os
from skbio import OrdinationResults
from skbio import DistanceMatrix
from skbio.stats.distance import permanova
from qiime2.plugins import feature_table
from qiime2.plugins import diversity
from qiime2 import Metadata
from qiime2 import Artifact
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from collections import defaultdict
import re


def qiime2_significance_test(distance_matrix,
                            metadata,
                            data_column,
                            output) -> None:
    column = metadata.get_column(data_column)
    permanova = diversity.visualizers.beta_group_significance(
            distance_matrix=distance_matrix,
            method='permanova',
            metadata=column,
            pairwise=True,
            permutations=999)
    permanova.visualization.save(f'{output}/significance_test.qzv')


# Signifcance test (PERMANOVA)
def significance_test(distance_matrix,
                     dataframe,
                     output,
                     metadata,
                     treatments,
                     data_column) -> pd.DataFrame:

    # Create empty dictionary to store results
    results_df = defaultdict(dict)

    # Transform qiime2 Distance Matrix object to skbio
    # DistanceMatrix object
    distance_matrix = distance_matrix.view(DistanceMatrix)
    print('Generating significance test...')

    all_ids = distance_matrix.ids

    for i in range(len(treatments)):
        # Set ith treatment
        treatment_a = treatments[i]

        # Get all samples from map file that relate to treatment_a
        a_ids = map_file.get_ids(f"[{data_column}]='{treatment_a}'")
        for j in range(i+1, len(treatments)):
            # Set jth treatment
            treatment_b = treatments[j]

            # Get all samples from map file that relate to treatment_b
            b_ids = map_file.get_ids(f"[{data_column}]='{treatment_b}'")

            # Union both b_ids and a_ids to get all samples
            # to be compared against
            compare_ids = list(a_ids.union(b_ids))

            for id in compare_ids:
                if id not in all_ids:
                    compare_ids.remove(id)

            # Filter DistanceMatrix to contain all samples to be compared
            # against
            filtered_dm = distance_matrix.filter(compare_ids)

            # Create a DataFrame object from filltered DistanceMatrix object
            dm_df = filtered_dm.to_data_frame()
            dm_df.index.names = ['IDs']

            # Create a DataFrame which maps IDs to treatments
            map_ids = map_file.get_column(f"{data_column}")
            map_ids = map_ids.filter_ids(compare_ids)
            map_ids = map_ids.to_dataframe()
            map_ids.index.names = ['IDs']
            # Merge to DataFrames to finish mapping IDs to treatments
            main_df = pd.merge(dm_df,
                               map_ids,
                               left_index=True,
                               right_index=True)

            # Run PERMANOVA test against ith and jth treatments
            results = permanova(filtered_dm,
                                main_df,
                                f'{data_column}')
            results_df[f"{treatment_a}"][f"{treatment_b}"] = {
                    "Sample size": results.get("sample size"),
                    "Permutations": results.get("number of permutations"),
                    "pseudo-F": round(results.get("test statistic"), 6),
                    "p-value": round(results.get("p-value"), 3)
                    }

    return pd.DataFrame.from_dict(results_df, orient='index')


# Generate statsics
def stats_generator(stats,
                    output,
                    sig_results) -> None:
    # Extract distance levels
    dists = stats.columns.to_list()

    # Drop all other levels and relabel columns
    dists_pts = stats.drop(columns=dists[5:])
    renum = dists_pts.columns.to_list()
    renum = [x + 1 for x in renum]
    dists_pts.columns = renum

    # Generate excel file filed with distance points/sig test
    print('Generating excel file...')
    dists_pts.to_excel(f'{output}beta_diversity_stats.xlsx')
    sig_results.to_excel(f'{output}significance_test_results.xlsx')

    # Generate markdown file
    print('Generating markdown file with table stats...')
    time_generated=datetime.now().strftime("%d/%m/%y %H:%M:%S")
    with open(f'{output}beta_diversity_stats.md', "w") as f:
        f.write(f'''#Beta diversity stats\n
                ## To find further sequence specific information, refer to table 03 generated previously\n
                **Please refer to the excel or csv file generated to perform further analysis.**\n
                Date file was generated: {time_generated}\n
                ## Distance points
                {dists_pts.to_markdown()}\n
                ## PERMANOVA results
                {sig_results.to_markdown()}''')

    # Generate html file
    print('Generating html file with table stats...')
    with open(f'{output}beta_diversity_stats.html', "w") as f:
        f.write(f'''<!doctype html>
    <html lang="en">
        <head>
            <meta charset="utf-8">
            <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css">
        </head>
        <body>
            <h1>Beta diversity stats</h1>
            <h2>To find further sequence specific information, refer to table 03 generated previously.</h2>
            <strong>Please refer to the excel file generated to perform further analysis. </strong>
            <p>Date file was generated: {time_generated}</p>
            <h2>Distance points</h2>
            {dists_pts.to_html()}
            <h2>PERMANOVA results</h2>
            {sig_results.to_html()}''')


def beta_diversity(asv_table,
                   map_file,
                   data_column,
                   treatments,
                   plot_tilte,
                   output) -> None:

    # Split treatments into list
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

    # Convert qiime2 distance martix object into skbio DistanceMatrix
    # https://forum.qiime2.org/t/load-distancematrix-artifact-to-dataframe/11660
    pcoa_results = diversity.methods.pcoa(distance_matrix=beta_diversity_table)
    pcoa_results = pcoa_results.pcoa
    pcoa_results = pcoa_results.view(OrdinationResults)

    # Output Ordination results and calculate sum of eigen values
    print(pcoa_results)
    eigen_values = pcoa_results.eigvals
    total_eigen_values = eigen_values.sum()

    # Extract distance points from pcoa results
    # https://medium.com/@conniezhou678/applied-machine-learning-part-12-principal-coordinate-analysis-pcoa-in-python-5acc2a3afe2d
    # https://www.tutorialspoint.com/numpy/numpy_matplotlib.htm
    pcoa_results = pcoa_results.samples

    # Preform significance test
    sig_results = significance_test(beta_diversity_table,
                                   pcoa_results,
                                   output,
                                   map_file,
                                   treatments,
                                   data_column)

    # Generate statsics
    stats_generator(pcoa_results,
                    output,
                    sig_results)

    # For testing purposes
    # qiime2_significance_test(beta_diversity_table,
    #                         map_file,
    #                         data_column,
    #                         output)

    fig, ax = plt.subplots(figsize=(15, 10))

    # Generate and assign color mapping
    mapping = {}
    for i in range(len(treatments)):
        match = re.search(r'Tm(\d+)', treatments[i])

        if not match:
            raise ValueError(f"Treatment name invalid: {treatments[i]}")
    
        Tm = int(match.group(1))
        if Tm == 0:
            mapping[treatments[i]] = 'blue'
        elif Tm == 154:
            mapping[treatments[i]] = 'orange'

    column = map_file.get_column(data_column)

    # Generate Scatter plot
    for row in pcoa_results.itertuples():
        label = column.get_value(row.Index)

        markers = [".", "o", "^", "s", "p", "P", "*", "H", "X", "D"]
        marker = re.search(r'T(\d+)', label)
        if not marker:
            raise ValueError(f"Label name invalid: {label}")
        marker = int(marker.group(1))
        ax.scatter(
            row[1],
            row[2],
            color=mapping[label],
            label=label,
            marker=markers[marker % len(markers)],
            s=150,
            zorder=2,
        )

    # Calculate distance axis
    plt.ylabel(f'Axis 2 [{(eigen_values[1]/total_eigen_values):.2%}]', fontsize='15')
    plt.xlabel(f'Axis 1 [{(eigen_values[0]/total_eigen_values):.2%}]', fontsize='15')

    # Filter out duplicates from legend table
    # https://stackoverflow.com/questions/13588920/stop-matplotlib-repeating-labels-in-legend
    handles, labels = plt.gca().get_legend_handles_labels()
    uniques = dict(zip(labels, handles))
    # maintain order of treatments in legend
    uniques = {k: uniques[k] for k in treatments if k in uniques}

    ax.legend(uniques.values(),
              uniques.keys(),
              bbox_to_anchor=(1, 1),
              frameon=False,
              title="Treatments",
              fontsize='15',
              title_fontsize='20',
              loc='upper left')

    # Save plot
    plt.title(f'{plot_tilte}', fontsize='20')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.tight_layout()
    plt.grid(True)
    fig.savefig(f"{output}beta_diversity.png", dpi=300)


def validate_data(asv_table) -> None:
    # Check if data is a qza type
    if '.qza' in asv_table:
        asv_table = Artifact.load(asv_table)
        return asv_table


if __name__ == '__main__':
    parser = argparse.ArgumentParser(add_help=False,
                                     prog="betsa-diversity-genator.py",
                                     description="Program to generate custom beta diversity scatter plots")

    parser.add_argument('-i',
                        "--input-file",
                        required=True,
                        help="Imported feature table",
                        type=str)

    parser.add_argument('-m',
                        "--map-file",
                        required=True,
                        help="Map file for data",
                        type=str)

    parser.add_argument('-c',
                        "--column",
                        required=True,
                        help="Colmun to parse for data",
                        type=str)

    parser.add_argument('-p',
                        "--plot-title",
                        help="Tilte for plot",
                        type=str)

    parser.add_argument('-l',
                        "--listing",
                        nargs='+',
                        type=str,
                        help="Set a preferred listing for x axis (Default is nothing)")

    parser.add_argument('-d',
                        "--output-dir",
                        required=True,
                        help="Output directory location",
                        type=str)

    parser.add_argument('-h',
                        '--help',
                        action='help',
                        default=argparse.SUPPRESS,
                        help='Display commands possible with this program.')

    args = parser.parse_args()
    data_file = args.input_file
    map_file = args.map_file
    data_column = args.column
    plot_tilte = args.plot_title
    treatments = args.listing
    output = os.path.join(args.output_dir, "beta-diversity/")

    # Load in ASV table and map file
    if ((asv_table := validate_data(data_file)) is not None) and ((map_file := Metadata.load(map_file)) is not None):
        if not os.path.exists(output):
            os.mkdir(output)
        beta_diversity(asv_table, map_file, data_column, treatments, plot_tilte, output)
    else:
        print('Invalid data type or map file')
        exit(1)
