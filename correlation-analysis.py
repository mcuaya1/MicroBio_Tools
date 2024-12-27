import argparse
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
import os
import scipy.stats as stats

from qiime2 import Metadata
from qiime2 import Artifact


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
            asv_list[i] = asv_list[i].split(';')[-1]

def rel_taxa_abundance(asv_table, map_file, data_column, treatments, num):
    treatments = treatments[0].split(',')
    print('Treatments to be processed...')
    for i in range(len(treatments)):
        print(treatments[i], end='\t')
    n = len(treatments)
    dataframe_list = []
    all_samples = asv_table.columns.to_list()
    for i in range(n):
        current_treatment = treatments[i]

        samples = list(map_file.get_ids(f"[{data_column}]='{current_treatment}'"))
        filtered_samples = []
        for sample in samples:
            if sample not in all_samples:
                print(f"{sample} is not in the ASV table, please check raw counts file for this sequence run")
            else:
                filtered_samples.append(sample)
        samples = filtered_samples
        temp_df = asv_table[samples]
        temp_df[f'{current_treatment}'] = temp_df[temp_df.columns].sum(axis=1)
        list_temp = temp_df.columns
        list_temp = list_temp[0:len(list_temp)-1]
        temp_df = temp_df.drop(columns=list_temp)
        dataframe_list.append(temp_df)

    merged_data = pd.concat(dataframe_list, axis=1)

    print("Merged, grouped, and filtered down table...")
    print(merged_data)

    counter = 0
    curr_row = 0
    treatments = merged_data.columns.to_list()
    top_n_taxa = []

    print(f"Finding top {num} ASVs...")
    while (counter < num):
        for i in range(len(treatments)):
            sorted_table = merged_data.sort_values(by=f'{treatments[i]}', ascending=False)
            taxa = sorted_table.iloc[curr_row].name
            if taxa not in top_n_taxa:
                top_n_taxa.append(taxa)
                counter += 1
            if counter >= num:
                break
        curr_row += 1

    other_df = merged_data.drop(top_n_taxa, axis=0)
    top_taxa_df = merged_data.drop(other_df.index, axis=0)
    other_total = other_df[other_df.columns].sum(axis=0)
    other_total = other_total.to_frame().T.rename(index={0: 'Other'})

    top_taxa_df = top_taxa_df.sort_values(by=top_taxa_df.columns.to_list()[0], ascending=False)
    top_n_taxa = top_taxa_df.index.to_list()
    raw_asv_strings = top_taxa_df.index.to_list()
    raw_asv_strings.append("Other")
    asv_label_formatter(top_n_taxa)

    top_n_taxa.append("Other")
    top_taxa_df = pd.concat([top_taxa_df, other_total], axis=0, ignore_index=True)
    top_taxa_df.index = top_n_taxa
    print(f"Found top {num} ASVs...")
    return top_taxa_df


def correlation_analysis(asv_table, map_file, data_column, treatments, against_column, plot_tilte, num, output):
    asv_table = asv_table.view(pd.DataFrame)
    asv_table = asv_table.T
    rel_abundance_table = rel_taxa_abundance(asv_table, map_file, data_column, treatments, num)
    print(rel_abundance_table)



def validate_data(asv_table) -> None:
    if '.qza' in asv_table:
        asv_table = Artifact.load(asv_table)
        return asv_table

    return None


if __name__ == '__main__':
    pd.options.mode.chained_assignment = None
    parser = argparse.ArgumentParser(add_help=False, prog="correlation-analysis.py", description="Program to generate custom correlation analysis plots")

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
    parser.add_argument('-a',
                        "--against",
                        required=True,
                        help="Colum to corrlate against",
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
    parser.add_argument('-n',
                        "--top-n-taxa",
                        required=True,
                        help="Filter for top N taxa",
                        type=int)
    args = parser.parse_args()

    data_file = args.input_file
    map_file = args.map_file
    data_column = args.column
    plot_tilte = args.plot_title
    treatments = args.listing
    against_column = args.against
    n_taxa = args.top_n_taxa
    output = os.path.join(args.output_dir, "correlation-output/")


    if ((asv_table := validate_data(data_file))) and ((map_file := Metadata.load(map_file))):
        if not os.path.exists(output):
            os.mkdir(output)
        correlation_analysis(asv_table,
                             map_file,
                             data_column,
                             treatments,
                             against_column,
                             plot_tilte,
                             n_taxa,
                             output)
    else:
        print('Invalid data type or map file')
        exit(1)
