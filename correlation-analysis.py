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
                         corr_col,
                         plot_title,
                         top_tax_file,
                         output_dir) -> None:

    # Extract top n asvs from csv file
    top_n = pd.read_excel(f"{top_tax_file}")
    top_n = top_n.rename(columns={'Unnamed: 0': 'Treatments'})
    top_n = top_n.set_index('Treatments')
    new_lables = top_n.columns.to_list()
    asv_label_formatter(new_lables)

    # Extract correlation values average them


def validate_data(asv_table) -> None:
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
    parser.add_argument('-c',
                        "--correlation-column",
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
    corr_col = args.correlation_column
    plot_title = args.plot_title
    taxa_file = args.taxa_file
    output = os.path.join(args.output_dir, "correlation-output/")

    if ((asv_table := validate_data(data_file))) and ((map_file := Metadata.load(map_file))):
        if not os.path.exists(output):
            os.mkdir(output)
        correlation_analysis(map_file,
                             corr_col,
                             plot_title,
                             taxa_file,
                             output)
    else:
        print('Invalid data type or map file')
        exit(1)
