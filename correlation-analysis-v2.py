import argparse
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
import os

from scipy import stats
from qiime2 import Metadata
from qiime2 import Artifact

def validate_data(asv_table):
    if '.qza' in asv_table:
        asv_table = Artifact.load(asv_table)
        return asv_table

    return None


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


def correlation_analysis(data_file,
                         map_file,
                         corr_col_0,
                         corr_col_1,
                         treatments,
                         plot_title,
                         taxa_file,
                         output):

    treatments = treatments.split(',')
    correlation_cols = corr_col_1.split(',') 

    # Map sample IDs to treatments specified
    sample_mapping = map_file.get_column(corr_col_0).drop_missing_values().to_dataframe()

    sample_mapping = sample_mapping.reset_index()

    sample_mapping = sample_mapping.set_index(f'{corr_col_0}')


    sample_mapping = sample_mapping[sample_mapping.index.isin(treatments)]


    # Map Correlation values to IDs
    sample_ids = sample_mapping['#SampleID'].to_list()
    correlation_mapping = map_file.to_dataframe()[correlation_cols]
    correlation_mapping = correlation_mapping[correlation_mapping.index.isin(sample_ids)]
    

    # Merging mappings
    merged_mapping = pd.merge(sample_mapping,
                      correlation_mapping,
                      right_index=True,
                      left_on=["#SampleID"])

    print(merged_mapping)
    asv_table = data_file.view(pd.DataFrame)
    
    # Extract top n asvs from csv file
    top_n = pd.read_excel(f"{taxa_file}")
    top_n = top_n.rename(columns={'Unnamed: 0': 'Treatments'})
    top_n = top_n.set_index('Treatments')
    top_n_list = top_n.columns.to_list()

    asv_table = asv_table[top_n_list[:len(top_n_list)-1]]
    asv_table = asv_table[asv_table.index.isin(sample_ids)]


    merged_mapping = pd.merge(merged_mapping,
                              asv_table,
                              right_index=True,
                              left_on=["#SampleID"])

    # Drop SampleIDs from DataFrame and merge treatments by mean of correlation values
    merged_mapping = merged_mapping.drop(columns=['#SampleID'])
    corr_dataframe = merged_mapping[correlation_cols].astype(float)
    corr_dataframe = corr_dataframe.groupby(level=0).mean()


    # Drop correlation columns from DataFrame and merge treatments by sum of taxanomic abundance
    taxa_dataframe = merged_mapping.drop(columns=correlation_cols)
    taxa_dataframe = taxa_dataframe.astype(int).groupby(level=0).sum()

    # Normalize taxa DataFrame
    treatment_total = taxa_dataframe[taxa_dataframe.columns].sum(axis=1)
    taxa_dataframe = taxa_dataframe.div(treatment_total,
                                        axis=0)

    new_lables = taxa_dataframe.columns.to_list()
    
    asv_label_formatter(new_lables)
    taxa_dataframe.columns = new_lables
    
    # Merge two DataFrames back together
    merged_mapping = pd.merge(corr_dataframe,
                              taxa_dataframe,
                              right_index=True,
                              left_index=True)

    results = {}
    for corr_col in corr_dataframe.columns.to_list():
        curr_col = corr_col
        results[corr_col] = {} 
        for taxa in taxa_dataframe.columns.to_list():
            spearman = stats.spearmanr(taxa_dataframe[taxa].to_list(),
                                                      corr_dataframe[curr_col].to_list())
            results[corr_col][taxa] = f'Statistic: {round(spearman[0], 6)}, pvalue: {round(spearman[-1], 3)}' 
    final_dataframe = pd.DataFrame.from_dict(results,
                                  orient='index')
    final_dataframe = final_dataframe.T
    print(final_dataframe)


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
    parser.add_argument('-t',
                        "--taxa-file",
                        required=True,
                        help="Top N taxa file location",
                        type=str)
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
    corr_col_0 = args.correlation_column_0
    corr_col_1 = args.correlation_column_1
    plot_title = args.plot_title
    taxa_file = args.taxa_file
    samples = args.samples
    output = os.path.join(args.output_dir, "correlation-output/")

    if ((asv_table := validate_data(data_file))) and ((map_file := Metadata.load(map_file))):
        if not os.path.exists(output):
            os.mkdir(output)
        correlation_analysis(asv_table,
                             map_file,
                             corr_col_0,
                             corr_col_1,
                             samples,
                             plot_title,
                             taxa_file,
                             output)
    else:
        print('Invalid data type or map file')
        exit(1)
