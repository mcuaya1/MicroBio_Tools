import argparse
from datetime import datetime
import math
import pandas as pd
import matplotlib.pyplot as plt
import os
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from scipy import stats
from qiime2 import Metadata
from qiime2 import Artifact
import numpy as np

def validate_data(asv_table):
    if '.qza' in asv_table:
        asv_table = Artifact.load(asv_table)
        return asv_table

    return None

def stats_generator(dataframe,
                    output):
    # Generate excel file filed with distance points/sig test
    print('Generating excel file...')
    dataframe.to_excel(f'{output}correlation_analysis_stats.xlsx')

    # Generate markdown file
    print('Generating markdown file with table stats...')
    time_generated=datetime.now().strftime("%d/%m/%y %H:%M:%S")
    with open(f'{output}correlation_analysis_stats.md', "w") as f:
        f.write(f'''#Correlation analysis stats\n
                ## To find further sequence specific information, refer to table 03 generated previously\n
                **Please refer to the excel or csv file generated to perform further analysis.**\n
                Date file was generated: {time_generated}\n
                ## Spearman results
                {dataframe.to_markdown()}''')

    # Generate html file
    print('Generating html file with table stats...')
    with open(f'{output}correlation_analysis_stats.html', "w") as f:
        f.write(f'''<!doctype html>
    <html lang="en">
        <head>
            <meta charset="utf-8">
            <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css">
        </head>
        <body>
            <h1>Correlation analysis stats</h1>
            <h2>To find further sequence specific information, refer to table 03 generated previously.</h2>
            <strong>Please refer to the excel file generated to perform further analysis. </strong>
            <p>Date file was generated: {time_generated}</p>
            <h2>Spearman results</h2>
            {dataframe.to_html()}''')



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

def visualizers(dataframe,
                plot_title,
                xcell_spacing=3.5,
                ycell_spacing=1.6,
                fontsize=14,
                linewidth=0):
    dataframe = dataframe.reindex(index=dataframe.index[::-1])
    stats = dataframe.to_numpy()
    fig, ax = plt.subplots(figsize=(12,12))


    patches = []
    colors = []
    for y in range(stats.shape[0]):
        for x in range(stats.shape[1]):
            circle = Circle((x*xcell_spacing, y*ycell_spacing),
                            radius=(stats[y][x]))
            colors.append(stats[y][x])
            patches.append(circle)
    

    # Add all circles as a single collection
    collection = PatchCollection(patches,
                                 cmap='coolwarm',
                                 edgecolor='black',
                                 linestyle='solid',
                                 linewidth=linewidth)
    collection.set_array(np.array(colors))  # Signed values for coloring
    collection.set_clim(vmin=-1, vmax=1)

    corr_labels = dataframe.columns.to_list()
    taxa = dataframe.index.to_list()
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)



    ax.add_collection(collection)
    cbar = fig.colorbar(collection,
                        ax=ax)
    cbar.ax.set_title(label='Legend',
                      pad=25.0,
                      fontsize=20)
    cbar.ax.tick_params(labelsize=fontsize)
    
    cbar.outline.set_visible(False)

    ax.set_aspect(aspect='equal',
                  adjustable='box',
                  anchor='SW',
                  share=False) 

    ax.set_xlim(xmin=0,
                xmax=stats.shape[1]*xcell_spacing,
                auto=False)
    ax.set_ylim(ymin=0,
                ymax=stats.shape[0]*ycell_spacing,
                auto=False)
    ax.set_xticks([x * xcell_spacing for x in range(stats.shape[1])])
    ax.set_yticks([y * ycell_spacing for y in range(stats.shape[0])])

    ax.autoscale(enable=True)
    ax.set_xticklabels(corr_labels,
                       fontsize=fontsize,
                       rotation=45,
                       ha='right')
    ax.set_yticklabels(taxa,
                       fontsize=fontsize)
    ax.yaxis.tick_right()


    

    plt.title(f'{plot_title}',
              fontsize=20,
              pad=25.0)
    fig.tight_layout()
    fig.savefig(f"{output}corr_analysis.png", dpi=300)


def correlation_analysis_scatter(data_file,
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
    print('Mapping sample IDs to treatments...')
    sample_mapping = map_file.get_column(corr_col_0).drop_missing_values().to_dataframe()

    sample_mapping = sample_mapping.reset_index()

    sample_mapping = sample_mapping.set_index(f'{corr_col_0}')


    sample_mapping = sample_mapping[sample_mapping.index.isin(treatments)]


    print(sample_mapping)


    # Map Correlation values to IDs
    print('Mapping correlation samples to Sample IDs...')
    sample_ids = sample_mapping['#SampleID'].to_list()
    correlation_mapping = map_file.to_dataframe()[correlation_cols]
    correlation_mapping = correlation_mapping[correlation_mapping.index.isin(sample_ids)]
    

    # Merging mappings
    print('Merging both mapping tables....')
    merged_mapping = pd.merge(sample_mapping,
                      correlation_mapping,
                      right_index=True,
                      left_on=["#SampleID"])

    print(merged_mapping)
    asv_table = data_file.view(pd.DataFrame)
    
    
    # Extract top n asvs from csv file
    #print('Extracting top n ASVs from ASV table...')
    #top_n = pd.read_excel(f"{taxa_file}")
    #top_n = top_n.rename(columns={'Unnamed: 0': 'Treatments'})
    #top_n = top_n.set_index('Treatments')
    #top_n_list = top_n.columns.to_list()

    #asv_table = asv_table[top_n_list[:len(top_n_list)-1]]
    #asv_table = asv_table[asv_table.index.isin(sample_ids)]
    print(asv_table) 
    print('Merging previous mapping table with ASV mapping table')
    merged_mapping = pd.merge(merged_mapping,
                              asv_table,
                              right_index=True,
                              left_on=["#SampleID"])

    # Drop SampleIDs from DataFrame and merge treatments by mean of correlation values
    print('Merging replicates from correlation table by getting the MEAN across all replicates...')
    merged_mapping = merged_mapping.drop(columns=['#SampleID'])
    corr_dataframe = merged_mapping[correlation_cols].astype(float)
    corr_dataframe = corr_dataframe.groupby(level=0).mean()
    print(corr_dataframe)

    # Drop correlation columns from DataFrame and merge treatments by sum of taxanomic abundance
    print('Merging replicates from ASV table and find the mean between counts...')
    taxa_dataframe = merged_mapping.drop(columns=correlation_cols)
    taxa_dataframe = taxa_dataframe.astype(float).groupby(level=0).mean()
    
#    # Normalize taxa DataFrame
#    treatment_total = taxa_dataframe[taxa_dataframe.columns].sum(axis=0)
#    print(treatment_total)
#    taxa_dataframe = taxa_dataframe.div(treatment_total,
#                                        axis=1)
    #new_lables = taxa_dataframe.columns.to_list()
    
    #asv_label_formatter(new_lables)
    #taxa_dataframe.columns = new_lables
    
    # Merge two DataFrames back together
    print('Merging both tables back together...')
    merged_mapping = pd.merge(corr_dataframe,
                              taxa_dataframe,
                              right_index=True,
                              left_index=True)
    print('Calculating spearman correlation value...')
    results = {}
    results_formatted = {}
    print(taxa_dataframe)
    for corr_col in corr_dataframe.columns.to_list():
        curr_col = corr_col
        results[curr_col] = {} 
        results_formatted[curr_col] = {}
        for taxa in taxa_dataframe.columns.to_list():
            print(taxa_dataframe[taxa])
            print(taxa)
            print(type(taxa_dataframe[taxa]))
            spearman = stats.spearmanr(taxa_dataframe[taxa].to_list(),
                                                      corr_dataframe[curr_col].to_list())
            results[corr_col][taxa] = f'Statistic: {spearman[0]}, pvalue: {spearman[-1]}'
            results_formatted[corr_col][taxa] = spearman[0]
    final_dataframe = pd.DataFrame.from_dict(results,
                                  orient='index')
    results_formatted = pd.DataFrame.from_dict(results_formatted,
                                               orient='index')
    results_formatted = results_formatted.T
    final_dataframe = final_dataframe.T
    print(final_dataframe)
    #print(results_formatted)
    exit(1)
    visualizers(results_formatted,
                plot_title=plot_title)
    stats_generator(final_dataframe,
                    output)
    print('Done!')



def correlation_analysis_nomralized(data_file,
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
    print('Mapping sample IDs to treatments...')
    sample_mapping = map_file.get_column(corr_col_0).drop_missing_values().to_dataframe()

    sample_mapping = sample_mapping.reset_index()

    sample_mapping = sample_mapping.set_index(f'{corr_col_0}')


    sample_mapping = sample_mapping[sample_mapping.index.isin(treatments)]


    print(sample_mapping)


    # Map Correlation values to IDs
    print('Mapping correlation samples to Sample IDs...')
    sample_ids = sample_mapping['#SampleID'].to_list()
    correlation_mapping = map_file.to_dataframe()[correlation_cols]
    correlation_mapping = correlation_mapping[correlation_mapping.index.isin(sample_ids)]
    

    # Merging mappings
    print('Merging both mapping tables....')
    merged_mapping = pd.merge(sample_mapping,
                      correlation_mapping,
                      right_index=True,
                      left_on=["#SampleID"])

    print(merged_mapping)
    asv_table = data_file.view(pd.DataFrame)
    
    
    # Extract top n asvs from csv file
    print('Extracting top n ASVs from ASV table...')
    top_n = pd.read_excel(f"{taxa_file}")
    top_n = top_n.rename(columns={'Unnamed: 0': 'Treatments'})
    top_n = top_n.set_index('Treatments')
    top_n_list = top_n.columns.to_list()

    asv_table = asv_table[top_n_list[:len(top_n_list)-1]]
    asv_table = asv_table[asv_table.index.isin(sample_ids)]
    print(asv_table) 
    print('Merging previous mapping table with ASV mapping table')
    merged_mapping = pd.merge(merged_mapping,
                              asv_table,
                              right_index=True,
                              left_on=["#SampleID"])

    # Drop SampleIDs from DataFrame and merge treatments by mean of correlation values
    print('Merging replicates from correlation table by getting the MEAN across all replicates...')
    merged_mapping = merged_mapping.drop(columns=['#SampleID'])
    corr_dataframe = merged_mapping[correlation_cols].astype(float)
    corr_dataframe = corr_dataframe.groupby(level=0).mean()
    print(corr_dataframe)

    # Drop correlation columns from DataFrame and merge treatments by sum of taxanomic abundance
    print('Merging replicates from ASV table and find the mean between counts...')
    taxa_dataframe = merged_mapping.drop(columns=correlation_cols)
    taxa_dataframe = taxa_dataframe.astype(float).groupby(level=0).mean()
    
#    # Normalize taxa DataFrame
#    treatment_total = taxa_dataframe[taxa_dataframe.columns].sum(axis=0)
#    print(treatment_total)
#    taxa_dataframe = taxa_dataframe.div(treatment_total,
#                                        axis=1)
    new_lables = taxa_dataframe.columns.to_list()
    
    asv_label_formatter(new_lables)
    taxa_dataframe.columns = new_lables
    print(taxa_dataframe)
    
    # Merge two DataFrames back together
    print('Merging both tables back together...')
    merged_mapping = pd.merge(corr_dataframe,
                              taxa_dataframe,
                              right_index=True,
                              left_index=True)
    print('Calculating spearman correlation value...')
    results = {}
    results_formatted = {}
    for corr_col in corr_dataframe.columns.to_list():
        curr_col = corr_col
        results[curr_col] = {} 
        results_formatted[curr_col] = {}
        for taxa in taxa_dataframe.columns.to_list():
            spearman = stats.spearmanr(taxa_dataframe[taxa].to_list(),
                                                      corr_dataframe[curr_col].to_list())
            results[corr_col][taxa] = f'Statistic: {round(spearman[0], 6)}, pvalue: {round(spearman[-1], 3)}'
            results_formatted[corr_col][taxa] = round(spearman[0], 6)
    final_dataframe = pd.DataFrame.from_dict(results,
                                  orient='index')
    results_formatted = pd.DataFrame.from_dict(results_formatted,
                                               orient='index')
    results_formatted = results_formatted.T
    final_dataframe = final_dataframe.T
    print(final_dataframe)
    print(results_formatted)
    visualizers(results_formatted,
                plot_title=plot_title)
    stats_generator(final_dataframe,
                    output)
    print('Done!')


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
    print('Mapping sample IDs to treatments...')
    sample_mapping = map_file.get_column(corr_col_0).drop_missing_values().to_dataframe()

    sample_mapping = sample_mapping.reset_index()

    sample_mapping = sample_mapping.set_index(f'{corr_col_0}')


    sample_mapping = sample_mapping[sample_mapping.index.isin(treatments)]


    print(sample_mapping)


    # Map Correlation values to IDs
    print('Mapping correlation samples to Sample IDs...')
    sample_ids = sample_mapping['#SampleID'].to_list()
    correlation_mapping = map_file.to_dataframe()[correlation_cols]
    correlation_mapping = correlation_mapping[correlation_mapping.index.isin(sample_ids)]
    print(correlation_mapping)
    

    # Merging mappings
    print('Merging both mapping tables....')
    merged_mapping = pd.merge(sample_mapping,
                      correlation_mapping,
                      right_index=True,
                      left_on=["#SampleID"])

    print(merged_mapping)
    asv_table = data_file.view(pd.DataFrame)
    
    # Extract top n asvs from csv file
    print('Extracting top n ASVs from ASV table...')
    top_n = pd.read_excel(f"{taxa_file}")
    top_n = top_n.rename(columns={'Unnamed: 0': 'Treatments'})
    top_n = top_n.set_index('Treatments')
    top_n_list = top_n.columns.to_list()

    asv_table = asv_table[top_n_list[:len(top_n_list)-1]]
    asv_table = asv_table[asv_table.index.isin(sample_ids)]


    print('Merging previous mapping table with ASV mapping table')
    merged_mapping = pd.merge(merged_mapping,
                              asv_table,
                              right_index=True,
                              left_on=["#SampleID"])
    print(merged_mapping)

    # Drop SampleIDs from DataFrame and merge treatments by mean of correlation values
    print('Merging replicates from correlation table by getting the MEAN across all replicates...')
    merged_mapping = merged_mapping.drop(columns=['#SampleID'])
    corr_dataframe = merged_mapping[correlation_cols].astype(float)
    corr_dataframe = corr_dataframe.groupby(level=0).mean()
    print(corr_dataframe)

    # Drop correlation columns from DataFrame and merge treatments by sum of taxanomic abundance
    print('Merging replicates from ASV table and normalizing raw count values...')
    taxa_dataframe = merged_mapping.drop(columns=correlation_cols)
    taxa_dataframe = taxa_dataframe.astype(int).groupby(level=0).sum()
    print(taxa_dataframe)
    
    # Normalize taxa DataFrame
    treatment_total = taxa_dataframe[taxa_dataframe.columns].sum(axis=0)
    print(treatment_total)
    taxa_dataframe = taxa_dataframe.div(treatment_total,
                                        axis=1)

    new_lables = taxa_dataframe.columns.to_list()
    
    asv_label_formatter(new_lables)
    taxa_dataframe.columns = new_lables
    print(taxa_dataframe)
    
    # Merge two DataFrames back together
    print('Merging both tables back together...')
    merged_mapping = pd.merge(corr_dataframe,
                              taxa_dataframe,
                              right_index=True,
                              left_index=True)

    print('Calculating spearman correlation value...')
    results = {}
    results_formatted = {}
    for corr_col in corr_dataframe.columns.to_list():
        curr_col = corr_col
        results[curr_col] = {} 
        results_formatted[curr_col] = {}
        for taxa in taxa_dataframe.columns.to_list():
            spearman = stats.spearmanr(taxa_dataframe[taxa].to_list(),
                                                      corr_dataframe[curr_col].to_list())
            results[corr_col][taxa] = f'Statistic: {round(spearman[0], 6)}, pvalue: {round(spearman[-1], 3)}'
            results_formatted[corr_col][taxa] = round(spearman[0], 6)
    final_dataframe = pd.DataFrame.from_dict(results,
                                  orient='index')
    results_formatted = pd.DataFrame.from_dict(results_formatted,
                                               orient='index')
    results_formatted = results_formatted.T
    final_dataframe = final_dataframe.T
    print(final_dataframe)
    print(results_formatted)
    visualizers(results_formatted,
                plot_title=plot_title)
    stats_generator(final_dataframe,
                    output)
    print('Done!')


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
        correlation_analysis_nomralized(asv_table,
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
