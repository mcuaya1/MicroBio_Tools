import argparse
from datetime import datetime
import pandas as pd
from qiime2.plugins import feature_table
from qiime2.plugins import diversity
from qiime2 import Metadata
from qiime2 import Artifact
import matplotlib.pyplot as plt
from qiime2.plugins.diversity.visualizers import alpha_group_significance
import os


def visualizer(data_dict, plot_title, outputdir):
    cmap = plt.get_cmap('tab20')
    fig, ax = plt.subplots(figsize = (15, 10))
    medianprops = dict(linestyle='-', linewidth=1.5, color='black')
    flierprops = dict(marker='o', markerfacecolor='blue', markersize=7,markeredgecolor='none')
    boxprops = dict(facecolor=cmap(0.09))

    plt.boxplot(data_dict.values(), labels=data_dict.keys(), 
                patch_artist=True, 
                boxprops=boxprops, medianprops=medianprops, flierprops=flierprops, 
                widths=0.7)

    plt.xticks(rotation=90,fontsize='13')
    plt.yticks(fontsize='13')

    #Old code to possiblly assign different colors to each boxplot
    #ax.boxplot(data_dict['T1Tm0_CAR'],labels=['T1Tm0_CAR'], patch_artist=True)
    #ax.boxplot(data_dict['T2Tm0_CAR'],labels=['T2Tm0_CAR'])

    #https://stackoverflow.com/questions/52273543/creating-multiple-boxplots-on-the-same-graph-from-a-dictionary
    #https://matplotlib.org/stable/gallery/statistics/boxplot.html#sphx-glr-gallery-statistics-boxplot-py
    #https://stackoverflow.com/questions/32443803/adjust-width-of-box-in-boxplot-in-python-matplotlib
    plt.ylabel('Shannon Diversity', fontsize='15') 
    plt.title(f'{plot_tilte}', fontsize='20') 
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    fig.savefig(f"{output}alpha_plot.png")

def stats_generator(stats, outputdir):
    #Get the amount of samples in a meta group
    length_of_samples = []
    for value in stats.items():
        length_of_samples.append(len(value))

    #Generating stats
    alpha_diversity_stats = pd.DataFrame({'Sample Tag':stats.keys(), 'Amount of samples': length_of_samples,'(Sample, Shannon Index)':stats_dict.values()})
    alpha_diversity_stats.index = alpha_diversity_stats['Sample Tag']
    alpha_diversity_stats.drop('Sample Tag', axis=1, inplace=True)

    #Saving them to an execel file
    print('Generating excel file...')
    alpha_diversity_stats.to_excel(f'{output}alpha_diversity_stats.xlsx')

    print('Generating markdown file with table stats...')
    time_generated=datetime.now().strftime("%d/%m/%y %H:%M:%S")
    with open(f'{output}alpha_diversity_stats.md', "w") as f:
        f.write(f'''#Alpha diversity stats\n
                ## To find further sequence specific information, refer to table 03 generated previously\n
                **Please refer to the excel or csv file generated to perform further analysis.**\n
                Date file was generated: {time_generated}\n
                {alpha_diversity_stats.to_markdown()}
                ''')
        
    print('Generating html file with table stats...')
    with open(f'{outputdir}alpha_diversity_stats.html', "w") as f:
        f.write(f'''<!doctype html>
    <html lang="en">
        <head>
            <meta charset="utf-8">
            <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css">
        </head>
        <body>
            <h1>Alpha diversity stats</h1>
            <h2 >To find further sequence specific information, refer to table 03 generated previously.</h2>
            <strong>Please refer to the excel file generated to perform further analysis. </strong>
            <p>Date file was generated: {time_generated}</p>
            {alpha_diversity_stats.to_html()}''')


def alpha_diversity(asv_table, map_file, data_column, treatments, plot_title, outputdir):
    pd.options.mode.chained_assignment = None
    #Further resources can be found at the following links below:
    #https://develop.qiime2.org/en/latest/intro.html
    #https://docs.qiime2.org/2024.5/plugins/
    
    #Filter feature table to only contain samples with a tag in the given column
    asv_table_filtered= feature_table.methods.filter_samples(table=asv_table, metadata=map_file, where=f"{data_column} NOT NULL")
    asv_table_filtered = asv_table_filtered.filtered_table
    #NOTE: Does filter properly! 
    
    #Calculate the alpha diversity of each sample 
    alpha_results = diversity.pipelines.alpha(table=asv_table_filtered, metric='shannon')
    alpha_diversity_table = alpha_results.alpha_diversity
    #NOTE: Does calculate alpha diversity values 
    alpha_diversity_table = pd.DataFrame(alpha_diversity_table.view(pd.Series))
    alpha_diversity_table.index.rename('samples',inplace=True)
    """
    The goal here is to the process only given treatments, for example if given CAR then 
    this code block should extract only these samples from the ASV table 

    NOTE: We should first get the alpha diversity table and parse out the samples we want
    """
    print("-------------------------------")
    #asv_table=asv_table.view(pd.DataFrame)
    #asv_table=asv_table.T
    
    treatments=treatments[0].split(',')
    print('Treatments to be processed...')
    for i in range(len(treatments)):
        print(treatments[i], end='\t')
    
    n = len(treatments)
    dataframe_list=[]
    all_samples=alpha_diversity_table.index.to_list()
    for i in range(n):
        #Get the current treatment
        current_treatment=treatments[i]

        #Extract the samples from map file that are labeled with the current treatement
        # *Uses qiime 2 Metdata function 'get_ids' to extract all samples from a treatment based on the map file 
        samples = list(map_file.get_ids(f"[{data_column}]='{current_treatment}'"))
        #Here we will go through each sample and find its corresponding alpha diversity value
        alpha_diversity_score=[]
        for sample in samples:
            if sample not in all_samples:
                print(f"{sample} is not in the ASV table, please check raw counts file for this sequence run")
                samples.remove(sample)
            else:
                alpha_diversity_score.append((sample,alpha_diversity_table.loc[sample,'shannon_entropy']))
        temp_df=pd.DataFrame({'score':[alpha_diversity_score]}, index=[current_treatment])
        temp_df.index.rename('treatment',inplace=True)
        print(f'current df {temp_df}')
        exit(1)
        #Append treatment dataframe to a list of dataframes
        dataframe_list.append(temp_df)
            
    
    #Concate each dataframe from the data frame list by columns
    print(dataframe_list[0])
    exit(1)
    asv_table_filtered=pd.concat(dataframe_list, axis=0, join='outer')

    print("Merged, grouped, and filtered down table...")
    print(asv_table_filtered.columns.to_list())
    exit(1)
    #Run Shannon alpha diversity metric on the filtered table
    alpha_results = diversity.pipelines.alpha(table=asv_table_filtered, metric='shannon')
    alpha_diversity_table = alpha_results.alpha_diversity
    

    #Turn the data into a pandas data frame for further parsing
    alpha_diversity_table = pd.DataFrame(alpha_diversity_table.view(pd.Series))
    alpha_diversity_table_sorted = alpha_diversity_table.sort_values('shannon_entropy', ascending=False)
    
    #Extract ids and turn the dataframe to dictoniary for further parsing.
    #This part could probably be done better, as I don't think you need to pull this information out from the data frame.
    sample_ids=alpha_diversity_table.index.to_list()
    sample_ids_sorted=alpha_diversity_table_sorted.index.to_list()
    shannon_score=alpha_diversity_table.to_dict()
    
    #Get the data associated with the column to parse
    #https://develop.qiime2.org/en/latest/plugins/references/metadata-api.html
    colum = map_file.get_column(f"{data_column}")

    data_dict={}
    stats_dict={}
    #Loop through the list and get the sample's meta tag from the column
    #and assign its Shannon index value to a dictionary where the key is the sample's meta tag and its value 
    # is a list containg all the samples associated with this meta tag and their Shannon index score
    for i in range(len(sample_ids)):
        if(colum.get_value(sample_ids[i]) not in data_dict):
            data_dict[colum.get_value(sample_ids[i])] = []
            data_dict[colum.get_value(sample_ids[i])].append(shannon_score['shannon_entropy'][sample_ids[i]])
        else:
            data_dict[colum.get_value(sample_ids[i])].append(shannon_score['shannon_entropy'][sample_ids[i]])

        if (colum.get_value(sample_ids_sorted[i]) not in stats_dict):
            stats_dict[colum.get_value(sample_ids_sorted[i])] = []
            stats_dict[colum.get_value(sample_ids_sorted[i])].append((sample_ids_sorted[i],shannon_score['shannon_entropy'][sample_ids_sorted[i]]))
            
        else:
            stats_dict[colum.get_value(sample_ids_sorted[i])].append((sample_ids_sorted[i],shannon_score['shannon_entropy'][sample_ids_sorted[i]]))

    visualizer(data_dict, plot_title, outputdir)
    stats_generator(stats_dict, outputdir)


def validate_data(asv_table) -> None:
    #Check if data is a qza type
    if '.qza' in asv_table:
        asv_table=Artifact.load(asv_table)
        return asv_table

    return None


if __name__ == '__main__':
    parser = argparse.ArgumentParser(add_help=False, prog="alpha-diversity-genator.py", description="Program to generate custom alpha diversity boxplots")

    parser.add_argument('-i',"--input-file", required=True, help="Imported feature table",type=str)
    parser.add_argument('-m',"--map-file", required=True, help="Map file for data",type=str)
    parser.add_argument('-c',"--column", required=True, help="Colmun to parse for data",type=str)
    parser.add_argument('-p', "--plot-title", help="Tilte for plot",type=str)
    parser.add_argument('-l', "--listing", nargs='+', type=str, help="Set a preferred listing for x axis (Default is nothing)")
    parser.add_argument('-d', "--output-dir", required=True, help="Output directory location",type=str)
    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Display commands possible with this program.')
    args = parser.parse_args()

    data_file=args.input_file
    map_file=args.map_file
    data_column=args.column
    plot_tilte=args.plot_title
    treatments=args.listing
    output=os.path.join(args.output_dir, "alpha-output/")

    if ((asv_table := validate_data(data_file)) != None) and ((map_file := Metadata.load(map_file)) != None):
        if not os.path.exists(output):
            os.mkdir(output)
        alpha_diversity(asv_table,map_file,data_column,treatments,plot_tilte,output)
    else:
        print('Invalid data type or map file')
        exit(1)
