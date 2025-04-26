#Python imports
import argparse
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
import os
import scipy.stats as stats

#Qiime2 imports
from qiime2.plugins import feature_table
from qiime2.plugins import diversity
from qiime2 import Metadata
from qiime2 import Artifact



def signifcance_test(dataframe, outputdir):
    #https://stackoverflow.com/questions/15943769/how-do-i-get-the-row-count-of-a-pandas-dataframe 
    n = dataframe[dataframe.columns[0]].count()
    dataframe_list=[]

    #Here we just find the statsical signifcance between treatments 
    #The algorithm isn't the most efficent as it is O(n^2) but if we add parallization then we 
    #can probably reduce the time spent calculating when using larger datasets
    for i in range(n-1):
        current_treatment = dataframe.iloc[i].item()
        current_treatment_name = dataframe.index[i]
        grouping=[]
        for j in range(i+1, n):
            jth_treatment = dataframe.iloc[j].item()
            jth_treatment_name = dataframe.index[j]
            kruskal_test= stats.kruskal(current_treatment, jth_treatment)
            grouping.append((jth_treatment_name,kruskal_test[0][1], kruskal_test[1][1]))
        temp_df=pd.DataFrame({'label, H value, P value': [grouping]}, index=[current_treatment_name])
        temp_df.index.rename('Treatment',inplace=True)
        #Append treatment dataframe to a list of dataframes
        dataframe_list.append(temp_df)
            
    
    #Concate each dataframe from the data frame list by columns
    signifcance_table=pd.concat(dataframe_list, axis=0, join='outer')
    return signifcance_table

def visualizer(dataframe, plot_title, outputdir):
    cmap = plt.get_cmap('tab20')
    fig, ax = plt.subplots(figsize = (15, 10))
    medianprops = dict(linestyle='-', linewidth=1.5, color='black')
    flierprops = dict(marker='o', markerfacecolor='blue', markersize=7,markeredgecolor='none')
    boxprops = dict(facecolor=cmap(0.09))
   
    data = dataframe.drop(columns=['labeled-scores'])
    plt.boxplot(data['raw-scores'],
                labels=data.index,
                patch_artist=True, 
                boxprops=boxprops, medianprops=medianprops, flierprops=flierprops, 
                widths=0.7)

    plt.xticks(rotation=90,fontsize='13')
    plt.yticks(fontsize='13')

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
    
    dataframe = stats.drop(columns=['raw-scores'])

    datatframe_stats = signifcance_test(dataframe, outputdir) 
    #Saving them to an execel file
    print('Generating excel file...')
    dataframe.to_excel(f'{output}alpha_diversity_stats.xlsx')

    print('Generating markdown file with table stats...')
    time_generated=datetime.now().strftime("%d/%m/%y %H:%M:%S")
    with open(f'{output}alpha_diversity_stats.md', "w") as f:
        f.write(f'''#Alpha diversity stats\n
                ## To find further sequence specific information, refer to table 03 generated previously\n
                **Please refer to the excel or csv file generated to perform further analysis.**\n
                Date file was generated: {time_generated}\n
                {dataframe.to_markdown()}\n
                {datatframe_stats.to_markdown()}''')
        
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
            {dataframe.to_html()}
            {datatframe_stats.to_html()}
            ''')


def alpha_diversity(asv_table, map_file, data_column, treatments, plot_title, outputdir):
    pd.options.mode.chained_assignment = None
    #Further resources can be found at the following links below:
    #https://develop.qiime2.org/en/latest/intro.html
    #https://docs.qiime2.org/2024.5/plugins/
    
    #Filter feature table to only contain samples with a tag in the given column
    asv_table_filtered= feature_table.methods.filter_samples(table=asv_table, metadata=map_file, where=f"{data_column} NOT NULL")
    asv_table_filtered = asv_table_filtered.filtered_table
    
    #Calculate the alpha diversity of each sample 
    alpha_results = diversity.pipelines.alpha(table=asv_table_filtered, metric='shannon')
    alpha_diversity_table = alpha_results.alpha_diversity
    alpha_diversity_table = pd.DataFrame(alpha_diversity_table.view(pd.Series))
    alpha_diversity_table.index.rename('samples',inplace=True)
    treatments=treatments[0].split(',')
    print('Treatments to be processed...')
    for i in range(len(treatments)):
        print(treatments[i], end='\t')
    n = len(treatments)
    dataframe_list=[]
    all_samples=alpha_diversity_table.index.to_list()
    #Here we are mapping treatments and samples together for data visualization and parsing later on 
    for i in range(n):
        #Get the current treatment
        current_treatment=treatments[i]

        #Extract the samples from map file that are labeled with the current treatement
        # *Uses qiime 2 Metdata function 'get_ids' to extract all samples from a treatment based on the map file 
        samples = list(map_file.get_ids(f"[{data_column}]='{current_treatment}'"))
        #Here we will go through each sample and find its corresponding alpha diversity value
        alpha_diversity_score=[]
        raw_scores=[]
        for sample in samples:
            if sample not in all_samples:
                print(f"{sample} is not in the ASV table, please check raw counts file for this sequence run")
            else:
                a_score = (sample,alpha_diversity_table.loc[sample,'shannon_entropy'])
                r_score = alpha_diversity_table.loc[sample,'shannon_entropy']
                alpha_diversity_score.append(a_score)
                raw_scores.append(r_score)

        #https://stackoverflow.com/questions/9376384/sort-a-list-of-tuples-depending-on-two-elements
        #Here we sort the labeled score by their shaonnon score, this is really just for easing viewing 
        if len(alpha_diversity_score) > 0:
            alpha_diversity_score=sorted(alpha_diversity_score, key=lambda scores: scores[-1])
        else:
            alpha_diversity_score.append(0)
            raw_scores.append(0)
        temp_df=pd.DataFrame({'labeled-scores':[alpha_diversity_score], 'raw-scores': [raw_scores]}, index=[current_treatment])
        temp_df.index.rename('treatment',inplace=True)
        #Append treatment dataframe to a list of dataframes
        dataframe_list.append(temp_df)
            
    
    #Concate each dataframe from the data frame list by columns
    asv_table_filtered=pd.concat(dataframe_list, axis=0, join='outer')

    print("Merged, grouped, and filtered down table...")
    print(asv_table_filtered)
    visualizer(asv_table_filtered, plot_title, outputdir)
    stats_generator(asv_table_filtered, outputdir)


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
