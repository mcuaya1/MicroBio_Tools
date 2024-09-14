import argparse
from datetime import datetime
import os
from qiime2 import Metadata
from qiime2 import Artifact
from qiime2.plugins import feature_table
from qiime2.plugins.taxa.visualizers import barplot
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd


#To clean up asv labels
def asv_label_formatter(asv_list):
    for i in range(len(asv_list)):
        if 'Other' in asv_list[i]:
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


def visualizer(top_taxa_table, plot_title, outputdir):
    treatment_total=top_taxa_table[top_taxa_table.columns].sum(axis=1)
    
    print('Values to be used to normalize')
    print(treatment_total)
    
    #Normalize values
    top_taxa_table=top_taxa_table.div(treatment_total,axis=0)
    print("Normalizing table...")
    print(top_taxa_table.T)
    
    top_taxa_table=top_taxa_table[top_taxa_table.columns].multiply(100,axis=0)
    print("Calculating percentage...")
    print(top_taxa_table.T)
    
    
        
    headers = top_taxa_table.columns.to_list()
    
    fig, ax = plt.subplots(figsize = (15, 10))
    
    #Generating color map for bar plot
    #https://stackoverflow.com/questions/16006572/plotting-different-colors-in-matplotlib
    
    cmap = plt.get_cmap('tab20')
    colors = [cmap(i) for i in np.linspace(0, 1, len(headers))]
    
    #Plot the "Other" column first.
    ax.bar(top_taxa_table.index, top_taxa_table["Other"], bottom=0, width=0.9, color='black', alpha=0.77, label='Other')
    current = top_taxa_table['Other']
    
    headers.pop(len(headers)-1)
    i=len(headers)-1
    
    #Loop through top N ASVs and plot them according to least abundant to most
    while i >= 0:
        ax.bar(top_taxa_table.index, top_taxa_table[headers[i]], bottom=current, width=0.9, label=headers[i], color=colors[i])
        current+=top_taxa_table[headers[i]]
        i-=1
    
    #Settings for barplot
    plt.xticks(rotation=90,fontsize='15')
    plt.yticks(fontsize='15')
    plt.ylabel("Relative Abundance %", fontsize='15')
    handles, labels = plt.gca().get_legend_handles_labels()
    labels.reverse()
    handles.reverse()


    #https://stackoverflow.com/questions/15637961/matplotlib-alignment-of-legend-title
    #https://stackoverflow.com/questions/4700614/how-to-put-the-legend-outside-the-plot
    ax.legend(handles, labels, bbox_to_anchor=(1, 1), frameon=False, title="ASV", alignment='left')

    #https://stackoverflow.com/questions/12402561/how-to-set-font-size-of-matplotlib-axis-legend
    plt.setp(plt.gca().get_legend().get_texts(), fontsize='15')
    plt.setp(plt.gca().get_legend().get_title(),fontsize='20')

    #https://matplotlib.org/stable/gallery/text_labels_and_annotations/titles_demo.html
    ax.set_title(f"{plot_title}",fontsize='20')

    #https://stackoverflow.com/questions/14908576/how-to-remove-frame-from-a-figure
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.tight_layout()
    print("Saving visualization...")
    fig.savefig(f'{outputdir}{plot_title}.png')

def stats_generator(asv_table: pd.DataFrame, outputdir: str, method:str, raw_asv_strings: list):
    time_generated=datetime.now().strftime("%d/%m/%y %H:%M:%S")
    
    #Normalize results
    treatment_total=asv_table[asv_table.columns].sum(axis=0)
    asv_table=asv_table.div(treatment_total,axis=1)
    
    asv_table.index = raw_asv_strings
    
    print('Generating excel file...')
    asv_table.T.to_excel(f'{outputdir}top_n_stats.xlsx')
    
    asv_table_normalized=asv_table.multiply(100, axis=1)
    print('Generating markdown file with table stats...')
    with open(f'{outputdir}top_n_stats.md', "w") as f:
        f.write(f'''# Top N stats\n## Method used: {method}\n## To find further sequence specific information, refer to table 03 generated previously\n**Please refer to the excel or csv file generated to perform further analysis.**\nDate file was generated: {time_generated}\n{asv_table_normalized.to_markdown()}''')
    
    print('Generating html file with table stats...')
    with open(f'{outputdir}top_n_stats.html', "w") as f:
        f.write(f'''<!doctype html>
    <html lang="en">
        <head>
            <meta charset="utf-8">
            <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css">
        </head>
        <body>
            <h1>Top N Stats</h1>
            <h2>Method used: {method}</h2>
            <h2 >To find further sequence specific information, refer to table 03 generated previously.</h2>
            <strong>Please refer to the excel file generated to perform further analysis. </strong>
            <p>Date file was generated: {time_generated}</p>
            {asv_table_normalized.to_html()}''')

def borneman_prism_formatter(asv_table, map_file: Metadata, data_column: str, treatments: list, num: int, outputdir: str):
    print("BORNEMAN PRISM FORMATTER")
    pd.options.mode.chained_assignment = None
    
    #Ensure correct data format
    if isinstance(asv_table, pd.DataFrame):
        print('Table to be processed is a txt biom file')
        asv_table=asv_table.set_index('#OTU ID')
        asv_table.index.name = None
    
    elif isinstance(asv_table, Artifact) and 'FeatureTable[RelativeFrequency]' in str(asv_table.view):
        print('Table to be processed is a Qiime 2 Artifact')
        asv_table=asv_table.T

    else:
        print('Invalid table type...')
        exit(1)
    
    print('Treatments to be processed...')
    for i in range(len(treatments)):
        print(treatments[i], end='\t')
    
    n = len(treatments)
    
    dataframe_list=[]
    for i in range(n):
        
        #Get the current treatment
        current_treatment=treatments[i]

        #Extract the samples from map file that are labeled with the current treatement
        # *Uses qiime 2 Metdata function 'get_ids' to extract all samples from a treatment based on the map file 
        samples = list(map_file.get_ids(f"[{data_column}]='{current_treatment}'"))

        
        #Create a temp dataframe which only contains samples related to current treatment
        temp_df=asv_table[samples]
        
        #Get each ASVs total abundance across all current samples
        # *Creates column in the temp dataframe with these values and label the column by the current treatment
        temp_df[f'{current_treatment}'] = temp_df[temp_df.columns].sum(axis=1)
        
        #Remove samples from temp dataframe by extracting them from the data frame and dropping them
        # *Ensures that 'treatment' column is the only one present
        list_temp = temp_df.columns
        list_temp = list_temp[0:len(list_temp)-1]
        temp_df=temp_df.drop(columns=list_temp)
        
        #Append treatment dataframe to a list of dataframes
        dataframe_list.append(temp_df)
    
    #Concate each dataframe from the data frame list by columns    
    merged_data=pd.concat(dataframe_list, axis=1)
    print("Grouped, and filtered down table...")
    print(merged_data)
    
    counter = 0
    curr_row = 0
    treatments=merged_data.columns.to_list()
    top_n_taxa = []
    
    #Finding top ASVs
    while (counter < num):
        #Extract treatment columns
        for i in range(len(treatments)):
            #Sort the merged data frame by current treatment
            sorted_table=merged_data.sort_values(by=f'{treatments[i]}', ascending=False)
            
            #Get the taxa from the current row
            taxa=sorted_table.iloc[curr_row].name
            
            #Check if taxa is not already in the top taxa list
            if taxa not in top_n_taxa:
                top_n_taxa.append(taxa)
                counter+=1
                
            #If we still don't have the top N taxa then keep looping
            if counter >= num:
                break
        
        #Increase row counter
        curr_row+=1
        
    #Create a 'Other' data frame which only has ASVs not in the top taxa list
    other_df=merged_data.drop(top_n_taxa,axis=0)
    
    #Create top taxa data frame which only has ASVs in the top taxa list
    top_taxa_df=merged_data.drop(other_df.index, axis=0)
    
    #Get the total of 'Other' ASVs for each treatment and create dataframe out of these values
    other_total=other_df[other_df.columns].sum(axis=0)
    other_total=other_total.to_frame().T.rename(index={0: 'Other'})
    
    #Group any duplicate ASVs and sort by the control sample (Currently hard coded to be the first column)
    top_taxa_df=top_taxa_df.groupby(top_taxa_df.index).sum()
    top_taxa_df=top_taxa_df.sort_values(by=top_taxa_df.columns.to_list()[0], ascending=False)
    
    top_n_taxa=top_taxa_df.index.to_list()
    
    #Format ASV lables
    asv_label_formatter(top_n_taxa)    
    top_n_taxa.append("Other")
    
    #Concat 'Other' dataframe to the top taxa dataframe
    top_taxa_df=pd.concat([top_taxa_df, other_total], axis=0, ignore_index=True)
    
    #Assign formatted labeling to dataframe
    top_taxa_df.index = top_n_taxa
    
    print(f"Found top {num} ASVs...")
    print(top_taxa_df)
    top_taxa_df=top_taxa_df.T   
    
    #Save as a csv file
    with open(f"{outputdir}{data_column}.csv","w") as file:
        file.write(top_taxa_df.to_csv())


def qiime_formatter(asv_table: Artifact, map_file: Metadata, data_column: str, output: str):
    #Ensure correct data format
    if 'FeatureTable[Frequency]' in asv_table.view:
        print('Table to be processed is a Qiime 2 Artifact')
    else:
        print('Invalid data type')
        exit(1)
    
    #Filter feature table
    asv_table_filtered = feature_table.methods.filter_samples(table=asv_table, metadata=map_file, where=f"{data_column} NOT NULL")
    asv_table_filtered = asv_table_filtered.filtered_table
    
    #Extract the column that data will be grouped by
    colum = map_file.get_column(f"{data_column}")

    #Group features under the same meta tag according to given map file
    asv_table_grouped_results = feature_table.methods.group(table=asv_table_filtered, axis='sample', metadata=colum, mode='sum')
    asv_table_grouped = asv_table_grouped_results.grouped_table
    
    #Create qiime 2 visualization
    asv_table_grouped_qzv = barplot(table=asv_table_grouped)
    asv_table_grouped_qzv = asv_table_grouped_qzv.visualization
    asv_table_grouped_qzv.save(f"{output}{data_column}")

def biime_formatter(asv_table : Artifact, map_file : Metadata , col ,treatments, num, outputdir, plot_title):
    print('BIIME FORMATTER')
    pd.options.mode.chained_assignment = None
    
    #Ensure correct data format
    if 'FeatureTable[Frequency]' in str(asv_table.view):
        print('Table to be processed is a Qiime 2 Artifact')
    else:
        print('Invalid data type')
        exit(1)

    asv_table=asv_table.view(pd.DataFrame)
    asv_table=asv_table.T

    treatments=treatments[0].split(',')
    print('Treatments to be processed...')
    for i in range(len(treatments)):
        print(treatments[i], end='\t')
    
    n = len(treatments)
    
    dataframe_list=[]
    for i in range(n):
        
        #Get the current treatment
        current_treatment=treatments[i]

        #Extract the samples from map file that are labeled with the current treatement
        # *Uses qiime 2 Metdata function 'get_ids' to extract all samples from a treatment based on the map file 
        samples = list(map_file.get_ids(f"[{col}]='{current_treatment}'"))
        
        
        #Create a temp dataframe which only contains samples related to current treatment
        temp_df=asv_table[samples]
        
        #Get each ASVs total abundance across all current samples
        #*Create column in the temp dataframe with these values and label the column by the current treatment
        temp_df[f'{current_treatment}'] = temp_df[temp_df.columns].sum(axis=1)
        
        #Remove samples from temp dataframe by extracting them from the data frame and dropping them
        # *Ensures that 'treatment' column is the only one present
        list_temp = temp_df.columns
        list_temp = list_temp[0:len(list_temp)-1]
        temp_df=temp_df.drop(columns=list_temp)
        
        
        #Append treatment dataframe to a list of dataframes
        dataframe_list.append(temp_df)
    
    #Concate each dataframe from the data frame list by columns    
    merged_data=pd.concat(dataframe_list, axis=1)

    print("Merged, grouped, and filtered down table...")
    print(merged_data)

    counter = 0
    curr_row = 0
    treatments=merged_data.columns.to_list()
    top_n_taxa = []
    
    #Get the top N taxa
    print(f"Finding top {num} ASVs...")
    while (counter < num):
        
        #Extract treatment columns
        for i in range(len(treatments)):
            #Sort the merged data frame by current treatment
            sorted_table=merged_data.sort_values(by=f'{treatments[i]}', ascending=False)
            
            #Get the taxa from the current row
            taxa=sorted_table.iloc[curr_row].name
            
            
            #Check if taxa is not already in the top taxa list
            if taxa not in top_n_taxa:
                top_n_taxa.append(taxa)
                counter+=1
            
            #If we still don't have the top N taxa then keep looping
            if counter >= num:
                break
    
        #Increase row counter
        curr_row+=1
    
    #Create a 'Other' data frame which only has ASVs not in the top taxa list
    other_df=merged_data.drop(top_n_taxa,axis=0)
    
    
    #Create top taxa data frame which only has ASVs in the top taxa list
    top_taxa_df=merged_data.drop(other_df.index, axis=0)
    
    #Get the total of 'Other' ASVs for each treatment and create dataframe out of these values
    other_total=other_df[other_df.columns].sum(axis=0)
    other_total=other_total.to_frame().T.rename(index={0: 'Other'})

    
    #Sort by the control sample (Currently hard coded to be the first column)
    top_taxa_df=top_taxa_df.sort_values(by=top_taxa_df.columns.to_list()[0], ascending=False)
    
    top_n_taxa=top_taxa_df.index.to_list()
    
    #Extract raw ASV labels for stats file
    raw_asv_strings=top_taxa_df.index.to_list()
    raw_asv_strings.append("Other")
    
    #Format ASV lables
    asv_label_formatter(top_n_taxa)

    top_n_taxa.append("Other")
    
    #Concat 'Other' dataframe to the top taxa dataframe
    top_taxa_df=pd.concat([top_taxa_df, other_total], axis=0, ignore_index=True)

    #Assign formatted labeling to dataframe
    top_taxa_df.index = top_n_taxa
    
    print(f"Found top {num} ASVs...")
    print(top_taxa_df)
    
    print("Generating visualization...")
    visualizer(top_taxa_df.T, plot_title, outputdir)

    

    print("Generating stats files...")
    stats_generator(top_taxa_df, outputdir, 'Beth Raw Counts Method', raw_asv_strings)
    

def validate_data(asv_table) -> None:
    
    #Check if data is a qza type
    if '.qza' in asv_table:
        asv_table=Artifact.load(asv_table)
        return asv_table
    
    #Check if data is biom txt file
    elif '.txt' in asv_table:
        asv_table=pd.read_table(asv_table, comment='~', sep = "\t", skiprows = 1)
        return asv_table

    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(add_help=False, prog="taxa-bar-genator.py", description="Program to generate custom taxaonmy barplots")
    parser.add_argument('-i',"--input-file", required=True, help="Imported qza file",type=str)
    parser.add_argument('-m',"--map-file", required=True, help="Map file for data",type=str)
    parser.add_argument('-c',"--column", required=True, help="Colmun to parse for data",type=str)
    parser.add_argument('-p', "--plot-title", help="Tilte for plot",type=str)
    parser.add_argument('-n', "--top-n-taxa", required=True, help="Filter for top N taxa",type=int)
    parser.add_argument('-f', "--filter", action="store_true", help="Filter out any taxa (Default is viruses)")
    parser.add_argument('-t', "--formatter-type", required=True, help="Type of formatter to process data with\nb = Biime Formatter\nj=Borneman prism formatter\nq=Qiime 2 Formatter", type=str)
    parser.add_argument('-l', "--treatments", nargs='+', type=str, help="Treatments to process")
    parser.add_argument('-d', "--output-dir", required=True, help="Output directory location",type=str)
    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Display commands possible with this program.')
    args = parser.parse_args()

    data_file=args.input_file
    map_file=args.map_file
    data_column=args.column
    n_taxa=args.top_n_taxa
    treatments=args.treatments
    output=os.path.join(args.output_dir,"taxanomic-output/")
    formatter_type=args.formatter_type
    title=args.plot_title
    
    if ((asv_table := validate_data(data_file)) != None) and ((map_file := Metadata.load(map_file)) != None):
        if not os.path.exists(output):
            os.mkdir(output)
        if formatter_type == 'b':
            biime_formatter(asv_table, map_file, data_column, treatments, n_taxa, output, title)
        elif formatter_type == 'j':
            borneman_prism_formatter(asv_table, map_file, data_column, treatments, n_taxa, output)
        elif formatter_type == 'q':
            qiime_formatter(asv_table, map_file, data_column, output)
    
        print(f"Output directory: {output}")
    else:
        print('Invalid data type or map file')
        exit(1)

