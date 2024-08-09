import argparse
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from qiime2 import Artifact
from qiime2 import Metadata
from qiime2.plugins import feature_table
from qiime2.plugins.taxa.visualizers import barplot
import numpy as np 
import os

#Further resources can be found at the following links below:
#https://develop.qiime2.org/en/latest/intro.html
#https://docs.qiime2.org/2024.5/plugins/

parser = argparse.ArgumentParser(add_help=False, prog="taxa-bar-genator.py", description="Program to generate custom taxaonmy barplots")
parser.add_argument('-i',"--input-file", required=True, help="imported qza file",type=str)
parser.add_argument('-m',"--map-file", required=True, help="Map file for data",type=str)
parser.add_argument('-c',"--column", required=True, help="Colmun to parse for data",type=str)
parser.add_argument('-p', "--plot-title", help="Tilte for plot",type=str)
parser.add_argument('-n', "--top-n-taxa", required=True, help="Filter for top N taxa",type=int)
parser.add_argument('-f', "--filter", action="store_true", help="Filter out any taxa (Default is viruses)")
parser.add_argument('-t', "--data-type", required=True, help="Type of sequence data program will be working with\nf = Fungal\nb=Bacteria", type=str)
parser.add_argument('-l', "--listing", nargs='+', type=str, help="Set a preferred listing for x axis (Default is nothing)")
parser.add_argument('-d', "--output-dir", required=True, help="Output directory location",type=str)
parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Display commands possible with this program.')
args = parser.parse_args()

data_file=args.input_file
map_file=args.map_file
data_column=args.column
n=args.top_n_taxa
index_listing=args.listing
output=os.path.join(args.output_dir,"taxanomic-output/")
seq_type=args.data_type
title=args.plot_title

if not os.path.exists(output):
    os.mkdir(output)
print(f"Output directory: {output}")

#Load in the map file and artifact file. 
#This is all specific to qiime2, and in order to work with the qiime2 API, these needed to be loaded in
artifact = Artifact.load(data_file)
print("Artifact loaded..")

#Depending on the version of panads a user has, a warning will be thrown if the map file has any empty cells
map_file=Metadata.load(map_file)
print("Meta data loaded...")


#Qiime2 filter feature table(ASV table) by samples command
#For more information, refer to the qiime2 documentation.
#https://docs.qiime2.org/2024.5/plugins/available/feature-table/filter-features/
asv_table_filtered= feature_table.methods.filter_samples(table=artifact, metadata=map_file, where=f"{data_column} NOT NULL")
print("Table filtered...")
asv_table_filtered = asv_table_filtered.filtered_table

#Extract the column that data will be grouped by
colum = map_file.get_column(f"{data_column}")

#Qiime2 group feature table(ASV table) that have the same meta tag in the map file
#For more options to test refer to the qiime2 documentation
#https://docs.qiime2.org/2024.5/plugins/available/feature-table/group/
asv_table_grouped_results = feature_table.methods.group(table=asv_table_filtered, axis='sample', metadata=colum, mode='sum')
asv_table_grouped = asv_table_grouped_results.grouped_table

#Creating qiime2 visualization file
print("Generating raw data as a qiime 2 qzv file...")
asv_table_grouped_qzv = barplot(table=asv_table_grouped)
asv_table_grouped_qzv = asv_table_grouped_qzv.visualization
asv_table_grouped_qzv.save(f"{output}Taxa_barplot_qiime")


#Take the grouped ASV table and convert it to a pandas Data Frame
df = asv_table_grouped.view(pd.DataFrame)

#Reording DF if listing is provided
if index_listing:
    print('Listing provided')
    df=df.reindex(index_listing[0].split(','))


print("=================================")
print("Extracting the top N taxa")
print("=================================")

#Summing all the columns in the DF to get a taxaonmy's overall abundance
abundance = df[df.columns].sum(axis=0)


#Sorting this DF by Greatest to least
taxa = abundance.sort_values(ascending=False)
#Extracting the top N taxa and their column headindgs
top_n = taxa.head(n)
top_n_columns = top_n.index.tolist()


#Currently, we just filter out any unwanted taxa from the top N list, such as viruses, which could be a result of phix containments.
#What this loop should do is basically loop through the top N list of taxa present and if an unwanted taxa is present.
#then it will just replace this taxa with the next abundant taxa that is considered "wanted".

#If the next tax in the top N list isn't a wanted tax, then it should keep looping until it finds
#one. It should do this for as many unwanted taxa present as possible, though I haven't encountered a situation like this, so 
#I haven't finalized the logic to address such a situation. The logic in theory seems simple enough to implement; I just need to get around to doing so.
if args.filter == True:
    counter=1
    taxa_removed =[]
    for i in range(len(top_n_columns)):
        if 'k__Virus' in top_n_columns[i]:
            print('Top N contains an unwanted taxa')
            new_top_n=taxa.head(n+1).index.tolist()
            if 'k_Virus' in new_top_n[-1]:
                counter+=1
            else:
                print('Found new taxa to replace previous')
                taxa_removed.append(top_n_columns[i])
                new_top_n.remove(top_n_columns[i])
                top_n_columns=new_top_n
                break
else:
    print('No filtering')
    
print("=================================")
print("Creating Other DF")
print("=================================")

#Creating a new DF by dropping the top N columns from the original DF
other_df=df.drop(top_n_columns, axis=1)


#Extracting the headers from this DF
old_headers = other_df.columns.tolist()

#Creating a new column titled other that has total amount of taxa present in each sample that aren't in the top N list
other_df['Other']=other_df[other_df.columns].sum(axis=1)

#Dropping all the taxa not in the top N list so the data frame only contains the "Other" column
other_df=other_df.drop(old_headers, axis=1)
print("=================================")
print("Creating Top N Taxa DF")
print("=================================")

#Creating a new DF by dropping the columns not in the Top N from the original DF
top_taxa_df=df.drop(old_headers,axis=1)
#Creating a new column order to reindex top_taxa_by
new_col_order = []

#We reserve the list to start from least abundant to most abundant to make it easier to graph
top_n_columns.reverse()
new_col_order+=top_n_columns

top_taxa_df=top_taxa_df.reindex(columns=new_col_order)
print("=================================")
print("Combining Top N DF and Other DF")
print("=================================")

#Merging df containing taxas grouped under other and top N taxa
combined_df=pd.merge(top_taxa_df,other_df, left_index=True, right_index=True)
#Caluclating the total percentage of taxa in each row
total=combined_df[combined_df.columns].sum(axis=1)

#Extracting column headings from the data frame to determine what data to manipulate later on 
data_columns=combined_df.columns.tolist()

#Dividing each row in a column by the total amount of taxa present in a row and multiplying each value by 100
#to get its percentage
combined_df[data_columns]=combined_df[data_columns].div(total,axis=0)
combined_df[data_columns]=combined_df[data_columns].multiply(100,axis=0)
headers = combined_df.columns[:len(combined_df.columns)-1].tolist()
#Extracting header columns while ignoring the 'Other' column
headers = combined_df.columns[:len(combined_df.columns)-1].tolist()

#Pulling out stats and reogrinizing them to be printed to stats file
data_stats = combined_df.T
data_stats_index = data_stats.index.to_list()
data_stats_index = data_stats_index[:len(data_stats_index)-1]
data_stats_index.reverse()
data_stats_index.append('Other')
data_stats=data_stats.reindex(data_stats_index)

time_generated=datetime.now().strftime("%d/%m/%y %H:%M:%S")
with open(f'{output}top_n_stats.txt', "w") as f:
    print(f"==========================================",file=f)
    print(f"Top_N_Stats", file=f)
    print("To find further sequence specific information, refer to table 03 generated previously.", file=f)
    print("Please refer to the excel file generated to perform further analysis.", file=f)
    print(f"==========================================",file=f)
    print(f"Date file was generated:{time_generated}\n",file=f)
    print(data_stats.to_markdown(),file=f)
    print(f"\n==========================================",file=f)
    print('Taxa filtered out',file=f)
    if(args.filter == True):
        print(taxa_removed,file=f)

#Saving to excel
data_stats.to_excel(f'{output}top_n_stats.xlsx')


#Formatting column headers
for i in range(len(top_n_columns)):
    if 'g_' not in top_n_columns[i]:
        if 'f_' in top_n_columns[i]:
            list = top_n_columns[i].split(';')
            for j in range(len(list)):
                if 'f__' in list[j]:
                    top_n_columns[i] = list[j]
        elif 'o_' in top_n_columns[i]:
            list = top_n_columns[i].split(';')
            for j in range(len(list)):
                if 'o__' in list[j]:
                    top_n_columns[i] = list[j]
        elif 'c_' in top_n_columns[i]:
            list = top_n_columns[i].split(';')
            for j in range(len(list)):
                if 'c__' in list[j]:
                    top_n_columns[i] = list[j]
        elif 'p_' in top_n_columns[i]:
            list = top_n_columns[i].split(';')
            for j in range(len(list)):
                if 'p__' in list[j]:
                    top_n_columns[i] = list[j]
        else:
            if seq_type == 'b':
                top_n_columns[i] = 'Bacterial ASV'
            elif seq_type == 'f':
                top_n_columns[i] = 'Fungal ASV'
    else:
        top_n_columns[i] = top_n_columns[i].split(';')[-1].split('_')[-1]


#Generating Barplot
fig, ax = plt.subplots(figsize = (15, 10))

#Generating color map for bar plot
#https://stackoverflow.com/questions/16006572/plotting-different-colors-in-matplotlib
cmap = plt.get_cmap('tab20')
colors = [cmap(i) for i in np.linspace(0, 1, len(headers))]

#Plot the "Other" column first.
ax.bar(combined_df.index, combined_df["Other"], bottom=0, width=0.9, color='black', alpha=0.8, label='Other')
current = combined_df['Other']

#Plotting each column from the combined_df into the barplot by looping through
#each top N list which starts from least abundant to most abundant
for i in range(len(headers)):
    ax.bar(combined_df.index, combined_df[headers[i]], bottom=current, width=0.9, label=top_n_columns[i], color=colors[i])
    current+=combined_df[headers[i]]

#Settings for barplot
plt.xticks(rotation=70,fontsize='10')
plt.yticks(fontsize='15')
plt.ylabel("Relative Abundance %", fontsize='15')
handles, labels = plt.gca().get_legend_handles_labels()
labels.reverse()
handles.reverse()


#https://stackoverflow.com/questions/15637961/matplotlib-alignment-of-legend-title
#https://stackoverflow.com/questions/4700614/how-to-put-the-legend-outside-the-plot
ax.legend(handles, labels, bbox_to_anchor=(1, 1), frameon=False, title="Genus", alignment='left')

#https://stackoverflow.com/questions/12402561/how-to-set-font-size-of-matplotlib-axis-legend
plt.setp(plt.gca().get_legend().get_texts(), fontsize='15')
plt.setp(plt.gca().get_legend().get_title(),fontsize='20')

#https://matplotlib.org/stable/gallery/text_labels_and_annotations/titles_demo.html
ax.set_title(title,fontsize='20')

#https://stackoverflow.com/questions/14908576/how-to-remove-frame-from-a-figure
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(f'{output}Taxa_barplot.png')