#!/bin/bash/ env python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 
import argparse
from datetime import datetime
"""
Resource list:
https://saturncloud.io/blog/how-to-divide-multiple-columns-by-another-column-in-pandas/
https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.bar.html
https://stackoverflow.com/questions/41968732/set-order-of-columns-in-pandas-dataframe
https://www.reddit.com/r/learnpython/comments/zl6ggo/how_to_sort_pandas_data_frame_columns_based_on_a/

Not used but might be good to reference
https://stackoverflow.com/questions/22263807/how-is-order-of-items-in-matplotlib-legend-determined
https://www.geeksforgeeks.org/how-to-change-order-of-items-in-matplotlib-legend/
https://www.reddit.com/r/learnpython/comments/8zbbql/changing_order_of_items_in_legend_matplotlibpyplot/


Qiime commands

Used to filter out  samples not labeled with Treatment
    *Used original metadata table
    *Must meet qiime standards of having no empty headers and no empty values (Maybe, still need to test this)
This creates a new ASV table containg only samples with this specified filter
qiime feature-table filter-samples --i-table asv_table_02_add_taxa_L6_norm.qza --m-metadata-file bact-map.txt --p-where "PHASE_3_Liquid_Greenhouse_SOIL IN ('T1Tm0','T2Tm0', 'T3Tm0', 'T4Tm0', 'T5Tm0', 'T6Tm0', 'T7Tm0', 'T8Tm0', 'T9Tm0', 'T10Tm0', 'T1Tm154', 'T2Tm154', 'T3Tm154', 'T4Tm154', 'T5Tm154', 'T6Tm154', 'T7Tm154', 'T8Tm154')" --o-filtered-table treatment_filtered.qza

From there we generating a visualization of the filtered asv table and download the CSV file generated from this filtering
    *This CSV file is what is used in this program
    *If you would like to use a metadata file with this qzv then you'll need a new map file with only these samples present
qiime taxa barplot --i-table treatment_filtered.qza --o-visualization treatment_filtered.qzv


Here are some extra qiime commands used not needed for this program
These were used to group ALL ASVs in the table by a metatag in the table.
For example a table like the one below, would group all samples with the same tag into one ASV
                Treatment
B001            Soil
B002            Root
B003            Soil
B004            Root

Used to make this type of table
    *Again must follow standard qiime metadata table rules
    *ACTUALLY I THINK I did infact use this 
qiime feature-table group --i-table treatment_filtered.qza --p-axis sample --m-metadata-file treatment_filtered_map.tsv --m-metadata-column PHASE_3_Liquid_Greenhouse_SOIL --p-mode sum --o-grouped-table treatment_filtered_group.qza 

Used to make a visulation of the table
qiime taxa barplot --i-table treatment_filtered.qza --m-metadata-file treatment_filtered_map.tsv --o-visualization treatment_filtered.qzv
"""

parser = argparse.ArgumentParser(add_help=False, prog="filter_data.py", description="Program to filter qzv data")

parser.add_argument('-i',"--input-file", required=True, help="CSV file needed containing data to filter",type=str)
parser.add_argument('-n', "--top-n-taxa", required=True, help="Filter for top N taxa",type=int)
parser.add_argument('-f', "--filter", action="store_true", help="Filter out any taxa(Default is viruses)")
parser.add_argument('-l', "--listing", nargs='+', type=str, help="Set a preferred listing for x axis (Default is whatever is listed in the CSV file)")
parser.add_argument('-d', "--output-dir", required=True, help="Output directory location",type=str)
parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Display commands possible with this program.')
args = parser.parse_args()

data_file=args.input_file
n=args.top_n_taxa
index_listing=args.listing
output=args.output_dir

df = pd.read_csv(data_file,index_col=False)

df.set_index('index', inplace = True)

#Reording DF if listing is provided
if index_listing:
    print('Listing provided')
    df=df.reindex(index_listing[0].split(','))


print("=================================")
print("Extracting the top N taxa")
print("=================================")

#Summing all the columns in the DF
abundance = df[df.columns].sum(axis=0)


#Soring DF by Greatest to least
taxa = abundance.sort_values(ascending=False)
#Extracting the top N taxa and their columns
top_n = taxa.head(n)
top_n_columns = top_n.index.tolist()


#Might need to add ability to keep looping until we find a replace for all unwanted taxa
#Same process honestly, just with a while looping.
#Keep a running list of all taxa to remove and keep looping until we have enough new taxa to replace old taxa
    #Also save their index (Tuple, mayhaps)
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

time_generated=datetime.now().strftime("%d/%m/%y %H:%M:%S")
with open(f'{output}/top_n_stats.txt', "w") as f:
    print(f"==========================================",file=f)
    print(f"Top_N_Stats", file=f)
    print(f"==========================================",file=f)
    print(f"Date file was generated:{time_generated}\n",file=f)
    print(top_n,file=f)
    print(f"==========================================",file=f)
    print('Taxa filtered out',file=f)
    print(taxa_removed,file=f)
    

print("=================================")
print("Creating Other DF")
print("=================================")

#Creating a new DF by dropping the top N columns from the original DF
other_df=df.drop(top_n_columns, axis=1)


#Extracting the headers from this DF

old_headers = other_df.columns.tolist()
#Creating a new column titled other that has total amount of taxa present in each sample
other_df['Other']=other_df[other_df.columns].sum(axis=1)
#Dropping all other headers
other_df=other_df.drop(old_headers, axis=1)
print("=================================")
print("Creating Top N Taxa DF")
print("=================================")

#Creating a new DF by dropping the columns not in the Top N from the original DF
top_taxa_df=df.drop(old_headers,axis=1)
#Creating a new column order to reindex top_taxa_by
new_col_order = []


top_n_columns.reverse()
new_col_order+=top_n_columns

top_taxa_df=top_taxa_df.reindex(columns=new_col_order)
print("=================================")
print("Combining Top N DF and Other DF")
print("=================================")

#Merging df containing taxas grouped under other
#and top N taxa
combined_df=pd.merge(top_taxa_df,other_df, left_index=True, right_index=True)
print(combined_df)
#Caluclating the total of each row
total=combined_df[combined_df.columns].sum(axis=1)

#Extracting just the data columns from the list
data_columns=combined_df.columns.tolist()

#Diving each row in a column by the total amount and multiplying it by 100
combined_df[data_columns]=combined_df[data_columns].div(total,axis=0)
combined_df[data_columns]=combined_df[data_columns].multiply(100,axis=0)
headers = combined_df.columns[:len(combined_df.columns)-1].tolist()

#Extracting header columns while ignoring the 'Other' column
headers = combined_df.columns[:len(combined_df.columns)-1].tolist()

#Formatting column headers
for i in range(len(top_n_columns)):
    if 'g_' not in top_n_columns[i]:
        if 'f_' in top_n_columns[i]:
            list = top_n_columns[i].split(';')
            for j in range(len(list)):
                if 'f__' in list[j]:
                    top_n_columns[i] = list[j]
        else:
            top_n_columns[i] = 'Unknown Bacteria'
    else:
        top_n_columns[i] = top_n_columns[i].split(';')[-1]


#Generating Barplot
fig, ax = plt.subplots(figsize = (15, 10))

#https://stackoverflow.com/questions/16006572/plotting-different-colors-in-matplotlib
cmap = plt.get_cmap('tab20')
colors = [cmap(i) for i in np.linspace(0, 1, len(headers))]

#ax.bar(combined_df["index"], combined_df["Other"], bottom=0, width=0.9, color='black', alpha=0.8, label='Other')
ax.bar(combined_df.index, combined_df["Other"], bottom=0, width=0.9, color='black', alpha=0.8, label='Other')
current = combined_df['Other']

#Plotting each column from the combined_df into the barplot
for i in range(len(headers)):
    #ax.bar(combined_df["index"], combined_df[headers[i]], bottom=current, width=0.9, label=top_n_columns[i], color=colors[i])
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
plt.tight_layout()
plt.savefig(f'{output}/Taxa_barplot.png')