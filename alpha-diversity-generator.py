import argparse
from datetime import datetime
import pandas as pd
from qiime2.plugins import feature_table
from qiime2.plugins import diversity
from qiime2 import Metadata
from qiime2 import Artifact
import matplotlib.pyplot as plt

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
index_listing=args.listing
output=args.output_dir

#Load in ASV table and map file
artifact = Artifact.load(data_file)
map_file=Metadata.load(map_file)

#Filter feature table to only contain samples with a tag in the given column
asv_table_filtered= feature_table.methods.filter_samples(table=artifact, metadata=map_file, where=f"{data_column} NOT NULL")
asv_table_filtered = asv_table_filtered.filtered_table

#Run Shannon alpha diversity metric on the filtered table
alpha_results = diversity.pipelines.alpha(table=asv_table_filtered, metric='shannon')
alpha_diversity_table = alpha_results.alpha_diversity

#Turn the data into a pandas data frame for further parsing
alpha_diversity_table = pd.DataFrame(alpha_diversity_table.view(pd.Series))
#Extract ids and turn dataframe to dictoniary for further parsing
id_list=alpha_diversity_table.index.to_list()
id_dict=alpha_diversity_table.to_dict()
#Get the data assiocated with the column to parse
colum = map_file.get_column(f"{data_column}")

stats_dict={}
data_dict={}

#Loop through the list and get a id's meta tag from column
#and assign its shannon index value to a dictonary where the key
#is the id's meta tag and its value a list containg all id's shannon index assoicated with
#the meta tag
for i in range(len(id_list)):
    if(colum.get_value(id_list[i]) not in data_dict):
        stats_dict[colum.get_value(id_list[i])] = []
        data_dict[colum.get_value(id_list[i])] = []
        data_dict[colum.get_value(id_list[i])].append(id_dict['shannon_entropy'][id_list[i]])
        stats_dict[colum.get_value(id_list[i])].append((id_list[i], id_dict['shannon_entropy'][id_list[i]]))
    else:
        data_dict[colum.get_value(id_list[i])].append(id_dict['shannon_entropy'][id_list[i]])
        stats_dict[colum.get_value(id_list[i])].append((id_list[i], id_dict['shannon_entropy'][id_list[i]]))

fig, ax = plt.subplots(figsize = (15, 10))
medianprops = dict(linestyle='-.', linewidth=3, color='black')
flierprops = dict(marker='o', markerfacecolor='blue', markersize=7,markeredgecolor='none')
boxprops = dict(facecolor='gray')

plt.boxplot(data_dict.values(), labels=data_dict.keys(), 
            patch_artist=True, 
            boxprops=boxprops, medianprops=medianprops, flierprops=flierprops, 
            widths=0.7)

plt.xticks(rotation=90,fontsize='13')
plt.yticks(fontsize='13')

#Old code to possible assign different colors to each boxplot
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
fig.savefig(f"{output}Alpha_diversity.png")

#Sort data for stats extraction
alpha_diversity_table = alpha_diversity_table.sort_values('shannon_entropy', ascending=False)
alpha_diversity_table = alpha_diversity_table.sort_values('shannon_entropy', ascending=False)

#Extract ids and turn dataframe to dictoniary for further parsing
id_list=alpha_diversity_table.index.to_list()
id_dict=alpha_diversity_table.to_dict()

#Get the data assiocated with the column to parse
colum = map_file.get_column(f"{data_column}")

stats_dict={}
data_dict={}
for i in range(len(id_list)):
    if(colum.get_value(id_list[i]) not in data_dict):
        stats_dict[colum.get_value(id_list[i])] = []
        data_dict[colum.get_value(id_list[i])] = []
        data_dict[colum.get_value(id_list[i])].append(id_dict['shannon_entropy'][id_list[i]])
        stats_dict[colum.get_value(id_list[i])].append((id_list[i], id_dict['shannon_entropy'][id_list[i]]))
    else:
        data_dict[colum.get_value(id_list[i])].append(id_dict['shannon_entropy'][id_list[i]])
        stats_dict[colum.get_value(id_list[i])].append((id_list[i], id_dict['shannon_entropy'][id_list[i]]))

#Get the amount of samples in a treatment
length_of_samples = []
for key, value in stats_dict.items():
  length_of_samples.append(len(value))

#Generating stats
alpha_diversity_stats = pd.DataFrame({'Sample Tag':stats_dict.keys(), 'Amount of samples': length_of_samples,'(Sample, Shannon Index)':stats_dict.values()})
alpha_diversity_stats.index = alpha_diversity_stats['Sample Tag']
alpha_diversity_stats.drop('Sample Tag', axis=1, inplace=True)

#Outputting stats
time_generated=datetime.now().strftime("%d/%m/%y %H:%M:%S")
with open(f'{output}/alpha_diversity_stats.txt', "w") as f:
    print(f"==========================================",file=f)
    print(f"Alpha_Diversity_Stats", file=f)
    print("To find further sequence specific information, refer to table 03 generated previously.", file=f)
    print("Please refer to the excel file generated to perform further analysis.", file=f)
    print(f"==========================================",file=f)
    print(f"Date file was generated:{time_generated}\n",file=f)
    print(alpha_diversity_stats.to_markdown(),file=f)
    print(f"\n==========================================",file=f)

alpha_diversity_stats.to_excel(f'{output}alpha_diversity_stats.xlsx')

