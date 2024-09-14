#TODO: Finish script
import argparse
from datetime import datetime
import pandas as pd
from skbio import DistanceMatrix
from skbio import OrdinationResults
from qiime2.plugins import feature_table
from qiime2.plugins import diversity
from qiime2 import Metadata
from qiime2 import Artifact
import matplotlib.pyplot as plt


#Further resources can be found at the following links below:
#https://develop.qiime2.org/en/latest/intro.html
#https://docs.qiime2.org/2024.5/plugins/

parser = argparse.ArgumentParser(add_help=False, prog="betsa-diversity-genator.py", description="Program to generate custom beta diversity scatter plots")

parser.add_argument('-i',"--input-file", required=True, help="Imported feature table",type=str)
parser.add_argument('-m',"--map-file", required=True, help="Map file for data",type=str)
parser.add_argument('-c',"--column", required=True, help="Colmun to parse for data",type=str)
parser.add_argument('-p', "--plot-title", help="Tilte for plot",type=str)
parser.add_argument('-l', "--listing", nargs='+', type=str, help="Set a preferred listing for x axis (Default is nothing)")
#parser.add_argument('-d', "--output-dir", required=True, help="Output directory location",type=str)
parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Display commands possible with this program.')
args = parser.parse_args()

data_file=args.input_file
map_file=args.map_file
data_column=args.column
plot_tilte=args.plot_title
index_listing=args.listing
#output=args.output_dir

#Load in ASV table and map file
artifact = Artifact.load(data_file)
map_file=Metadata.load(map_file)

#Filter asv table to include only samples from specified group
asv_table_filtered= feature_table.methods.filter_samples(table=artifact, metadata=map_file, where=f"{data_column} NOT NULL")
asv_table_filtered = asv_table_filtered.filtered_table

#beta_diversity_table = asv_table_filtered.view(pd.DataFrame)


#print(beta_diversity_table)

#Preform Braycurtis metric
beta_results = diversity.pipelines.beta(table=asv_table_filtered, metric='braycurtis')

beta_diversity_table = beta_results.distance_matrix

#https://forum.qiime2.org/t/load-distancematrix-artifact-to-dataframe/11660
#Turn table into a distance matrix 
beta_diversity_dm = beta_diversity_table.view(DistanceMatrix)

pcoa_results = diversity.methods.pcoa(distance_matrix=beta_diversity_table)
pcoa_results = pcoa_results.pcoa

matrix=pcoa_results.view(OrdinationResults)
print(matrix)
#print(type(pcoa_results))

#print(beta_diversity_dm)

beta_diversity_df = matrix.samples

fig, ax = plt.subplots(figsize = (15, 10))
values = beta_diversity_df.values
plt.scatter(values[:,0], values[:,1])
plt.tight_layout()
fig.savefig(f"Beta_diversity.png")

