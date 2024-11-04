import argparse
from datetime import datetime
import pandas as pd
from qiime2.plugins import feature_table
from qiime2.plugins import diversity
from qiime2 import Metadata
from qiime2 import Artifact
from qiime2.plugins.diversity.visualizers import alpha_group_significance
import os

#Further resources can be found at the following links below:
#https://develop.qiime2.org/en/latest/intro.html
#https://docs.qiime2.org/2024.5/plugins/

parser = argparse.ArgumentParser(add_help=False, prog="alpha-diversity-genator.py", description="Program to generate custom alpha diversity boxplots")

parser.add_argument('-i',"--input-file", required=True, help="Imported feature table",type=str)
parser.add_argument('-m',"--map-file", required=True, help="Map file for data",type=str)
parser.add_argument('-d', "--output-dir", required=True, help="Output directory location",type=str)
parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Display commands possible with this program.')
args = parser.parse_args()

data_file=args.input_file
map_file=args.map_file
output=os.path.join(args.output_dir, "alpha-output/")

#Creating output directory
if not os.path.exists(output):
    os.mkdir(output)
print(f"Output directory: {output}")

#Load in ASV table and map file
artifact = Artifact.load(data_file)
map_file=Metadata.load(map_file)

#Qiime2 alpha diversity file
alpha_results = diversity.pipelines.alpha(table=artifact, metric='shannon')
alpha_diversity_table = alpha_results.alpha_diversity
alpha_vis_file=alpha_group_significance(alpha_diversity=alpha_diversity_table, metadata=map_file)
alpha_vis_file = alpha_vis_file.visualization
alpha_vis_file.save(f"{output}Taxa_barplot_qiime")
