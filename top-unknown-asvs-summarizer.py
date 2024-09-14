import argparse
import pandas as pd
import os
import biom

parser = argparse.ArgumentParser(add_help=False, prog="taxa-bar-genator.py", description="Program to generate custom taxaonmy barplots")
parser.add_argument('-i',"--input-file", required=True, help="imported qza file",type=str)
parser.add_argument('-n', "--top-n-taxa", required=True, help="Filter for top N taxa",type=int)
#parser.add_argument('-d', "--output-dir", required=True, help="Output directory location",type=str)
parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Display commands possible with this program.')
args = parser.parse_args()

file  = args.input_file
n = args.top_n_taxa
#output=os.path.join(args.output_dir,"unknown-asvs-output/")

#if not os.path.exists(output):
#    os.mkdir(output)
#print(f"Output directory: {output}")

biom_table = biom.load_table(file)


