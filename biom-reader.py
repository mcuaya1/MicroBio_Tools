import qiime2
import pandas as pd
from qiime2.plugins import feature_table
from qiime2 import Artifact
from qiime2.plugins.taxa.visualizers import barplot

table = Artifact.load('../Fungal/FUN_L6.qza')

vis = qiime2.
df = table.view(pd.DataFrame)