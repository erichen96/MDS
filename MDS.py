from Bio import SeqIO
from Bio import pairwise2
from sklearn.cluster import AffinityPropagation, KMeans, DBSCAN
from sklearn.manifold import MDS, t_sne, Isomap
from sklearn.metrics import silhouette_score
# from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
from IPython.display import Image
from numpy import array
from sgt import Sgt

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# pairwise2 returns scores greater than 100, which causes divide by zero errors.
sequences = [s for s in SeqIO.parse('all.fasta', 'fasta')]
matrix = [[0] * (len(sequences)-1)] * (len(sequences)-1)
matrix2 = [[0] * (len(sequences)-1)] * (len(sequences)-1)
for i in range((len(sequences)) - 1):
    for j in range((len(sequences)) - 1):
        matrix[i][j] = pairwise2.align.globalxx((str(sequences[i].seq)), (str(sequences[j].seq)), score_only=True)
        matrix2[i][j] = pairwise2.align.localxx((str(sequences[i].seq)), (str(sequences[j].seq)), score_only=True)

print(matrix)
print(matrix2)

mds = MDS()
mds_coords = mds.fit_transform(matrix)
mds_coords