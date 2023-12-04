# python Script_4_cellEnrich_InNB.py --input Data_Cell_Clinical_Group_NB.csv --output output_cellInrichmentInNB.csv
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import scanpy as sc
import scanpy.external as sce
import anndata as ad
from joblib import parallel_backend
import warnings
import seaborn as sns
import matplotlib
import statsmodels.api as sm

# Argument parser setup
parser = argparse.ArgumentParser(description='Process data and perform analysis.')
parser.add_argument('--input', type=str, help='Input CSV file name')
parser.add_argument('--output', type=str, help='Output CSV file name')
args = parser.parse_args()

# Use the provided input and output file names
input_file = args.input
output_file = args.output

# Data Loading and Preprocessing
meta = pd.read_csv(input_file, index_col='Unnamed: 0')
meta = meta[meta['CT_final'] != "Artifact"]
ct_freq_1 = meta[["Image", "CT_final"]]

ct_freq = pd.crosstab(index=ct_freq_1['Image'],
                      columns=ct_freq_1['CT_final'],
                      normalize="index")

all_freqs_1 = meta[["neighborhood10", "Image", "CT_final"]]

set(all_freqs_1['Image'])

all_freqs_1["Image"].unique()

all_freqs = pd.crosstab(index=[all_freqs_1['neighborhood10'], all_freqs_1['Image']],
                        columns=all_freqs_1['CT_final'],
                        normalize="index")


all_freqs_2 = all_freqs.reset_index()
all_freqs_2.index = all_freqs_2['Image']

df_extracted = meta[['Image', 'Group']]

df_unique_images = df_extracted.drop_duplicates()

group_counts = df_unique_images['Group'].value_counts()
count_responder = group_counts['Responder']
count_non_responder = group_counts['Non_Responder']

responder_images = df_unique_images[df_unique_images['Group'] == 'Responder']['Image'].tolist()
non_responder_images = df_unique_images[df_unique_images['Group'] == 'Non_Responder']['Image'].tolist()

neighborhood_col = 'neighborhood10'
nbs = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
patients = responder_images + non_responder_images

group = pd.Series([1] * count_responder + [0] * count_non_responder)

cells = set(meta['CT_final'])


def normalize(X):
    arr = np.array(X.fillna(0).values)
    return pd.DataFrame(np.log2(1e-3 + arr / arr.sum(axis=1, keepdims=True)), index=X.index.values,
                        columns=X.columns).fillna(0)


X_cts = normalize(ct_freq.reset_index().set_index('Image').loc[patients, cells])

df_list = []

for nb in nbs:
    cond_nb = all_freqs_2.loc[all_freqs_2[neighborhood_col] == nb, cells].rename(
        {col: col + '_' + str(nb) for col in cells}, axis=1)
    df_list.append(normalize(cond_nb))

X_cond_nb = pd.concat(df_list, axis=1).loc[patients]

X_ct2 = X_cts.reset_index(drop=True)

changes = {}
nbs = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
for col in cells:
    for nb in nbs:
        X = pd.concat([X_ct2[col], group.astype('int'), pd.Series(np.ones(len(group)), index=group.index.values)],
                      axis=1).values
        if col + '_%d' % nb in X_cond_nb.columns:
            Y = X_cond_nb[col + '_%d' % nb].values
            X = X[~pd.isna(Y)]
            Y = Y[~pd.isna(Y)]
            results = sm.OLS(Y, X).fit()
            changes[(col, nb)] = (results.pvalues[1], results.params[1])

dat = (pd.DataFrame(changes).loc[1].unstack())
dat = pd.DataFrame(np.nan_to_num(dat.values), index=dat.index, columns=dat.columns).T.sort_index(ascending=True).loc[:,
      X_cts.columns]
pvals = (pd.DataFrame(changes).loc[0].unstack()).T.sort_index(ascending=True).loc[:, X_cts.columns]

f, ax = plt.subplots(figsize=(12, 6))
g = sns.heatmap(dat, cmap='coolwarm', vmin=-1, vmax=1, cbar=False, ax=ax)

for a, b in zip(*np.where(pvals < 0.1)):
    plt.text(b + .5, a + .55, '*', fontsize=30, ha='center', va='center')
plt.tight_layout()

resolution_value = 400
plt.savefig("EN_total0.05.svg", format="svg", dpi=resolution_value, bbox_inches="tight")

all_freqs_3 = all_freqs.reset_index()

map_group = dict(zip(meta['Image'], meta['Group']))

all_freqs_3['group'] = (
    all_freqs_3['Image']
        .map(map_group)
        .astype('category')
)

all_freqs_4 = all_freqs_3[all_freqs_3["group"] != 'Unknown']

all_freqs_4.to_csv(output_file)
