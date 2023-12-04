## usage > python Script_1_QC_Normalization.py --input experssion_fromQupath.csv --output anndata_qced_normalized.h5ad --exclude_core "core1,core2" --exclude_marker "marker1,marker2"



import argparse
from os import listdir, path, getcwd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import harmonypy
import warnings
from joblib import parallel_backend
import anndata as ad
import scanpy as sc
import pytometry as pm
import scanpy.external as sce
import os

warnings.filterwarnings('ignore')

def main(args):
    # Read CSV file
    df = pd.read_csv(args.input)

    # clean table
    data = df[list(filter(lambda x: ': Cell: Median' in x, set(df)))]
    meta = df[['Image', 'Object ID', 'Name', 'Class', 'Parent', 'ROI', 'Centroid X µm', 'Centroid Y µm',
               'Nucleus: Area µm^2', 'Nucleus: Length µm', 'Nucleus: Circularity', 'Nucleus: Solidity',
               'Nucleus: Max diameter µm', 'Nucleus: Min diameter µm', 'Cell: Area µm^2', 'Cell: Length µm',
               'Cell: Circularity', 'Cell: Solidity', 'Cell: Max diameter µm', 'Cell: Min diameter µm',
               'Nucleus/Cell area ratio']]
    data.columns = data.columns.str.replace(": Cell: Median", "")
    meta = meta.rename(columns={'Centroid X µm': "x", 'Centroid Y µm': "y"})


    # make anndata

    adata = ad.AnnData(data)
    adata.obs = meta
    spatial = pd.DataFrame(adata.obs[['x', 'y']])
    spatial = spatial.to_numpy()
    adata.obsm["spatial"] = spatial

    adata.obs['Image'] = adata.obs['Image'].str.replace(' ', '-')

    # QC -1
    image_counts = adata.obs['Image'].value_counts()
    adata = adata[adata.obs['Image'].isin(image_counts[image_counts > 500].index)].copy()

    # Exclude user-specified markers
    excluded_core = args.exclude_core.split(',')
    adata = adata[~adata.obs['Image'].isin(excluded_core)].copy()

    # QC -2
    #sc.pl.violin(adata, keys=['DAPI'])
    DAPI_threshold = np.percentile(adata[:, 'DAPI'].X, 10)
    adata = adata[adata[:, 'DAPI'].X > DAPI_threshold, :]

    # QC -3
    adata = adata[(adata.obs['Nucleus: Area µm^2'] >= adata.obs['Nucleus: Area µm^2'].quantile(0.01)) &
                  (adata.obs['Nucleus: Area µm^2'] <= adata.obs['Nucleus: Area µm^2'].quantile(0.99))]

    # QC -4
    # Exclude user-specified genes
    exclude_markers = args.exclude_marker.split(',')
    adata = adata[:, ~adata.var_names.isin(exclude_markers)]

    # Normalization and Zscale
    pm.tl.normalize_arcsinh(adata, cofactor=150, inplace=True)
    X = adata.X.copy()
    X = stats.zscore(X, axis=0)
    X = stats.zscore(X, axis=1)
    adata.X = X

    # PCa and PCa harmony
    sc.tl.pca(adata)
    with parallel_backend('threading', n_jobs=32):
        sc.external.pp.harmony_integrate(adata, key='Image', max_iter_harmony=40)

    # clustering
    random_seed = 49
    sc.settings.seed = random_seed
    sce.tl.phenograph(adata, clustering_algo='leiden', resolution_parameter=0.5, n_jobs=-1, k=20)
    sc.tl.rank_genes_groups(adata, 'pheno_leiden', method='t-test')

    # Write ouput
    del adata.obsp["pheno_jaccard_ig"]
    output_dir = "output_QC"
    os.makedirs(output_dir, exist_ok=True)
    output_file_path = os.path.join(output_dir, args.output)
    adata.write(filename=output_file_path)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process spatial data.")
    parser.add_argument("--input", required=True, help="Input CSV file path")
    parser.add_argument("--output", required=True, help="Output h5ad file name")
    parser.add_argument("--exclude_core", default="", help="Comma-separated list of cores to exclude during QC")
    parser.add_argument("--exclude_marker", default="", help="Comma-separated list of markers to exclude during QC")
    args = parser.parse_args()
    main(args)
