#### usage > python Script_5_Dist_intraction_scimap.py --data_file data/ counts_table.csv --meta_file data/Data_Cell_Clinical_Group_NB.csv --distance_output_name my_distance_data.csv --interaction_output_name my_interaction_data.csv --anndata_output_name anndata_distance.h5ad

import argparse
import os
import scimap as sm
import scanpy as sc
import pandas as pd
import anndata as ad
import seaborn as sns
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def perform_analysis(adata, output_dir, distance_csv_name, interaction_csv_name, anndata_output_name):
    # Spatial distance calculation
    sm.tl.spatial_distance(adata, x_coordinate='y', y_coordinate='y', z_coordinate=None,
                           phenotype='CT_final', subset=None, imageid='Image', label='spatial_distance')

    # Save spatial distances to CSV
    spatial_distance_csv_path = os.path.join(output_dir, distance_csv_name)
    adata.uns['spatial_distance'].to_csv(spatial_distance_csv_path)

    # Spatial interaction analysis
    adata = sm.tl.spatial_interaction(adata, method='knn', radius=30, pval_method='zscore',
                                      imageid='Image', x_coordinate='x', y_coordinate='y',
                                      label='spatial_interaction_knn', phenotype='CT_final')

    # Save spatial interaction results to CSV
    spatial_interaction_csv_path = os.path.join(output_dir, interaction_csv_name)
    adata.uns['spatial_interaction_knn'].to_csv(spatial_interaction_csv_path)

    # Save the modified adata to H5AD file with the changed filename
    output_path = os.path.join(output_dir, anndata_output_name)  # Changed variable name
    adata.write(filename=output_path)

def main(args):
    if args.data_file and args.meta_file:
        # Load data from CSV files
        data = pd.read_csv(args.data_file)
        meta = pd.read_csv(args.meta_file)

        # convert to anndata
        adata = ad.AnnData(data)
        adata.obs = meta
        spatial = pd.DataFrame(adata.obs[['x', 'y']])
        spatial = spatial.to_numpy()
        adata.obsm["spatial"] = spatial
    elif args.anndata_file:
        # Load anndata from file
        adata = ad.read_h5ad(args.anndata_file)
    else:
        print("Specify either CSV files or AnnData file.")
        return

    # Use the directory of the script as the output directory
    script_directory = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_directory, "output_dist_int")

    # Create the "output" directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Perform analysis and save in the "output" directory
    perform_analysis(adata, output_dir, args.distance_output_name, args.interaction_output_name, args.anndata_output_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process neighborhood analysis data.')
    parser.add_argument('--data_file', type=str, help='Path to the neighborhood counts table CSV file')
    parser.add_argument('--meta_file', type=str, help='Path to the metadata CSV file')
    parser.add_argument('--anndata_file', type=str, help='Path to the AnnData file')
    parser.add_argument('--distance_output_name', type=str, default='dataF_distance.csv', help='Name of the spatial distance CSV file')
    parser.add_argument('--interaction_output_name', type=str, default='dataF_intract.csv', help='Name of the spatial interaction CSV file')
    parser.add_argument('--anndata_output_name', type=str, default='anndata_dist_intact.h5ad', help='Name of the output AnnData H5AD file')

    args = parser.parse_args()

    # Run the main function
    main(args)
