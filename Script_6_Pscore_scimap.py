# python Script_6_Pscore_scimap.py --anndata output_dist_int/anndata_distance.h5ad --output_andata_name output_adata_pscore.h5ad --output_volume_name merged_volume.csv --output_density_name merged_density.csv

import os
import argparse
import pandas as pd
from anndata import read_h5ad
import scimap as sm

def calculate_spatial_pscore(adata, cells):
    for i_index, i_item in enumerate(cells):
        for j_index in range(i_index + 1, len(cells)):
            j_item = cells[j_index]
            sm.tl.spatial_pscore(adata, proximity=[i_item, j_item], score_by='Image', x_coordinate='x',
                                 y_coordinate='y', phenotype='CT_final', method='knn', radius=20, knn=3,
                                 imageid='Image', subset=None, label=f'spatial_pscore_{i_item}__{j_item}')
    
    return adata

def save_spatial_scores(adata, cells, output_pscore_dir):
    for i_index, i_item in enumerate(cells):
        for j_index in range(i_index + 1, len(cells)):
            j_item = cells[j_index]
            key = f'spatial_pscore_{i_item}__{j_item}'
            
            if key in adata.uns:
                adata.uns[key].to_csv(os.path.join(output_pscore_dir, f"{i_item}_{j_item}.csv"))
            else:
                print(f"Warning: {key} not found in adata.uns. Skipping.")

def process_csv_files(file_list, direct):

    csv_files = [file for file in file_list if file.endswith('.csv')]

    
    selected_data_frames_v = []
    selected_data_frames_d = []
    
    unique_columns_v = set()
    unique_columns_d = set()
    
    for csv_file in sorted(csv_files):
        file_path = os.path.join(direct, csv_file)
        
        try:
            df = pd.read_csv(file_path, index_col='Image')
        except pd.errors.EmptyDataError:
            print(f"Warning: Empty CSV file at {file_path}")
            continue
        
        new_name_v = f"{csv_file}_Proximity Volume"
        new_name_d = f"{csv_file}_Proximity Density"
        
        df.rename(columns={'Proximity Volume': new_name_v, 'Proximity Density': new_name_d}, inplace=True)
        
        selected_columns_v = [col for col in df.columns if col.endswith('Proximity Volume')]
        selected_columns_d = [col for col in df.columns if col.endswith('Proximity Density')]

 
        if selected_columns_v:
            unique_columns_v.add(selected_columns_v[0])
            selected_data_frames_v.append(df[selected_columns_v])
        else:
            print(f"No columns selected for file {csv_file}")

        if selected_columns_d:
            unique_columns_d.add(selected_columns_d[0])
            selected_data_frames_d.append(df[selected_columns_d])
        else:
            print(f"No columns selected for file {csv_file}")

    return selected_data_frames_v, selected_data_frames_d

def main(args):
    adata = read_h5ad(args.anndata)

    cells = list(set(adata.obs['CT_final']))
    adata = calculate_spatial_pscore(adata, cells)

    output_dir = "output_pscore"
    os.makedirs(output_dir, exist_ok=True)

    adata.write(filename=os.path.join(output_dir, args.output_andata_name))

    save_spatial_scores(adata, cells, output_dir)

    direct = os.getcwd()
    direct = os.path.dirname(__file__)
    direct = os.path.join(os.path.dirname(__file__), output_dir)
    
    file_list = os.listdir(direct)
    

    selected_data_frames_v, selected_data_frames_d = process_csv_files(file_list, direct)

    merged_df_volume = pd.concat(selected_data_frames_v, axis=1)
    merged_df_density = pd.concat(selected_data_frames_d, axis=1)
    
    concatenated_df_v= pd.merge(merged_df_volume, 
                           adata.obs[['Image','Group']].drop_duplicates(subset=['Image','Group']).set_index('Image'),
                           on='Image', how='inner')

    concatenated_df_d= pd.merge(merged_df_density, 
                           adata.obs[['Image','Group']].drop_duplicates(subset=['Image','Group']).set_index('Image'),
                           on='Image', how='inner')

    concatenated_df_v.to_csv(os.path.join(output_dir, args.output_volume_name))
    concatenated_df_d.to_csv(os.path.join(output_dir, args.output_density_name))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process spatial data.")
    parser.add_argument("--anndata", type=str, help="Path to the Anndata file.")
    parser.add_argument("--output_andata_name", type=str, default="output_adata.h5ad", help="Output file name for processed Anndata file.")
    parser.add_argument("--output_volume_name", type=str, default="merged_volume.csv", help="Output file name for merged volume CSV.")
    parser.add_argument("--output_density_name", type=str, default="merged_density.csv", help="Output file name for merged density CSV.")
    args = parser.parse_args()
    
    main(args)
