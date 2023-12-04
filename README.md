# From-Pixels-to-Insights

Readme

Step 1: This script performs QC and normalization. The input is CSV from Qupath, and the output is anndata, which is normalised and QCed.

Script_1_QC_Normalization.py --input experssion_fromQupath.csv --output anndata_qced_normalized.h5ad --exclude_core "core1,core2" --exclude_marker "marker1,marker2"



Step 2: Jupyter Notebook for creating heatmaps and manual cell typing using a Python dictionary.


Normalised anndata is headed in this notebook. And cell typing applied manually using python dictionary.


Step 3: Neighbourhood analysis The input is metadata with cell type columns (CT_final).

python Script_3_NB.py --input_file Data_Cell_Clinical.csv --output_file_test output_forTest.csv --output_file_metadata Data_Cell_Clinical_Group_NB.csv


Step 4: Cell enrichments in NB The input is metadata, which includes cell type and NB.

python Script_4_cellEnrich_InNB.py --input Data_Cell_Clinical_Group_NB.csv --output output_cellInrichmentInNB.csv


Step 5: Input is count and meta data, or anndata.

python Script_5_Dist_intraction_scimap.py --data_file data/ counts_table.csv --meta_file data/Data_Cell_Clinical_Group_NB.csv --distance_output_name my_distance_data.csv --interaction_output_name my_interaction_data.csv --anndata_output_name anndata_distance.h5ad

Step 6: The input file is the data from step 5.

python Script_6_Pscore_scimap.py --anndata output_dist_int/anndata_distance.h5ad --output_andata_name output_adata_pscore.h5ad --output_volume_name merged_volume.csv --output_density_name merged_density.csv



Step 7: Test the data frame, including features, and group (y label) in columns. Groups include responders and non-responders.

python Script_7_tTest.py --input Tumor.csv --output results.csv --group1 Responder --group2 Non_Responder --label Group
![image](https://github.com/clinicalomx/From-Pixels-to-Insights/assets/152827690/c059c9ea-0c16-4c84-a3cb-9f40a0835fa7)
