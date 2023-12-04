# usage > python Script_3_NB.py --input_file Data_Cell_Clinical.csv --output_file_test output_forTest.csv --output_file_metadata Data_Cell_Clinical_Group_NB.csv
import os
from os import listdir, path, getcwd
import glob
import re
from joblib import parallel_backend
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import seaborn as sns
import matplotlib
from scipy import stats
import argparse
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import phenograph
import copy
from scipy import sparse
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import NearestNeighbors
import time
import sys
from sklearn.cluster import MiniBatchKMeans
import seaborn as sns


warnings.filterwarnings('ignore')

###################################################

def parse_arguments():
    parser = argparse.ArgumentParser(description="Spatial Transcriptomics Analysis Script")
    parser.add_argument("--input_file", required=True, help="Path to the input CSV file (Data_Cell_Clinical_Group_headneck.csv)")
    parser.add_argument("--output_file_test", required=True, help="Path to the output CSV file for ML and survival analysis results")
    parser.add_argument("--output_file_metadata", required=True, help="Path to the output CSV file for metadata results")
    return parser.parse_args()


# Function for identifying the windows
def get_windows(job, n_neighbors):
    '''
    For each region and each individual cell in dataset, return the indices of the nearest neighbors.

    'job:  meta data containing the start time,index of region, region name, indices of region in original dataframe
    n_neighbors:  the number of neighbors to find for each cell
    '''
    start_time, idx, tissue_name, indices = job
    job_start = time.time()

    print("Starting:", str(idx + 1) + '/' + str(len(exps)), ': ' + exps[idx])

    # tissue_group: a grouped data frame with X and Y coordinates grouped by unique tissue regions
    tissue = tissue_group.get_group(tissue_name)

    to_fit = tissue.loc[indices][['x', 'y']].values

    # Unsupervised learner for implementing neighbor searches.
    fit = NearestNeighbors(n_neighbors=n_neighbors).fit(tissue[['x', 'y']].values)

    # Find the nearest neighbors

    m = fit.kneighbors(to_fit) 

    m = m[0], m[1]

    ## sort_neighbors
    args = m[0].argsort(axis=1)

    add = np.arange(m[1].shape[0]) * m[1].shape[1]

    sorted_indices = m[1].flatten()[args + add[:, None]]

    neighbors = tissue.index.values[sorted_indices]

    end_time = time.time()

    print("Finishing:", str(idx + 1) + "/" + str(len(exps)), ": " + exps[idx], end_time - job_start,
          end_time - start_time)
    return neighbors.astype(np.int32)



if __name__ == "__main__":
    # Parse command-line arguments
    args = parse_arguments()

    # Load data from the specified input file
    input_file_path = args.input_file
    data = pd.read_csv(input_file_path)

    # Create output directory if it doesn't exist
    output_directory = "output_NB"
    os.makedirs(output_directory, exist_ok=True)

####################################



input_file_path = args.input_file
data = pd.read_csv(input_file_path)

data=data.reset_index(drop=True)
data['cell_id'] = data.index

# make dummy variables
data2 = pd.concat([data,pd.get_dummies(data['CT_final'])], axis = 1)
# Extract the cell types with dummy variables
sum_cols = data2['CT_final'].unique()
values = data2[sum_cols].values




###########################################################################

# Keep the X and Y coordianates + the tissue regions >> then group by tissue regions 
tissue_group = data2[['x','y','Image']].groupby('Image')

# Create a list of unique tissue regions
exps = list(data2['Image'].unique())


tissue_chunks = [(time.time(),exps.index(t),t,a) for t,indices in tissue_group.groups.items() for a in np.array_split(indices,1)]

# Get the window (the 10 closest cells to each cell in each tissue region)
tissues = [get_windows(job,10) for job in tissue_chunks]

########################################
ks = [10]
out_dict = {}
for k in ks:
    for neighbors, job in zip(tissues, tissue_chunks):
        chunk = np.arange(len(neighbors))  # indices
        tissue_name = job[2]
        indices = job[3]
        window = values[neighbors[chunk, :k].flatten()].reshape(len(chunk), k, len(sum_cols)).sum(axis=1)
        out_dict[(tissue_name, k)] = (window.astype(np.float16), indices)



#########
        
keep_cols = ['x','y','Image','CT_final']
windows = {}
for k in ks:
    window = pd.concat(
        [pd.DataFrame(out_dict[(exp, k)][0], index=out_dict[(exp, k)][1].astype(int), columns=sum_cols) for exp in
         exps], 0)
    window = window.loc[data2.index.values]
    window = pd.concat([data2[keep_cols], window], 1)
    windows[k] = window


neighborhood_name = "neighborhood"+str(k)
k_centroids = {}

windows2 = windows[10]

##################################
km = MiniBatchKMeans(n_clusters = 10,random_state=0)

labelskm = km.fit_predict(windows2[sum_cols].values)
k_centroids[10] = km.cluster_centers_
data2['neighborhood10'] = labelskm
data2[neighborhood_name] = data2[neighborhood_name].astype('category')

cell_order=set(data2['CT_final'])

specific_item = 'Unclassified'
cell_order.discard(specific_item)

##############################################
sns.set(font_scale=1)

niche_clusters = (k_centroids[10])
tissue_avgs = values.mean(axis = 0)

fc = np.log2(((niche_clusters+tissue_avgs)/(niche_clusters+tissue_avgs).sum(axis = 1, keepdims = True))/tissue_avgs)
fc = pd.DataFrame(fc,columns = sum_cols)



s=sns.clustermap(fc.loc[[0,1,2,3,4,5,6,7,8,9],cell_order], vmin =-3,vmax = 3,cmap = 'coolwarm',row_cluster = False)






heatmap_figure_path = os.path.join(output_directory, "NB_heatmp.svg")
s.savefig(heatmap_figure_path, dpi=300, bbox_inches="tight")
print(f"Heatmap figure saved to {heatmap_figure_path}")



##########################################


fc = data2.groupby(['Image','Group']).apply(lambda x: x['neighborhood10'].value_counts(sort = False,normalize = True))

fc.columns = range(10)
melt = pd.melt(fc.reset_index(),id_vars = ['Image','Group'])
melt = melt.rename(columns = {'variable':'neighborhood','value':'frequency of neighborhood'})
melt['neighborhood'] = melt['neighborhood'].map(
    {
        0: 'NB 0',
        1: 'NB 1',
        2: 'NB 2',
        3: 'NB 3',
        4: 'NB 4',
        5: 'NB 5',
        6: 'NB 6',
        7: 'NB 7',
        8: 'NB 8',
        9: 'NB 9'
    }
)

f,ax = plt.subplots(figsize = (15,7))
#sns.stripplot(data = melt, hue = 'Group',dodge = True,alpha = 0.1,x ='neighborhood', y ='frequency of neighborhood')
sns.pointplot(data = melt,hue = 'Group',dodge = .1,join = False,x ='neighborhood', y ='frequency of neighborhood')



handles, labels = ax.get_legend_handles_labels()


plt.xticks(rotation=0, fontsize="15", ha="center")


ax.legend(handles[:2], labels[:2], title="Groups",
          handletextpad=0, columnspacing=2,
          loc="upper left", ncol=3, frameon=True)
plt.tight_layout()

plt.xticks(rotation = 45)


frequency_figure_path = os.path.join(output_directory, "NB_feq.svg")
plt.savefig(frequency_figure_path, dpi=300, bbox_inches="tight")
print(f"Frequency figure saved to {frequency_figure_path}")

##########

df=fc.reset_index()


output_file_test_path = os.path.join(output_directory, args.output_file_test)
df.to_csv(output_file_test_path, index=False)
print(f"table results saved to {output_file_test_path}")








df=pd.read_csv(input_file_path)

data3=data2[['Object ID', 'neighborhood10']]
merged_df = pd.merge(df, data3, on='Object ID', how='inner') 

output_file_meta_path = os.path.join(output_directory, args.output_file_metadata)

merged_df.to_csv(output_file_meta_path, index=False)
print(f"ML and survival analysis results saved to {output_file_meta_path}")





