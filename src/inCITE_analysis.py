# Daniel Zhu
# March 19th, 2021
# Exploring inCITE-seq data and the AnnData format.
import numpy as np
import scanpy as sc
from harmony import harmonize
import matplotlib.pyplot as plt
import os

# Specify the path to the Anndata object:
path = os.getcwd() + '\\data\\inCITE.h5ad'
adata = sc.read_h5ad(path)

# Visualize the completely unprocessed UMAP:
sc.pp.neighbors(adata)
sc.tl.umap(adata)
# Color code by batch, and then uncomment line 21 to view an assignment-color coded version (Kainic acid treatment or PBS for this dataset).
sc.pl.umap(adata, color='batch', title="Unprocessed, colored by batch", show=True, save='_unprocessed_batch.png')
#sc.pl.umap(adata, color='assignment', title="Unprocessed, colored by assignment (PBS or kainic acid)", show=True, save='unprocessed_assignment.png')

# Step 1: (Optional preprocessing) select highly variable genes:
# Step 1, part 1: log-scale the data:
adata_copy = adata.copy()
sc.pp.log1p(adata_copy)
# Step 1, part 2: select highly variable genes:
sc.pp.highly_variable_genes(adata_copy, min_mean=0.001, max_mean=0.15, min_disp=0.25)
# And then actually do the filtering:
adata_copy = adata_copy[:, adata_copy.var.highly_variable]
print("Selected highly-variable genes checkpoint: \n", adata_copy)

# Step 2: Further processing on highly-variable genes:
# Step 2, part 1: regress out UMI counts (transcript counts) and mitochondrial content (frac_mito):
adata_regr = sc.pp.regress_out(adata_copy, keys=['n_counts', 'frac_mito'], copy=True)
# Step 2, part 2: scale highly-variable genes:
sc.pp.scale(adata_regr, max_value=10)  # scales to unit variance, and removes genes with stdev higher than 10.
# Step 2, part 3: PCA (not sure if this is necessary, but the harmony demo used the PCA representation for the batch correction function):
sc.pp.pca(adata_regr)
print("Regressed out UMI and mito checkpoint: \n", adata_regr)


# Step 3: Harmony for batch correction as well as assignment correction.
# use_rep is an obsm entry; the actual correction is done with the code below and adds the corrected representation to obsm:
adata_regr.obsm['X_harmony'] = harmonize(np.array(adata_regr.obsm['X_pca']),
									adata_regr.obs,
									batch_key = ['batch','assignment'])
print("Batch correction checkpoint: \n", adata_regr)

# Visualize UMAP of the processed data (pre-batch correction):
sc.pp.neighbors(adata_regr, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata_regr)
sc.tl.leiden(adata_regr, resolution=1)   # cluser cells using the Leiden algorithm; use the default resolution for now (higher values lead to more clusters).
sc.pl.umap(adata_regr, color='batch', title='Processed, pre-batch correction, colored by batch', show=True, save='_processed_batch.png')
#sc.pl.umap(adata_regr, color='assignment', title='Processed, pre-batch correction, colored by assignment' show=True, save='processed_assignment.png')

# Visualize UMAP of the processed data (post-batch correction):
sc.pp.neighbors(adata_regr, use_rep='X_harmony', n_neighbors=10, n_pcs=30)
sc.tl.umap(adata_regr)
sc.tl.leiden(adata_regr, resolution=1)   # cluser cells using the Leiden algorithm; use the default resolution for now (higher values lead to more clusters).
sc.pl.umap(adata_regr, color='batch', title='Processed, post-batch correction, colored by batch', show=True, save='_harmony_batch.png')
#sc.pl.umap(adata_regr, color='assignment', title='Processed, post-batch correction, colored by assignment', show=True, save='harmony_assignment.png')