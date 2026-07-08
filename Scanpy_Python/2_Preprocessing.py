"""
File: 2_Preprocessing.py
Author: Celine I. Garcia Rodriguez
Description:

This script performs Quality Control, Normalization and identifies Highly 
Variable Genes (HVGs) of the AnnData object using Scanpy framework:
<https://scanpy.scverse.org/en/stable/tutorials/basics/clustering.html>


"""

import anndata as ad
import scanpy as sc

print("\n=== Phase 2: PREPROCESSING ===")

# 1. Loading Anndata Object
print("\nLoading Checkpoint 1...")
adata = ad.read_h5ad("CHECKPOINT_1_Raw_.h5ad")


# 2. Quality Control
print("\nFiltering mitochondrial genes...")
# mitochondrial genes, "MT-" for human
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, log1p=True)
adata = adata[adata.obs['pct_counts_mt'] < 5, :].copy()

print("\nFiltering cells with less than 200 genes expressed...")
sc.pp.filter_cells(adata, min_genes=200)

# 3. Saving count data
adata.layers["counts"] = adata.X.copy() # Save raw data
adata.raw = adata # Keep all genes

# 4. Doublet detection
print("\nDetecting and filtering doublets...")
sc.pp.scrublet(adata, batch_key="ID")
adata = adata[adata.obs['predicted_doublet'] == False].copy()

# 5. Normalization
print("\nNormalizing to media total counts and logarithmize the data...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["log_normalization"] = adata.X.copy()


# 6. Feature selection (method vst)
print("\nSelecting highly variable genes (HVGs)...")
sc.pp.highly_variable_genes(
    adata, 
    n_top_genes=2000, 
    subset=True, 
    layer="counts",              
    flavor="seurat_v3_paper",    
    batch_key="Database"         
)

# 7. Saving the processed Anndata Object
print("\nSaving Checkpoint 2 ...")
adata.write_h5ad("CHECKPOINT_2_Filtered.h5ad")
print("\nPhase 2 Completed")
