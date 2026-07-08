"""
File: 3_DataIntegration.py
Author: Celine I. Garcia Rodriguez
Description:

This script performs:
- Dimensionality Reduction (Principal Components Analysis) 
- 2 different Integration Methods:
    - Harmony
    - scVI
"""

import anndata as ad
import scanpy as sc
import scvi
import harmonypy as hm 


print("\n=== PHASE 3: DATA INTEGRATION ===")

print("\nLoading Checkpoint 2...")
adata = ad.read_h5ad("CHECKPOINT_2_Filtered.h5ad")

# 1. PCA
print("\nPerforming Principal Components Analysis...")
sc.pp.scale(adata, max_value=10, zero_center=False)
sc.pp.pca(adata, svd_solver='covariance_eigh', n_comps=50, zero_center=True)

# 2. Harmony
print("\nPerforming Harmony Integration...")
pcs = adata.obsm['X_pca'] # (n_cells, n_pcs)
harmony_out = hm.run_harmony(pcs, adata.obs, "Database")
adata.obsm['X_harmony'] = harmony_out.Z_corr

# 3. scVI model
scvi.settings.dl_num_workers = 30 
print("\nSetting scVI model...")
# layer = "counts" because is where raw datas (unnormalization ) is.
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="Database")

model = scvi.model.SCVI(adata)

print("\nTraining scVI...")
model.train(max_epochs=100, 
            early_stopping=True, 
            early_stopping_patience=10, 
            train_size=0.9, 
            batch_size=512) 
 
print("\nExtracting latent space...")
adata.obsm["X_scVI"] = model.get_latent_representation()
model.save("scvi_model_scz/", overwrite=True)

# 4. Saving Checkpoint
print("\nSaving Checkpoint 3...")
adata.write_h5ad("CHECKPOINT_3_Integrated.h5ad")
print("\nPhase 3 completed")

