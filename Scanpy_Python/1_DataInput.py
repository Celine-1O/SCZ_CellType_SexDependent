"""
File: Data input.py
Author: Celine I. Garcia Rodriguez
Description:

This script allows the input of 3 different single-cell datasets:
- BrainScope:<https://brainscope.gersteinlab.org/data/snrna_expr_matrices_zip/SZBDMulti-Seq.zip> 
- Zenodo: <https://zenodo.org/records/6921620>
- GEO: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254569>

and the creation of a unique Anndata.

"""


import anndata as ad
import pandas as pd
import gc
import scipy.sparse as sp
import numpy as np

print("\n=== PHASE 1: DATA LOADING AND FORMATTING ===")

# 1. ZENODO input
print("\nProcessing Zenodo...")
adata_zenodo = ad.read_h5ad("zenodo.h5ad")
adata_zenodo.obs['Database'] = 'Zenodo'
ids_zenodo = adata_zenodo.obs_names.str.split('_').str[0]

zenodo_dict = {
    "MB6": ["Male", "Schizophrenia"],
    "MB7": ["Female", "Control"],
    "MB8": ["Male", "Schizophrenia"], 
    "MB8-2": ["Male", "Schizophrenia"],
    "MB9": ["Female", "Control"],
    "MB10": ["Female", "Schizophrenia"],
    "MB11": ["Male", "Control"],
    "MB12": ["Female", "Schizophrenia"],
    "MB13": ["Male", "Control"],
    "MB14": ["Female", "Schizophrenia"],
    "MB15": ["Female", "Control"],
    "MB16": ["Male", "Control"],
    "MB17": ["Male", "Control"],
    "MB18-2": ["Male", "Control"],
    "MB19": ["Male", "Control"],
    "MB21": ["Female", "Control"],
    "MB22": ["Female", "Schizophrenia"],
    "MB23": ["Male", "Schizophrenia"],
    "MB50": ["Female", "Schizophrenia"], 
    "MB51": ["Male", "Control"],
    "MB52": ["Female", "Schizophrenia"], 
    "MB53": ["Female", "Control"],
    "MB54": ["Female", "Schizophrenia"],
    "MB55": ["Male", "Control"],
    "MB56": ["Male", "Schizophrenia"],
    "MB57": ["Female", "Control"]
}

adata_zenodo.obs['ID'] = ids_zenodo
adata_zenodo.obs['Sex'] = ids_zenodo.map(
    lambda x: zenodo_dict.get(x, ["Unknown", "Unknown"])[0])
adata_zenodo.obs['Disease'] = ids_zenodo.map(
    lambda x: zenodo_dict.get(x, ["Unknown", "Unknown"])[1])
adata_zenodo.obs = adata_zenodo.obs[['ID', 'Sex', 'Disease', 'Database']]


# 2. GEO input
print("\nProcessing GEO...")
adata_geo = ad.read_h5ad("GSE254569_adata_RNA.h5ad")

# Extracting raw counts from GEO Anndata
adata_geo.X = adata_geo.layers['counts'].copy()

# Deleting processed layers and embeddings to free up RAM
for key in list(adata_geo.layers.keys()): del adata_geo.layers[key]
for key in list(adata_geo.obsm.keys()): del adata_geo.obsm[key]

adata_geo.obs['Database'] = 'GEO'
adata_geo.obs['ID'] = adata_geo.obs['Donor']
adata_geo.obs.rename(columns={'Classification': 'Disease'}, inplace=True)
adata_geo.obs = adata_geo.obs[['ID', 'Sex', 'Disease', 'Database']]

# Filtering diseases
diseases = ['Control', 'Schizophrenia']
adata_geo = adata_geo[adata_geo.obs['Disease'].isin(diseases)].copy()
adata_geo.obs['Disease'] = (
    adata_geo.obs['Disease'].cat.remove_unused_categories()
    )

# 3. BRAINSCOPE input
print("\nProcessing BrainScope...")
adata_bs = ad.read_h5ad("Brainscope.h5ad")
ids_bs = adata_bs.obs_names.str.split('_').str[0]
adata_bs.obs['Database'] = 'BrainScope'
adata_bs.obs['ID'] = ids_bs
adata_bs.obs = adata_bs.obs[['ID', 'Sex', 'Disease', 'Database']]


# 4. Ensure Sparse format and Integers as anti-error system for scVI.
print("\nOptimizing matrices to sparse and integer...")
for adata in [adata_zenodo, adata_geo, adata_bs]:
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    adata.X.data = np.round(adata.X.data)


# 5. Fusion
print("\nMerging datasets...")
adata_scz = ad.concat([adata_zenodo, adata_geo, adata_bs], join='inner')

# Fill any null value
for col in adata_scz.obs.columns:
    adata_scz.obs[col] = (
        adata_scz.obs[col]
        .fillna("Unknown")
        .astype(str)
        .astype('category')
        )

# 6. Saving the new Anndata Object
print("\nClearing RAM...")
del adata_zenodo
del adata_geo
del adata_bs
gc.collect() 

print("\nSaving Checkpoint 1...")
adata_scz.write_h5ad("CHECKPOINT_1_Raw_.h5ad")
print("\nPhase 1 completed")
