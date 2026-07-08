"""
File: 4_Clustering_Benchmark.py
Author: Celine I. Garcia Rodriguez
Description:

This script applies nearest neighbors algorithm, Leiden clustering and UMAP
visualization for every layers. Then, evaluates integration performance via 
benchmark. 

"""

import anndata as ad
import scanpy as sc
from scib_metrics.benchmark import Benchmarker

print("\n=== PHASE 4: CLUSTERING & BENCHMARKING ===")

print("\nLoading Checkpoint 3...")
adata = ad.read_h5ad("CHECKPOINT_3_Integrated.h5ad")

def NearestNeighbors_Clustering_UMAP(anndata_obj, key_label):
    """
    This function computes Nearest Neighbors Algorithm, Leiden clustering and
    UMAP visualization.
    Input:
        anndata(anndata_obj): Data matrix
        key_label = Layer where method was saved.
    Output:
        Saves a picture of UMAP
        
    """
    print(f"\nComputing nearest neighbors using: {key_label}")
    sc.pp.neighbors(anndata_obj,
                    use_rep=key_label, 
                    key_added=f"neighbors_{key_label}",
		            n_neighbors=50)
    
    print(f"\nClustering {key_label}")
    sc.tl.leiden(anndata_obj, 
                flavor="igraph", 
                key_added=f"leiden_{key_label}",
                neighbors_key=f"neighbors_{key_label}",
                resolution=0.5)
    
    print(f"\nComputing UMAP for {key_label}")
    sc.tl.umap(anndata_obj, neighbors_key=f"neighbors_{key_label}")
    anndata_obj.obsm[f"_UMAP_{key_label}"] = anndata_obj.obsm["X_umap"].copy()

    sc.pl.umap(anndata_obj, 
               color=["Database", f"leiden_{key_label}"], 
               save=f"_{key_label}.png", 
               show=False)
    
    return anndata_obj

# 1. PCA
adata = NearestNeighbors_Clustering_UMAP(adata, "X_pca")
# 2. Harmony
adata = NearestNeighbors_Clustering_UMAP(adata, "X_harmony")
# 3. scVI
adata = NearestNeighbors_Clustering_UMAP(adata, "X_scVI")

# 4. Saving Checkpoint
print("\nSaving Checkpoint 4 ...")
adata.write_h5ad("CHECKPOINT_4_Clusters.h5ad")

# 4. Benchmark
print("\n=== BENCHMARK ANALYSIS ===")
# As the datasets is enormous, 100k cells were subsampled
print("\nSubsetting 100,000 cells for benchmark...")
adata_bench = sc.pp.subsample(adata, n_obs=100000, copy=True)

print("Starting benchmark...")
bm = Benchmarker(
    adata_bench, 
    batch_key="Database",
    label_key="leiden_X_pca",
    embedding_obsm_keys=["X_pca", "X_scVI", "X_harmony"],
    n_jobs=-1
    )

bm.benchmark()

# Saving results
bm.plot_results_table(min_max_scale=False, save_dir=".")
df_benchmark = bm.get_results(min_max_scale=False)
print("\nBenchmark summary:")
print(df_benchmark)
df_benchmark.to_csv("Benchmark_Results_1.csv")

print("\nPhase 4 Completed")