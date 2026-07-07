#' ############################################################################
#' File: 4_Clustering.R
#' Author: Celine I. Garcia Rodriguez
#' Date: 27/05/2026
#' 
#' Description:
#'  This script computates UMAps, nearest neighbors algorithm and clustering.
#' 
#' Libraries: Seurat, qs2
#' 
#' ############################################################################

#| 1. Packages
library(Seurat)
library(qs2)


SCZ_Seurat <- qs_read("/home/cgarcia/SCZ_Seurat_integrated.qs2")
options(future.globals.maxSize = 600 * 1024^3) # Increase RAM of R program

#| 2. UMAPs
SCZ_Seurat <- RunUMAP(SCZ_Seurat, 
                        reduction = "pca", 
                        dims = 1:20, 
                        reduction.name = "umap.unintegrated"
)

SCZ_Seurat <- RunUMAP(SCZ_Seurat, 
                        reduction = "integrated.rpca", 
                        dims = 1:20,
                        reduction.name = "umap.rpca"
)

SCZ_Seurat <- RunUMAP(SCZ_Seurat, 
                        reduction = "harmony_SI", 
                        dims = 1:20,
                        reduction.name = "umap.harmonySI"
)

SCZ_Seurat <- RunUMAP(SCZ_Seurat, 
                        reduction = "harmony", 
                        dims = 1:20,
                        reduction.name = "umap.harmony"
)

#| 3. Nearest Neighbors Algorithm
SCZ_Seurat <- FindNeighbors(SCZ_Seurat, 
                            reduction = "harmony", 
                            dims = 1:20, 
                            k.param = 30
)

#| 4. Clustering
SCZ_Seurat <- FindClusters(SCZ_Seurat, resolution = 0.5, verbose = FALSE)

#| 5. Fusion
SCZ_Seurat <- JoinLayers(SCZ_Seurat)

# Error:vector::reserve
