#' ############################################################################
#' File: 3_IntegrationData.R
#' Author: Celine I. Garcia Rodriguez
#' Date: 27/05/2026
#' 
#' Description:
#'  This script integrates the datasets across the "Database" batch variable 
#' using RPCA and Harmony, and calculates benchmarking scores to evaluate the
#' best integration method.
#' 
#' Packages: Seurat, qs2, SeuratIntegrate, SeuratWrappers, harmony
#' 
#' ############################################################################

#| 1. Packages
library(Seurat)
library(qs2)
library(SeuratIntegrate)
library(SeuratWrappers)
library(harmony)

SCZ_Seurat <- qs_read("/home/cgarcia/SCZ_Seurat_preprocessed.qs2")

options(future.globals.maxSize = 600 * 1024^3) # Increase RAM of R program

#| 2. PCA
SCZ_Seurat <- RunPCA(SCZ_Seurat, 
                    features = VariableFeatures(SCZ_Seurat), 
                    npcs = 50,
                    verbose = FALSE
)

#| 3. IntegrateLayers

# ---------------------------------------------------------------------------
# RPCA (SeuratIntegrate)
# ---------------------------------------------------------------------------
SCZ_Seurat <- IntegrateLayers(
  object = SCZ_Seurat, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca", verbose = FALSE
)

# ---------------------------------------------------------------------------
# Harmony(SeuratIntegrate)
# ---------------------------------------------------------------------------
SCZ_Seurat <- IntegrateLayers(
  object = SCZ_Seurat, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony_SI", verbose = FALSE
)

# ---------------------------------------------------------------------------
# Harmony (native)
# ---------------------------------------------------------------------------
SCZ_Seurat <- RunHarmony(SCZ_Seurat, group.by.vars = "Database")

#| 4. Benchmark
batch.var <- "Database"

# ---------------------------------------------------------------------------
# RegressPC 
# ---------------------------------------------------------------------------
SCZ_Seurat <- AddScoreRegressPC(SCZ_Seurat, 
                                integration = "unintegrated", 
                                batch.var = batch.var, 
                                reduction = "pca"
)

SCZ_Seurat <- AddScoreRegressPC(SCZ_Seurat, 
                                integration = "rpca",
                                batch.var = batch.var,
                                reduction = "integrated.rpca"
)

SCZ_Seurat <- AddScoreRegressPC(SCZ_Seurat, 
                                integration = "harmonySI", 
                                batch.var = batch.var, 
                                reduction = "harmony_SI"
)

SCZ_Seurat <- AddScoreRegressPC(SCZ_Seurat, 
                                integration = "R_harmony", 
                                batch.var = batch.var, 
                                reduction = "harmony"
)


# ---------------------------------------------------------------------------
# DensityPC
# ---------------------------------------------------------------------------
SCZ_Seurat <- AddScoreDensityPC(SCZ_Seurat, 
                                integration = "unintegrated", 
                                batch.var = batch.var, reduction = "pca"
)

SCZ_Seurat <- AddScoreDensityPC(SCZ_Seurat, integration = "rpca", 
                                batch.var = batch.var, 
                                reduction = "integrated.rpca"
)

SCZ_Seurat <- AddScoreDensityPC(SCZ_Seurat, 
                                integration = "harmony", 
                                batch.var = batch.var, 
                                reduction = "harmony_SI"
)

SCZ_Seurat <- AddScoreDensityPC(SCZ_Seurat, 
                                integration = "R_harmony", 
                                batch.var = batch.var, 
                                reduction = "harmony"
)

# ---------------------------------------------------------------------------
# CellCycleScoring
# ---------------------------------------------------------------------------

SCZ_Seurat <- CellCycleScoringPerBatch(
  SCZ_Seurat, 
  batch.var = batch.var, 
  s.features = cc.genes$s.genes, 
  g2m.features = cc.genes$g2m.genes
)

SCZ_Seurat <- AddScoreRegressPC.CellCycle(SCZ_Seurat, 
                                            integration = "unintegrated", 
                                            batch.var = batch.var, 
                                            what = "pca", 
                                            compute.cc = FALSE, 
                                            dims.use = 1:20
)

SCZ_Seurat <- AddScoreRegressPC.CellCycle(SCZ_Seurat, 
                                            integration = "rpca", 
                                            batch.var = batch.var, 
                                            what = "integrated.rpca", 
                                            compute.cc = FALSE, 
                                            dims.use = 1:20
)

SCZ_Seurat <- AddScoreRegressPC.CellCycle(SCZ_Seurat, 
                                            integration = "harmony_SI", 
                                            batch.var = batch.var, 
                                            what = "harmony_SI", 
                                            compute.cc = FALSE, 
                                            dims.use = 1:20
)

SCZ_Seurat <- AddScoreRegressPC.CellCycle(SCZ_Seurat, 
                                            integration = "R_harmony", 
                                            batch.var = batch.var, 
                                            what = "harmony", 
                                            compute.cc = FALSE, 
                                            dims.use = 1:20
)

SCZ_Seurat <- ScaleScores(SCZ_Seurat)

#| Benchmark Summary

scores_seurat <- data.frame(IntegrationScores(SCZ_Seurat, scaled = TRUE))

PlotScores(SCZ_Seurat, rescale = FALSE)
PlotScores(SCZ_Seurat, plot.type = "lollipop")

#| 5. Saving file
qs_save(SCZ_Seurat, "/home/cgarcia/SCZ_Seurat_integrated.qs2")
