"""
File: 5_CellTypeAnnotation_DEGs.py
Author: Celine I. Garcia Rodriguez
Date: 25/05/2026
Description:

This script extracts the marker genes of each Leiden Harmony cluster and saves
them in a CSV file. After identify each cluster using CellKB, cellular annota-
tion and stratified differential gene expression(DGE) and differential 
abundance (DA) analysis between cell-types and sex were performed.

"""

import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import glob
import os

print("\n=== PHASE 5: CELL-TYPE ANNOTATION & DEG ANALYSIS ===")

print("\nLoading Checkpoint 4...")
adata = ad.read_h5ad("CHECKPOINT_4_Clusters.h5ad")

# 1m. Searching marker genes.
print("Computing marker genes from Harmony Clusters...")
sc.tl.rank_genes_groups(adata, groupby="leiden_X_harmony", method="wilcoxon", 
                        layer='log_normalization',use_raw=False, pts=True)

print("Extracting results into a CSV...")
df_markers = sc.get.rank_genes_groups_df(adata, group=None)
df_markers.to_csv("Harmony_markers.csv", index=False)
print("CSV Dataframe saved")

# 2. Cell-type Annotation at different depth levels.

adata.obs["celltype_deep"] = adata.obs["leiden_X_harmony"].map(
    {
        "0": "L3_4_5",
        "1": "LAMP5",
        "2": "L6",
        "3": "L2_3",
        "4": "VIP",
        "5": "Basket",
        "6": "L2",
        "7": "L4_5",
        "8": "L5_6",
        "9": "SST",
        "10": "L4_6",
        "11": "L2_3",
        "12": "L4_6",
        "13": "Reelin",
        "14": "L5_6",
        "15": "LAMP5",
        "16": "L5_6",
        "17": "PVALB_Chandelier",
        "18": "L5_6",
        "19": "Oligodendrocyte",
        "20": "L2_3",
        "21": "Astrocyte",
        "22": "Endothelial",
        "23": "Microglia",
        "24": "OPC"
    }
)


adata.obs["celltype_medium"] = adata.obs["leiden_X_harmony"].map(
    {
        "0": "Excitatory_deep",
        "1": "LAMP5",
        "2": "Excitatory_deep",
        "3": "Excitatory_superficial",
        "4": "VIP",
        "5": "Basket",
        "6": "Excitatory_superficial",
        "7": "Excitatory_deep",
        "8": "Excitatory_deep",
        "9": "SST",
        "10": "Excitatory_deep",
        "11": "Excitatory_superficial",
        "12": "Excitatory_deep",
        "13": "Reelin",
        "14": "Excitatory_deep",
        "15": "LAMP5",
        "16": "Excitatory_deep",
        "17": "PVALB_Chandelier",
        "18": "Excitatory_deep",
        "19": "Oligodendrocyte",
        "20": "Excitatory_superficial",
        "21": "Astrocyte",
        "22": "Endothelial",
        "23": "Microglia",
        "24": "OPC"
    }
)

adata.obs["celltype_superficial"] = adata.obs["leiden_X_harmony"].map(
    {
        "0": "Excitatory_neuron",
        "1": "Inhibitory_neuron",
        "2": "Excitatory_neuron",
        "3": "Excitatory_neuron",
        "4": "Inhibitory_neuron",
        "5": "Inhibitory_neuron",
        "6": "Excitatory_neuron",
        "7": "Excitatory_neuron",
        "8": "Excitatory_neuron",
        "9": "Inhibitory_neuron",
        "10": "Excitatory_neuron",
        "11": "Excitatory_neuron",
        "12": "Excitatory_neuron",
        "13": "Inhibitory_neuron",
        "14": "Excitatory_neuron",
        "15": "Inhibitory_neuron",
        "16": "Excitatory_neuron",
        "17": "Inhibitory_neuron",
        "18": "Excitatory_neuron",
        "19": "Oligodendrocyte",
        "20": "Excitatory_neuron",
        "21": "Astrocyte",
        "22": "Endothelial",
        "23": "Microglia",
        "24": "OPC"
    }
)



# 3. Differential Gene Expression Analysis per Sex and Condition
output_dir = "DGE_results_SCZ"
os.makedirs(output_dir, exist_ok=True)

depth_layers = ['celltype_superficial','celltype_medium','celltype_deep'] 
genders = ["Male", "Female"]


print("\nStarting DGEs Analysis...")

for layer in depth_layers:        
    print(f"\nProcessing at {layer}")
    cell_types = adata.obs[layer].unique()
    layer_deg_list = []
    
    for ct in cell_types:
        obj_ct = adata[adata.obs[layer] == ct].copy()
        
        for gender in genders:
            gender_ct = obj_ct[obj_ct.obs['Sex'] == gender].copy()
            cases = gender_ct.obs['Disease'].unique()

            print(f"Computing: {ct} | Sex: {gender}")
                
            sc.tl.rank_genes_groups(
                gender_ct, 
                groupby='Disease', 
                groups=['Schizophrenia'], 
                reference='Control', 
                method='wilcoxon',
                use_raw=False
                )
            
            df_complete = sc.get.rank_genes_groups_df(gender_ct, 
                                                      group='Schizophrenia')
            
            df_complete['celltype'] = str(ct)
            df_complete['Sex'] = str(gender)
            df_complete['Layer'] = layer
            
            layer_deg_list.append(df_complete)

            # Saving individual CSV
            ct_name = str(ct).replace(' ', '_').replace('/', '_')
            filename = f"DEG_{layer}_{ct_name}_{gender}.csv"            
            df_complete.to_csv(filename, index=False)

                                
    if layer_deg_list:
        df_general_layer= pd.concat(layer_deg_list, ignore_index=True)

        # storage in adata uns
        adata.uns[f"DEG_{layer}"] = df_general_layer

        # Layer specific csv
        df_general_layer.to_csv(f"_DEG_{layer}.csv", index=False)

        print(f" adata.uns['DEG_{layer}'] and CSV file")
    else:
        print(f"No {layer}.")


# 4. Saving Checkpoint
print("\nSaving Checkpoint 5...")
adata.write_h5ad("CHECKPOINT_5_DGE.h5ad")


# 5. Differential Abundance Analysis
print("\nStarting DAs Analysis...")

metadata = adata.obs[['ID', 'Disease', 'Sex']].drop_duplicates().set_index('ID')
comparisons = [
    ('Sex', 'Female', 'Disease', 'Schizophrenia', 
     'Control', 'Schizophrenia_vs_Control_in_Females'),
    ('Sex', 'Male', 'Disease', 'Schizophrenia', 
     'Control', 'Schizophrenia_vs_Control_in_Males'),
    ('Disease', 'Control', 'Sex', 'Male', 
     'Female', 'Male_vs_Female_in_Control'),
    ('Disease', 'Schizophrenia', 'Sex', 'Male', 
     'Female', 'Male_vs_Female_in_Schizophrenia')
]

for layer in depth_layers:        
    print(f"\nEvaluating compositional variance at: {layer}")

    # Formulate empirical compositional matrix
    counts_per_donor = pd.crosstab(adata.obs['ID'], adata.obs[layer])
    proportions_per_donor = counts_per_donor.div(counts_per_donor.sum(axis=1), 
                                                axis=0) * 100
    proportions_per_donor = proportions_per_donor.join(metadata)
        
    cell_types = counts_per_donor.columns.tolist()    
    layer_da_list = []
        
    # Execute non-parametric assessments 
    for filter_col, filter_val, test_col, group_a, group_b, comp_name in comparisons:
        subset_df = proportions_per_donor[proportions_per_donor[filter_col] == filter_val]
        
        for cell in cell_types:
            props_group_a = subset_df[subset_df[test_col] == group_a][cell]
            props_group_b = subset_df[subset_df[test_col] == group_b][cell]                 
                
            stat, p_val = mannwhitneyu(props_group_a, props_group_b, 
                                       alternative='two-sided')
            mean_a = props_group_a.mean()
            mean_b = props_group_b.mean()
            
            log2fc = np.log2((mean_a + 1e-5) / (mean_b + 1e-5))
                
            layer_da_list.append({
                'Annotation_Level': layer,
                'Comparison': comp_name,
                'Cell_Type': cell,
                f'Mean_{group_a}': mean_a,
                f'Mean_{group_b}': mean_b,
                'Log2FC': log2fc,
                'p_value': p_val
            })

    results_df = pd.DataFrame(layer_da_list)
    
    # Benjamini-Hochberg False Discovery Rate correction
    reject, p_adj, _, _ = multipletests(results_df['p_value'], alpha=0.05, method='fdr_bh')
    results_df['p_adj'] = p_adj
    results_df['Significant'] = reject
    
    # Sorting 
    results_df = results_df.sort_values(by=['Comparison', 'p_adj'])
    
    # Export comprehensive statistical outcomes
    csv_filename = os.path.join(output_dir, f"DA_comprehensive_results_{layer}.csv")
    results_df.to_csv(csv_filename, index=False)
    
    # Isolate and report statistically robust findings
    significant_df = results_df[results_df['Significant'] == True]
    print(f"Identified {len(significant_df)} significant DA shifts in {layer}.")


print("Phase 5 Completed")
