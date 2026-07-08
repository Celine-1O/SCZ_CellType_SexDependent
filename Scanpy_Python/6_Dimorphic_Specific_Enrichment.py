"""
File: 6_Dimorphic_Specific_Enrichment.py
Author: Celine I. Garcia Rodriguez

Description:
Sex-specific and dimorphic genes are selected from an AnnData object, specifica-
lly, Differential Gene Expression (DGE) matrices. 
It performs a Gene Set Enrichment Analysis (GSEA) across all cell types, and 
restricts Over Representation Analysis (ORA) exclusively to the computationally
specific and dimorphic gene subsets. 

Functions:
- Sex_Specific_Dimorphic
- Gsea_Analysis
- Ora_Analysis
"""

import anndata as ad
import pandas as pd
import numpy as np
import gseapy as gp
from gseapy import dotplot, barplot, enrichment_map
import os
import time
import networkx as nx
import matplotlib.pyplot as plt


PVAL_THR = 0.05
LFC_THR = 0.5
DATABASE = ['GO_Biological_Process_2023']

OUTPUT_ORA = "./ORA_results_SCZ/"
OUTPUT_GSEA = "./GSEA_results_SCZ/"

os.makedirs(OUTPUT_ORA, exist_ok=True)
os.makedirs(OUTPUT_GSEA, exist_ok=True)

print("\n=== PHASE 4: CLUSTERING & BENCHMARKING ===")

# 1. Sex-specific and sex dimorphic genes

def Sex_Specific_Dimorphic(df, layer_name, pval_thr, lfc_thr):
    """
    Extracts sex-specific and dimorphic gene signatures from a global DGE matrix 
    and exports the subset matrice.

    Parameters
    ----------
    df(pandas.DataFrame) : DGE matrix 
    layer_name (str) : hierarchical annotation 
    pval_thr (float) 
    lfc_thr (float)
        
    Returns
    -------
    tuple
        A tuple containing two pandas DataFrames: (df_dimorphic, df_specific).
    """
    print(f"\n--- Extracting Signatures: {layer_name.upper()} ---")
    
    df_male = df[df['Sex'] == 'Male'].copy()
    df_female = df[df['Sex'] == 'Female'].copy()
    
    df_merged = pd.merge(
        df_male, df_female, 
        on=['names', 'celltype', 'Layer'], 
        how='outer', suffixes=('_Male', '_Female')
    )

    sig_male = (df_merged['pvals_adj_Male'] < pval_thr) & (
        df_merged['logfoldchanges_Male'].abs() > lfc_thr)
    sig_female = (df_merged['pvals_adj_Female'] < pval_thr) & (
        df_merged['logfoldchanges_Female'].abs() > lfc_thr)

    # Dimorphic Genes
    sex_logfc = (df_merged['logfoldchanges_Male'] * df_merged['logfoldchanges_Female']) < 0
    df_dimorphic = df_merged[sig_male & sig_female & sex_logfc].copy()
    
    df_dimorphic['Upregulated_In'] = np.where(df_dimorphic['logfoldchanges_Male'] > 0, 'Male', 'Female')
    
    # Sex-Specific Genes
    not_sig_female = ((df_merged['pvals_adj_Female'] >= pval_thr) | (
        df_merged['pvals_adj_Female'].isna()))
    specific_male_mask = sig_male & not_sig_female
    
    not_sig_male = ((df_merged['pvals_adj_Male'] >= pval_thr) | (
        df_merged['pvals_adj_Male'].isna()))
    specific_female_mask = sig_female & not_sig_male

    df_merged.loc[specific_male_mask, 'Specific_To'] = 'Male'
    df_merged.loc[specific_female_mask, 'Specific_To'] = 'Female'
    
    df_specific = df_merged[specific_male_mask | specific_female_mask].copy()

    # Saving
    df_dimorphic.to_csv(f"Dimorphic_{layer_name}.csv", index=False)
    df_specific.to_csv(f"Specific_{layer_name}.csv", index=False)

    return df_dimorphic, df_specific

# 2. GSEA

def Gsea_Analysis(df, layer_name, database, output_dir):
    """
    Executes a Pre-ranked Gene Set Enrichment Analysis (GSEA) iteratively across 
    all discrete cell types within a specified hierarchical laye


    Parameters
    ----------
    df (pandas.DataFrame) : DGE matrix
    layel_name (str) : specific nomenclature identifier
    database : ['GO_Biological_Process_2023']).
    output_dir (str) : The root directory path 
    """
    print(f"\nInitiating GSEA for layer: {layer_name} ...")
    
    lfc_col = 'logfoldchanges' 
    pval_col = 'pvals'         
    gene_col = 'names'

    for celltype, df_cell in df.groupby('celltype'):
        
        # Sanitize data specifically for the current cell type
        df_clean = df_cell.dropna(subset=[lfc_col, pval_col, gene_col]).copy()
        df_clean = df_clean.drop_duplicates(subset=[gene_col])
        
        # Up to 10 genes
        if len(df_clean) < 9:
            continue
        
        # Algorithmic ranking metric
        df_clean['Rank'] = np.sign(df_clean[lfc_col]) * -np.log10(
            df_clean[pval_col] + 1e-300)
        rnk_df = df_clean[[gene_col, 'Rank']].sort_values(
            by='Rank', ascending=False)
        
        cell_str = str(celltype).replace("/", "_")
        cell_out_path = os.path.join(output_dir, layer_name, cell_str)
        os.makedirs(cell_out_path, exist_ok=True)
        
        try:
            pre_res = gp.prerank(
                rnk=rnk_df, gene_sets=database, 
                organism='human', outdir=cell_out_path,           
                permutation_num=1000, min_size=15, max_size=500, weight=1,
                ascending=False, threads=4, no_plot=True, seed=42, 
                verbose=False
            )
            
            res2d = pre_res.res2d.copy()
            
            if not res2d.empty:
                prefix = f"{layer_name}_{cell_str}"
                
                # 1. Dotplot 
                dotplot_file = os.path.join(cell_out_path, f"{prefix}_dotplot.png")
                dotplot(res2d, column="FDR q-val", title=f'GSEA: {prefix}',
                        cmap=plt.cm.viridis, size=5, figsize=(6, 7), 
                        cutoff=0.05, ofname=dotplot_file)
                plt.close('all') 
                
                # 2. Split Barplot
                sig_res = res2d[res2d['FDR q-val'] < 0.05].copy()
                if len(sig_res) >= 2:
                    sig_res = sig_res.sort_values(by='NES', ascending=False)
                    top_combined = pd.concat([sig_res.head(10), sig_res.tail(10)])
                    top_combined['Direction'] = np.where(top_combined['NES'] > 0, 
                                                         'Upregulated', 
                                                         'Downregulated')
                    
                    barplot_file = os.path.join(cell_out_path, f"{prefix}_barplot.png")
                    barplot(top_combined, column="FDR q-val", group='Direction', 
                            title=f'Top Pathways: {prefix}', figsize=(6, 10), 
                            color={'Upregulated': 'salmon', 'Downregulated': 'darkblue'}, ofname=barplot_file)
                    plt.close('all')

                # 3. Enrichment Network Map
                if len(sig_res) > 2:
                    nodes, edges = enrichment_map(res2d, cutoff=0.05)
                    if not edges.empty:
                        G = nx.from_pandas_edgelist(edges, source='src_idx', 
                                                    target='targ_idx', 
                                                    edge_attr=['jaccard_coef'])
                        fig, ax = plt.subplots(figsize=(8, 8))
                        pos = nx.layout.spring_layout(G) 
                        nx.draw_networkx_nodes(G, pos=pos, cmap=plt.cm.RdYlBu, 
                                               node_color=list(nodes.NES), 
                                               node_size=list(nodes.Hits_ratio * 1000))
                        nx.draw_networkx_labels(G, pos=pos, labels=nodes.Term.to_dict(), 
                                                font_size=8)
                        edge_weight = nx.get_edge_attributes(G, 'jaccard_coef').values()
                        nx.draw_networkx_edges(G, pos=pos, width=list(map(lambda x: x*10, edge_weight)), edge_color='#CDDBD4')
                        plt.title(f"Enrichment Map: {prefix}")
                        plt.savefig(os.path.join(cell_out_path, f"{prefix}_enrichmap.png"), bbox_inches='tight', dpi=300)
                        plt.close('all')
            
            print(f"    -> GSEA exported to: {cell_out_path}")
            
        except Exception as e:
            print(f"    -> Critical Error processing GSEA for {layer_name} - {celltype}. Exception: {e}")

# 3. ORA
def ORA_Analysis(gene_list, label, database, output_dir, max_retries=3):
    """
    Executes Over Representation Analysis (ORA) from sex-specific and sex-
    dimorphic genes.

    Parameters
    ----------
    gene_list (List) : significant gene identifiers 
    label (str) : specific nomenclature identifier used for output file naming 
    database : "GO_biological_2023")
    output_dir (str): Root directory path for exporting ORA results.
    max_retries (int, optional): Maximum number of API connection attempts
    
    """
    gene_list = pd.Series(gene_list).dropna().unique().tolist()
    if len(gene_list) < 9:
        return 
        
    out_path = os.path.join(output_dir, label)
    os.makedirs(out_path, exist_ok=True)
    
    enr = None
    
    for attempt in range(max_retries):
        try:
            time.sleep(3) # API rate limit 
            enr = gp.enrichr(gene_list=gene_list, gene_sets=database, 
                             organism='human', outdir=out_path, 
                             cutoff=0.05, no_plot=True)
            break 
            
        except Exception as e:
            print(f"    ! Critical API Error for {label}: {str(e)}")
    
    if enr is None or not hasattr(enr, 'res2d'):
        print(f"    ! Failed to complete ORA for {label} after {max_retries} attempts.")
        return
        
    res2d = enr.res2d[enr.res2d['Adjusted P-value'] < 0.05]
    
    if not res2d.empty:
        dotplot(res2d, column="Adjusted P-value", title=f"ORA: {label}", 
                figsize=(6, 7), ofname=os.path.join(out_path, "dotplot.png"))
        
        barplot(res2d.sort_values('Adjusted P-value').head(10), 
                column="Adjusted P-value", title=f"Top 10: {label}", 
                figsize=(6, 8), color='salmon', ofname=os.path.join(out_path, 
                                                                    "barplot.png"))
        
        if len(res2d) > 2:
            try:
                nodes, edges = enrichment_map(res2d, cutoff=0.05)
                if not edges.empty:
                    G = nx.from_pandas_edgelist(edges, source='src_idx', 
                                                target='targ_idx', 
                                                edge_attr=['jaccard_coef'])
                    fig, ax = plt.subplots(figsize=(8, 8))
                    pos = nx.layout.spring_layout(G)
                    nx.draw_networkx_nodes(G, pos=pos, node_size=500, 
                                           node_color='skyblue')
                    nx.draw_networkx_labels(G, pos=pos, 
                                            labels=nodes.Term.to_dict(), 
                                            font_size=8)
                    nx.draw_networkx_edges(G, pos=pos, edge_color='gray', alpha=0.5)
                    plt.title(f"Map: {label}")
                    plt.savefig(os.path.join(out_path, "enrichmap.png"), dpi=300)
                    plt.close('all')
            except Exception as e:
                print(f"    - Warning: Could not generate Enrichment Map for {label}: {e}")
                
        print(f"    -> {label}: ORA successfully completed.")


# 4. Workflow
if __name__ == "__main__":
    
    print("Loading CHECKPOINT 5...")
    adata = ad.read_h5ad("CHECKPOINT_5_DGE.h5ad")

    deg_keys = [key for key in adata.uns.keys() if key.startswith("DEG_")]

    for key in deg_keys:
        layer_name = key.replace("DEG_", "")
        print(f"\n Processing Global Layer: {layer_name.upper()}")
        
        # Isolate layer dataframe
        df_full = adata.uns[key].copy()

        # 1. Sex-specific and sex-dimorphic genes
        df_dimorphic, df_specific = Sex_Specific_Dimorphic(df_full, layer_name, 
                                                           PVAL_THR, LFC_THR, 
                                                           )
        
        # 2. GSEA Analysis 
        Gsea_Analysis(df_full, layer_name, DATABASE, OUTPUT_GSEA)
        
        
        
        # 3. ORA Analysis 
        print(f"\nInitiating ORA for isolated signatures in layer: {layer_name} ...")
        
        for data_matrix in [df_dimorphic, df_specific]:

            if 'Specific_To' in data_matrix.columns:
                group_col = 'Specific_To'
                paradigm = 'Specific'
                
            elif 'Upregulated_In' in data_matrix.columns:
                group_col = 'Upregulated_In'
                paradigm = 'Dimorphic_UpIn'
            else:
                continue

            for celltype, df_cell in data_matrix.groupby('celltype'):
                for sex_group, df_subset in df_cell.groupby(group_col):
                    
                    if len(df_subset) > 10:
                        label = f"{layer_name}_{str(celltype).replace('/', '_')}_{paradigm}_{sex_group}"
                        ORA_Analysis(df_subset['names'], label, DATABASE, OUTPUT_ORA)

    del adata
    print("Phase 6 Completed")