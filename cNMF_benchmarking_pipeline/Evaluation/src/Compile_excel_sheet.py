import pandas as pd
import muon as mu
from scipy import sparse
import sys
import os
import numpy as np
from scipy import stats
from mygene import MyGeneInfo
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import entropy



# Change path to wherever you have repo locally
sys.path.append('/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline')


from Plotting.src import rename_adata_gene_dictionary ,compute_gene_waterfall_cor


# Compile simple sheets 

#-------------- helper methods for loading sheets--------------

def compile_Program_loading_score_sheet_long(mdata, num_gene = 300):

    print('Load program loadings data in long form')

    program_loading_df = pd.DataFrame(data=mdata['cNMF'].varm["loadings"], columns = mdata['rna'].var_names)
    num_gene = 300

    top_df = program_loading_df.apply(
    lambda row: row.nlargest(num_gene).index.tolist(),
    axis=1
    )

    result_df = pd.DataFrame(top_df.tolist(), columns=range(1, num_gene+1))
    result_df.index = program_loading_df.index

    result_df.index.name = "Program"

    mg = MyGeneInfo()
    long_data = []

    for program_idx, row in result_df.iterrows():

        # Query all genes for this program at once
        genes_list = row.dropna().astype(str).str.strip().tolist()
        genes_list = [g for g in genes_list if g]  # Remove empty strings
        
        if not genes_list:
            continue
        
        # Get results for all genes in this program
        results = mg.querymany(genes_list, scopes='symbol', fields='summary', species='human', verbose=False)
        
        # Create a lookup dictionary for faster access
        summary_dict = {}
        for r in results:
            q = r.get('query')
            if q not in summary_dict:  # keep first (best) match
                summary_dict[q] = r.get('summary', 'N/A')
        
        # Now iterate through the ranked genes
        for rank, gene in enumerate(row, 1):
            
            # Skip NaN values
            if pd.notna(gene) and gene != '':
                gene_clean = str(gene).strip()
                annotation = summary_dict.get(gene_clean, 'N/A')
                
                long_data.append({
                    'Program': program_idx,
                    'Rank': rank,
                    'Gene': gene_clean,
                    'Annotation': annotation
                })

    # Convert to DataFrame
    annotation_df = pd.DataFrame(long_data)

    return annotation_df

def compile_Program_loading_score_sheet_flat(mdata, num_gene = 300):

    print('Load program loadings data in flat form')


    program_loading_df = pd.DataFrame(data=mdata['cNMF'].varm["loadings"], columns = mdata['rna'].var_names)

    top_df = program_loading_df.apply(
    lambda row: row.nlargest(num_gene).index.tolist(),
    axis=1
    )

    result_df = pd.DataFrame(top_df.tolist(), columns=range(1, num_gene+1))
    result_df.index = program_loading_df.index

    result_df.index.name = "Program"

    return result_df

def Compile_GO_sheet(GO_path, gene_num = 5):

    print('Load GO data')

    df = pd.read_csv(GO_path, sep = "\t", index_col = 0)
    df = df.reset_index().set_index("Term")

    df["Genes"] = df["Genes"].str.split(';').str[:gene_num].str.join(";")

    return df

def Compile_Geneset_sheet(Geneset_path, gene_num = 5):

    print('Load geneset data')


    df = pd.read_csv(Geneset_path, sep = "\t", index_col = 0)
    df = df.reset_index().set_index("Term")

    df["Genes"] = df["Genes"].str.split(';').str[:gene_num].str.join(";")

    return df

def Compile_Trait_sheet(Trait_path, gene_num = 5):

    print('Load trait data')

    df = pd.read_csv(Trait_path, sep = "\t", index_col = 0)
    df = df.reset_index().set_index("Term")

    df["Genes"] = df["Genes"].str.split(';').str[:gene_num].str.join(";")

    return df

def Compile_Perturbation_sheet(Perturbation_path, Sample = ["D0", "sample_D1","sample_D2","sample_D3" ]):

    print('Load perturbation data')


    combined_conditions = []
    for samp in Sample: 
        df = pd.read_csv(f"{Perturbation_path}_{samp}.txt", sep = "\t")
        df['Sample'] = samp
        combined_conditions.append(df)

    df = pd.concat(combined_conditions)
    
    return  df

def Compile_Association_sheet(Association_path, gene_num = 5):

    print('Load categorical associa data')

    df = pd.read_csv(Association_path, sep = "\t", index_col = 0)
    df = df.reset_index(drop=True)
    df.index.name = "program_name"
    df.sort_values(by = df.index.name, inplace = True, ascending=True)

    return df

def Compile_Explained_variance(Explained_Variance_path):

    print('Load explained variance data')

    df = pd.read_csv(Explained_Variance_path, sep = "\t", index_col = 0)
    df.reset_index(drop=False, inplace=True)
    df.drop("ProgramID", axis = 1, inplace=True)
    df.index = range(len(df))
    df.index.name = 'program_name'

    return df

#-------------- helper methods for loading sheets--------------

# combine methods for simple sheets
def load_simple_sheets(mdata, out_dir, run_name, k, sel_thresh, num_gene = 300, perturbation_file_name = "perturbation_association_results", Sample = ['1', '2', '3']):

    print('Load simple sheet')

    GO_path = f'{out_dir}/{run_name}/Eval/{k}_{str(sel_thresh).replace(".","_")}/{k}_GO_term_enrichment.txt'
    Geneset_path = f'{out_dir}/{run_name}/Eval/{k}_{str(sel_thresh).replace(".","_")}/{k}_geneset_enrichment.txt'
    Trait_path = f'{out_dir}/{run_name}/Eval/{k}_{str(sel_thresh).replace(".","_")}/{k}_trait_enrichment.txt'
    Perturbation_path_base = f'{out_dir}/{run_name}/Eval/{k}_{str(sel_thresh).replace(".","_")}/{k}_{perturbation_file_name}'
    Association_path = f'{out_dir}/{run_name}/Eval/{k}_{str(sel_thresh).replace(".","_")}/{k}_categorical_association_results.txt'
    Explained_Variance_path = f'{out_dir}/{run_name}/Eval/{k}_{str(sel_thresh).replace(".","_")}/{k}_Explained_Variance.txt'

    # compile program loadings
    df_Program_loading_long = compile_Program_loading_score_sheet_long(mdata, num_gene = num_gene)
    df_Program_loading_flat = compile_Program_loading_score_sheet_flat(mdata, num_gene = num_gene)

    # compile GO
    if os.path.exists(GO_path):
        df_GO = Compile_GO_sheet(GO_path, gene_num = num_gene)
    else:
        print(f'GO file not found: {GO_path}')
        df_GO = None

    # compile Genesets
    if os.path.exists(Geneset_path):
        df_Geneset = Compile_Geneset_sheet(Geneset_path, gene_num = num_gene)
    else:
        print(f'Geneset file not found: {Geneset_path}')
        df_Geneset = None

    # compile Trait
    if os.path.exists(Trait_path):
        df_Trait = Compile_Trait_sheet(Trait_path, gene_num = num_gene)
    else:
        print(f'Trait file not found: {Trait_path}')
        df_Trait = None

    # compile perturbation
    perturbation_files = [f"{Perturbation_path_base}_{samp}.txt" for samp in Sample]
    if any(os.path.exists(f) for f in perturbation_files):
        df_Perturbation = Compile_Perturbation_sheet(Perturbation_path_base, Sample = Sample)
    else:
        print(f'No perturbation files found for base: {Perturbation_path_base}')
        df_Perturbation = None

    # compile association
    if os.path.exists(Association_path):
        df_Association = Compile_Association_sheet(Association_path, gene_num = num_gene)
    else:
        print(f'Association file not found: {Association_path}')
        df_Association = None

    # compile explained variance
    if os.path.exists(Explained_Variance_path):
        df_Explained_Variance = Compile_Explained_variance(Explained_Variance_path)
    else:
        print(f'Explained variance file not found: {Explained_Variance_path}')
        df_Explained_Variance = None

    return df_Program_loading_long, df_Program_loading_flat, df_GO, df_Geneset, df_Trait, df_Perturbation, df_Association, df_Explained_Variance





# calcaulte specicicity scores

#-------------- helper methods for specicicity scores--------------

## helper methods for specicicity score calcualtion 
def Compute_regulator_zscore(perturb_df):

        # Pivot to get program × gene log2FC table
        pivot_table = perturb_df.pivot_table(
            index='target_name',
            columns='program_name',
            values='log2FC',
            aggfunc='first'
        )

        # convert to abs values
        # pivot_table_abs = pivot_table.abs()

        # Compute z-scores row-wise
        median = pivot_table.median(axis=0)
        mad = (pivot_table - median).abs().median(axis=0)
        epsilon = 1e-10

        z_score_table = (pivot_table - median) / (1.4826 * mad + epsilon)

        return z_score_table

def softmax_with_temperature(effects_df, T=1.0):

    # eq = e^{z/T}/ sum[ e^{z/T}]
    
    # Convert pandas to numpy and abs
    effects = effects_df.abs().values
    
    effects_scaled = effects / T
    effects_scaled = effects_scaled - effects_scaled.max(axis=1, keepdims=True) # for overflow prevention 
    
    exp_effects = np.exp(effects_scaled)
    weights = exp_effects / exp_effects.sum()
    
    # Convert back to DataFrame
    return pd.DataFrame(weights, 
                       index=effects_df.index, 
                       columns=effects_df.columns)

def compute_joint_distribution(weights):

    weights_array = weights.values
    
    # Joint: divide entire matrix by total sum
    P_tp = weights_array / weights_array.sum()
    
    # Marginals: sum by rows/columns
    P_t = P_tp.sum(axis=1)  # Row sums
    P_p = P_tp.sum(axis=0)  # Column sums

    # Convert back to DataFrame
    P_tp = pd.DataFrame(P_tp, 
                       index=weights.index, 
                       columns=weights.columns)

    
    return P_tp, P_t, P_p

def compute_pointwise_mutual_information(P_tp, P_t, P_p):

    # Reshape marginals for broadcasting
    P_t_col = P_t.reshape(-1, 1)  # Column vector
    P_p_row = P_p.reshape(1, -1)  # Row vector
    
    # P(t) * P(p) for all combinations
    P_independent = P_t_col @ P_p_row  # Outer product
    
    # Avoid log(0)
    epsilon = 1e-10
    
    # PMI = log(P_tp / P_independent)
    PMI = np.log((P_tp.values + epsilon) / (P_independent + epsilon))
    
    # Convert back to DataFrame for readability
    PMI_df = pd.DataFrame(PMI, 
                    index=P_tp.index, 
                    columns=P_tp.columns)
    
    return PMI_df

#-------------- helper methods specicicity scores--------------

# get specific programs and scores 
def get_specificity_program(Perturbation_path, Sample = ["D0", "sample_D1", "sample_D2", "sample_D3"], T=0.5, save_path = None):

    print("Calculate specificity scores for program")

    all_specificity = []

    for samp in Sample:
        df = pd.read_csv(f"{Perturbation_path}_{samp}.txt", sep="\t")
        # df = df[df['adj_pval'] <= adj_pval_threshold].copy()   

        z_score_table = Compute_regulator_zscore(df)                # compute into robust row z-score
        weights = softmax_with_temperature(z_score_table, T=T)    # compute softmax with T
        P_tp, P_t, P_p = compute_joint_distribution(weights)        # compute joint/marginal probability
        PMI = compute_pointwise_mutual_information(P_tp, P_t, P_p)  # compute pointwise MI

        if save_path is not None:
            PMI.to_csv(f'{save_path}/specificity_score_{samp}.txt', sep = '\t')

        program_specificity = []
        for i in range(len(PMI)):
            gene_name = PMI.index[i]
            row = PMI.iloc[i]

            top_5 = row.dropna().nlargest(5)
            top_5_programs = ', '.join(top_5.index.astype(str)) if len(top_5) > 0 else ''
            top_5_scores = ', '.join([f"{score:.4f}" for score in top_5.values]) if len(top_5) > 0 else ''

            program_specificity.append({
                'target_name': gene_name,
                f'top 5 specific programs {samp}': top_5_programs,
                f'top5 specificity scores {samp}': top_5_scores,
            })

        df_specificity = pd.DataFrame(program_specificity).set_index('target_name')
        all_specificity.append(df_specificity)

    # Merge all results by index (gene)
    merged_df = pd.concat(all_specificity, axis=1)
    
    return merged_df






# load target summary

#-------------- helper methods target summary--------------

# Get cells with the guide per day
def get_guide_cells_per_days(mdata, categorical_key="sample", data_key='rna', prog_key='cNMF',
                            guide_assignment_key="guide_assignment", guide_targets_key="guide_targets"):
                                                                                                                                                                                                                                                    
    print('Compute guides cells per days')
                                                                                                                                                                                                                                                    
    guide_assignment = mdata[prog_key].obsm[guide_assignment_key]  # (n_cells × n_targets), no copy                                                                                                                                              
    guide_targets = mdata[prog_key].uns[guide_targets_key]

    if not sparse.issparse(guide_assignment):
        guide_assignment = sparse.csc_matrix(guide_assignment)

    # Build sparse group indicator matrix (n_cells × n_groups)
    groups = mdata[data_key].obs[categorical_key]
    unique_groups = sorted(groups.unique())
    group_to_idx = {g: i for i, g in enumerate(unique_groups)}
    col_idx = groups.map(group_to_idx).values.astype(np.int32)

    indicator = sparse.csc_matrix(
        (np.ones(len(groups), dtype=np.float32), (np.arange(len(groups)), col_idx)),
        shape=(len(groups), len(unique_groups))
    )

    # Sparse matmul: (n_groups × n_cells) @ (n_cells × n_targets) = (n_groups × n_targets)
    group_sums = indicator.T @ guide_assignment

    # Build result: (n_targets × n_groups)
    df_merge = pd.DataFrame(
        np.asarray(group_sums.todense()).T,
        index=guide_targets,
        columns=[f'# Cells {g}' for g in unique_groups]
    )

    # Collapse duplicate target names
    df_merge = df_merge.groupby(df_merge.index).sum()
    df_merge.index.name = "target_name"

    return df_merge


# Get the gene expression average per day 
def get_guide_mean_expr_per_day(mdata, categorical_key = "sample", prog_key = 'cNMF', data_key = 'rna', guide_targets_key = "guide_targets"):
                                                                                                                                                                            
    print('Get targeted gene mean expressuion per day')                                                                                                                                      
                                                                                                                                                                            
    adata = mdata[data_key]  # no .copy() — avoid duplicating the whole object                                                                                           

    perturbed_genes = mdata[prog_key].uns[guide_targets_key]

    gene_mask = adata.var_names.isin(perturbed_genes)
    gene_names = adata.var_names[gene_mask]
    X_sub = adata[:, gene_mask].X  # keep sparse, only subset columns

    if not sparse.issparse(X_sub):
        X_sub = sparse.csc_matrix(X_sub)

    # Compute mean per group directly on the sparse matrix
    groups = adata.obs[categorical_key]
    mean_per_group = {}
    for group_name in groups.unique():
        row_mask = (groups == group_name).values
        mean_per_group[f'mean_expression_{group_name}'] = np.asarray(X_sub[row_mask].mean(axis=0)).ravel()

    mean_expr_per_day = pd.DataFrame(mean_per_group, index=gene_names)
    mean_expr_per_day.sort_index(inplace=True)
    mean_expr_per_day.index.name = 'target_name'

    return mean_expr_per_day

# Function to get significant programs for each gene across all days
def get_significant_programs(Perturbation_path, Sample = ["D0", "sample_D1","sample_D2","sample_D3"], adj_pval_threshold=0.05):

    print('Compute significant programs')
    
    significant_programs = {}
    
    for samp in Sample:
        # Read perturbation results for each day
        df = pd.read_csv(f"{Perturbation_path}_{samp}.txt", sep="\t")
        
        # Filter for significant associations
        significant = df[df['adj_pval'] < adj_pval_threshold]
        
        # Group by target_name and collect significant programs
        for target in significant['target_name'].unique():
            target_data = significant[significant['target_name'] == target]
            programs = target_data['program_name'].tolist()
            
            if target not in significant_programs:
                significant_programs[target] = {}
            
            significant_programs[target][samp] = programs
    
    return significant_programs

# get signifcant programs 
def get_significant_programs_df(Perturbation_path, Sample = ["D0", "sample_D1","sample_D2","sample_D3"], adj_pval_threshold=0.05):

    print('Compile significant programs')

    # Get significant programs
    significant_programs = get_significant_programs(Perturbation_path, Sample, adj_pval_threshold)

    # Convert to DataFrame format for easier viewing
    sig_prog_data = []

    for gene, days_data in significant_programs.items():
        row = {'target_name': gene}
        for samp in Sample:
            programs = days_data.get(samp, [])
            row[f'significant programs {samp}'] = ', '.join(map(str, programs)) if programs else ''
        sig_prog_data.append(row)

    df_significant_programs = pd.DataFrame(sig_prog_data)
    df_significant_programs = df_significant_programs.set_index('target_name')

    for samp in Sample:
        df_significant_programs[f'# programs {samp}'] = df_significant_programs[f'significant programs {samp}'].apply(
            lambda x: str(len(x.split(','))) if x and x != '' else 0
        )


    return df_significant_programs

# Get top correlation terms 
def get_correlation_df(perturbation_path, Sample=["D0", "sample_D1", "sample_D2", "sample_D3"], top_n=5, save_path = None):

    correlation_results_all_days = {}

    for day in Sample:
        print(f"Compute correlation for {day}")
        perturb_path = f"{perturbation_path}_{day}.txt"

        # Read and pivot directly
        df = pd.read_csv(perturb_path, sep='\t', index_col=0)
        pivot_df = df.pivot_table(index='target_name', columns='program_name', values='log2FC')

        # Compute correlation matrix (genes x genes) — vectorized
        corr_matrix = pivot_df.T.corr()

        # Zero out diagonal so self-correlation doesn't appear in top/bottom
        np.fill_diagonal(corr_matrix.values, np.nan)

        if save_path is not None:
            corr_matrix.to_csv(f"{save_path}/corr_matrix_{day}.txt", sep="\t")
            corr_matrix.to_csv(f"{save_path}/corr_matrix_{day}.txt.gz", sep="\t", compression="gzip")

        genes = corr_matrix.index.values
        vals = corr_matrix.values  # numpy array

        # Use argpartition for O(n) top/bottom selection instead of full sort
        rows = []
        for i in range(len(genes)):

            # for each row, remove NaN and get valid genes 
            row = vals[i]
            valid = ~np.isnan(row)
            valid_genes = genes[valid]

            # Top positive
            top_idx = np.argpartition(row[valid], -(top_n-1))[-top_n:] # get top 5 value's index (-1 for accessing k-th value)
            top_idx = top_idx[np.argsort(row[valid][top_idx])[::-1]]   # get top 5 values and sort them in right otder 

            # Bottom negative
            bot_idx = np.argpartition(row[valid], (top_n-1))[:top_n]  # get bottom 5 value's index (-1 for accessing k-th value)
            bot_idx = bot_idx[np.argsort(row[valid][bot_idx])]        # get bottom 5 values and sort them in right otder 

            rows.append({
                'target_name': genes[i],
                f'top 5 pos correls targets (program log2fc) {day}': '; '.join(str(g) for g in valid_genes[top_idx]),
                f'top 5 pos correls (program log2fc) {day}': '; '.join(f"{v:.3f}" for v in row[valid][top_idx]),
                f'top 5 neg correls targets (program log2fc) {day}': '; '.join(str(g) for g in valid_genes[bot_idx]),
                f'top 5 neg correls (program log2fc) {day}': '; '.join(f"{v:.3f}" for v in row[valid][bot_idx]),
            })

        correlation_results_all_days[day] = pd.DataFrame(rows).set_index('target_name')



    return pd.concat(correlation_results_all_days.values(), axis=1)


#-------------- helper methods starget summary--------------

# final function to compile target Summary sheet 
def Compile_Target_Summary_sheet(mdata, perturbation_path, Sample = ["D0", "sample_D1","sample_D2","sample_D3"], adj_pval_threshold= 0.05, 
top_n=5, T= 0.5 , categorical_key = "sample", prog_key = 'cNMF', data_key = 'rna', guide_targets_key = "guide_targets", guide_assignment_key = "guide_assignment", save_path = None):

    print('Generate target summary sheet')


    df_mean_expr_per_day = get_guide_mean_expr_per_day(mdata, categorical_key=categorical_key, prog_key=prog_key, data_key=data_key, 
    guide_targets_key=guide_targets_key)
    df_guide_days = get_guide_cells_per_days(mdata,  categorical_key=categorical_key, guide_assignment_key=guide_assignment_key, prog_key=prog_key, 
    data_key=data_key, guide_targets_key=guide_targets_key)
    df_significant_program = get_significant_programs_df(perturbation_path ,Sample=Sample, adj_pval_threshold=adj_pval_threshold)
    df_specificity_program =  get_specificity_program(perturbation_path, Sample=Sample, T = T, save_path = save_path) # save_path for saving specificity scores 
    df_correlation = get_correlation_df(perturbation_path ,Sample=Sample, top_n=top_n, save_path = save_path)

    final_merged_df = pd.merge(
      df_mean_expr_per_day,
      df_guide_days,
      left_index=True,
      right_index=True,
      how='outer'
    )

    final_merged_df = pd.merge(                                                                                                                                                                            
        final_merged_df,                                                                                                                                                                                  
        df_significant_program,                                                                                                                                                                          
        left_index=True,                                                                                                                                                                                 
        right_index=True,                                                                                                                                                                                 
        how='outer'                                                                                                                                                                                       
     ) 

    final_merged_df = pd.merge(                                                                                                                                                                            
      final_merged_df,                                                                                                                                                                                  
      df_specificity_program,                                                                                                                                                                          
       left_index=True,                                                                                                                                                                                 
       right_index=True,                                                                                                                                                                                 
       how='outer'                                                                                                                                                                                       
    )
    
    final_merged_df = pd.merge(                                                                                                                                                                            
      final_merged_df,                                                                                                                                                                                  
      df_correlation,                                                                                                                                                                          
       left_index=True,                                                                                                                                                                                 
       right_index=True,                                                                                                                                                                                 
       how='outer'                                                                                                                                                                                       
    ).fillna('')

    return final_merged_df






# load summary

#-------------- helper methods summary--------------

# get simply items ready in the summary sheet
def simple_Summary_cols(df, k, df_GO, df_Perturbation, df_Program_loading, df_Explained_Variance = None, Sample = ["D0", "sample_D1","sample_D2","sample_D3" ]
, non_tagerting_key = None):

    # create GO summary col
    if df_GO is not None:
        df_GO_enriched = df_GO.loc[df_GO['Adjusted P-value']<=0.05]
        df['Total Enriched GO Terms'] = [df_GO_enriched[df_GO_enriched['program_name']==i].shape[0] for i in range(k)] 

    # remove non-targeting off the list of perturbed genes 
    if df_Perturbation is not None:
        if non_tagerting_key is not None:
            df_Perturbation = df_Perturbation[~df_Perturbation['target_name'].isin(non_tagerting_key)]

        # compute # regulator per condition for + and -
        conditions = (df_Perturbation['Sample']).unique()
        for condition in conditions:
            df_Perturbation_ = df_Perturbation[df_Perturbation['Sample'] == condition]

            # create perturbation program summary col
            df_Perturbation_enriched = df_Perturbation_.loc[df_Perturbation_['adj_pval']<=0.05]
            df_Perturbation_positive = df_Perturbation_enriched.loc[df_Perturbation_enriched['log2FC'] > 0 ]
            df_Perturbation_negative = df_Perturbation_enriched.loc[df_Perturbation_enriched['log2FC'] < 0 ]

            df[f'Regulators with positive effect {condition}'] = [df_Perturbation_positive[df_Perturbation_positive['program_name']==i].shape[0] for i in range(k)] 
            df[f'Regulators with negative effect {condition}'] = [df_Perturbation_negative[df_Perturbation_negative['program_name']==i].shape[0] for i in range(k)] 

        # create perturbation gene summary col
        for condition in conditions:
            df_Perturbation_enriched = df_Perturbation.loc[df_Perturbation['adj_pval']<=0.05]
            df_Perturbation_D = df_Perturbation_enriched.loc[df_Perturbation_enriched['Sample'] == condition] # make for each condition

            targets_list = []
            for i in range(k):
                matching = df_Perturbation_D.loc[df_Perturbation_D['program_name'] == i] # for each program
                unique_indices = matching["target_name"].unique()                                  # find unique programs
                joined = ';'.join([str(x) for x in unique_indices])                      # join them 
                targets_list.append(joined)

            df[f'sigfdr0.05_targets_sorted_abslog2fcd_{condition}'] = targets_list

    # create motif summary col
    df['Total Enriched Enhancer Motifs'] = [''] * k
    df['Total Enriched Promoter Motifs'] = [''] * k


    # create top gene summary col
    if df_Program_loading is not None:
        df['top10_loaded_genes'] = [';'.join(df_Program_loading.iloc[i][:10]) for i in range(k)]

    # create explained variance summary col
    if df_Explained_Variance is not None:
        df['variance_explained'] = df_Explained_Variance

# make the program info in summary sheet
def get_program_info_Summary_cols(mdata, categorical_key = "sample"):

    # create cell info col summary
    df_cell = pd.DataFrame(data=mdata['cNMF'].X, index=mdata['cNMF'].obs_names)
    results = []

    # program #
    k = df_cell.shape[1]

    # Loop through all k values
    for i in range(k):
        df_cell_program = pd.DataFrame({
            "expression": df_cell.iloc[:, i],
            "cell_type": mdata['rna'].obs[categorical_key].values  
        })
        
        df_mean = df_cell_program.groupby("cell_type")["expression"].mean()
        df_frac = df_cell_program.groupby("cell_type")["expression"].apply(
            lambda x: (x > x.mean()).mean()
        )
        
        # Store as row with program as index
        results.append({
            'program_name': i,
            **{f'Mean program score {ct}': df_mean[ct] for ct in df_mean.index},
            **{f'Fra cells above mean program score {ct}': df_frac[ct] for ct in df_frac.index}
        })

    df_cell_info = pd.DataFrame(results).set_index('program_name')

    return df_cell_info

# make top terms for the summray sheet
def get_top_terms_Summary_cols(df_GO,df_Geneset):

    top_GO = None
    top_Geneset = None

    if df_GO is not None:
        top_GO = df_GO.groupby('program_name').apply(
            lambda x: ';'.join(x.sort_values('Adjusted P-value').index[:10])
        )
    if df_Geneset is not None:
        top_Geneset = df_Geneset.groupby('program_name').apply(
                lambda x: ';'.join(x.sort_values('Adjusted P-value').index[:10])
            )

    top_terms = pd.DataFrame({
        'top10_enriched_genesets': top_Geneset,
        'top10_enriched_go_terms': top_GO
    })

    return top_terms

#-------------- helper methods summary--------------

# compile summry sheet
def Compile_Summary_sheet(mdata, df_GO, df_Geneset, df_Perturbation, df_Program_loading, df_Explained_Variance , Sample = ["D0", "sample_D1","sample_D2","sample_D3"],
categorical_key = "sample",non_tagerting_key=None):

    print('complie summary sheet')

    # set program # from mdata (works even if df_GO is None)
    k = mdata['cNMF'].varm["loadings"].shape[0]

    df = pd.DataFrame({
    'manual_annotation_label': [''] * k,
    'manual_timepoint': [''] * k,
    'Notes': [''] * k,
    'Automatic Timepoint': [''] * k }, index=pd.Index(range(k), name='program_name'))

    simple_Summary_cols(df, k, df_GO, df_Perturbation, df_Program_loading, df_Explained_Variance,  Sample = Sample, non_tagerting_key=non_tagerting_key)
    df_cell_info_cols = get_program_info_Summary_cols(mdata,categorical_key)
    df_top_terms = get_top_terms_Summary_cols(df_GO, df_Geneset)

    # refill automatic time point
    col_names = [f'Mean program score {samp}' for samp in Sample]
    col_mapping = {f'Mean program score {samp}': samp for samp in Sample}      # Create mapping from column names to condition A, B, C
    df_mean = df_cell_info_cols[col_names]
    
    df['Automatic Timepoint'] = df_mean.idxmax(axis=1).map(col_mapping)

    merged_df = pd.merge(
        df,
        df_cell_info_cols,
        left_index=True,
        right_index=True,
        how='outer'
    )

    final_merged_df = pd.merge(
        merged_df,
        df_top_terms,
        left_index=True,
        right_index=True,
        how='outer'
    )

    # move top 10 loaded gene left of notes
    cols = final_merged_df.columns.tolist()                                                                                                                                                          
    cols.remove('top10_loaded_genes')
    i = cols.index('Notes')
    cols.insert(i, 'top10_loaded_genes')
    final_merged_df = final_merged_df[cols]

    return final_merged_df
   

# merge to have specificity score added for perturbation results
def add_specificity_scores(save_path, Perturbation_path_base, samp):                                                                                                                                                    
                                                                                                                                                                                                                     
      PMI = pd.read_csv(f'{save_path}/specificity_score_{samp}.txt', sep='\t', index_col=0)                                                                                                                          
      df_perturbation = pd.read_csv(f'{Perturbation_path_base}_{samp}.txt', sep='\t')                                                                                                                                     
                                                                                                                                                                                                                     
      # convert PMI to long form
      df_PMI = PMI.reset_index().melt(id_vars="target_name", var_name="program_name", value_name="specificity_scores")

      # make sure program name is in str form
      df_PMI["program_name"] = df_PMI["program_name"].astype(str)
      df_perturbation["program_name"] = df_perturbation["program_name"].astype(str)

      # merge to have specificity score added for perturbation results
      df_perturbation_merged = df_PMI.merge(df_perturbation, on=["target_name", "program_name"])

      df_perturbation_merged.to_csv(f'{save_path}/perturbation_merged_{samp}.txt', sep='\t')

      return df_perturbation_merged
