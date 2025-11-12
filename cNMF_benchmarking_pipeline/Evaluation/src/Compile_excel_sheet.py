import pandas as pd
import muon as mu
from scipy import sparse
import sys
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
def compile_Program_loading_score_sheet_long(mdata, num_gene = 300):

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
        summary_dict = {r.get('query'): r.get('summary', 'N/A') for r in results}
        
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

    df = pd.read_csv(GO_path, sep = "\t", index_col = 0)
    df = df.reset_index().set_index("Term")

    df["Genes"] = df["Genes"].str.split(';').str[:gene_num].str.join(";")

    return df

def Compile_Geneset_sheet(Geneset_path, gene_num = 5):

    df = pd.read_csv(Geneset_path, sep = "\t", index_col = 0)
    df = df.reset_index().set_index("Term")

    df["Genes"] = df["Genes"].str.split(';').str[:gene_num].str.join(";")

    return df

def Compile_Trait_sheet(Trait_path, gene_num = 5):

    df = pd.read_csv(Trait_path, sep = "\t", index_col = 0)
    df = df.reset_index().set_index("Term")

    df["Genes"] = df["Genes"].str.split(';').str[:gene_num].str.join(";")

    return df

def Compile_Perturbation_sheet(Perturbation_path, Sample = ["D0", "sample_D1","sample_D2","sample_D3" ]):

    combined_conditions = []
    for samp in Sample: 
        df = pd.read_csv(f"{Perturbation_path}_{samp}.txt", sep = "\t", index_col = 0)
        df['Sample'] = samp
        combined_conditions.append(df)

    df = pd.concat(combined_conditions)
    
    return  df

def Compile_Association_sheet(Association_path, gene_num = 5):

    df = pd.read_csv(Association_path, sep = "\t", index_col = 0)
    df = df.reset_index(drop=True)
    df.index.name = "program_name"
    df.sort_values(by = df.index.name, inplace = True, ascending=True)

    return df

def Compile_Explained_variance(Explained_Variance_path):
    df = pd.read_csv(Explained_Variance_path, sep = "\t", index_col = 0)
    df.reset_index(drop=False, inplace=True)
    df.drop("ProgramID", axis = 1, inplace=True)
    df.index = range(len(df))
    df.index.name = 'program_name'

    return df



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
    #effects_scaled = effects_scaled - effects_scaled.max(axis=1, keepdims=True) # for overflow prevention 
    
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
def get_specificity_program(Perturbation_path, Sample = ["D0", "sample_D1", "sample_D2", "sample_D3"], T=0.5):

    all_specificity = []

    for samp in Sample:
        df = pd.read_csv(f"{Perturbation_path}_{samp}.txt", sep="\t")
        # df = df[df['adj_pval'] <= adj_pval_threshold].copy()   

        z_score_table = Compute_regulator_zscore(df)                # compute into robust row z-score
        weights = softmax_with_temperature(z_score_table, T=T)    # compute softmax with T
        P_tp, P_t, P_p = compute_joint_distribution(weights)        # compute joint/marginal probability
        PMI = compute_pointwise_mutual_information(P_tp, P_t, P_p)  # compute pointwise MI

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
def get_guide_cells_per_days(adata):

    X = adata.obsm['guide_assignment'].T
    df = pd.DataFrame(X.toarray if hasattr(X, "toarray") else X,
                    index =  adata.uns["guide_targets"],
                    columns = adata.obs["sample"]
                    )
                    
    df_merge = df.groupby(df.index).sum()
    df_merge = df_merge.groupby(df_merge.columns, axis = 1).sum()
    df_merge.index.name = "target_name"

    df_merge.columns = [f'# Cells D{col}' for col in range(0,df_merge.shape[1])]
    
    return df_merge

# Get the gene expression average per day 
def get_guide_mean_expr_per_day(adata):

    # Extract perturbed genes and calculate mean expression across days
    perturbed_genes = adata.uns["guide_targets"]

    # Get expression data for perturbed genes
    gene_mask = adata.var_names.isin(perturbed_genes)
    perturbed_gene_expr = adata[:, gene_mask].X

    # Convert to dense if sparse
    if hasattr(perturbed_gene_expr, "toarray"):
        perturbed_gene_expr = perturbed_gene_expr.toarray()

    # Create DataFrame with gene expression
    expr_df = pd.DataFrame(
        perturbed_gene_expr.T,
        index=adata.var_names[gene_mask],
        columns=adata.obs.index
    )

    # Add sample information and calculate mean per day
    expr_df_with_samples = expr_df.T
    expr_df_with_samples['sample'] = adata.obs['sample'].values

    # Calculate mean expression per gene per day
    mean_expr_per_day = expr_df_with_samples.groupby('sample').mean().T

    # Rename
    mean_expr_per_day.columns = [f'mean_expression_{col}' for col in mean_expr_per_day.columns]   
    mean_expr_per_day.sort_index(inplace = True)
    mean_expr_per_day.index.name = 'target_name'

    return mean_expr_per_day

# Function to get significant programs for each gene across all days
def get_significant_programs(Perturbation_path, Sample = ["D0", "sample_D1","sample_D2","sample_D3"], adj_pval_threshold=0.05):
    
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

    # Get significant programs
    significant_programs = get_significant_programs(Perturbation_path, Sample, adj_pval_threshold)

    # Convert to DataFrame format for easier viewing
    sig_prog_data = []

    for gene, days_data in significant_programs.items():
        row = {'target_name': gene}
        for day in ["D0", "sample_D1", "sample_D2", "sample_D3"]:
            programs = days_data.get(day, [])
            row[f'significant programs {day}'] = ', '.join(map(str, programs)) if programs else ''
        sig_prog_data.append(row)

    df_significant_programs = pd.DataFrame(sig_prog_data)
    df_significant_programs = df_significant_programs.set_index('target_name')

    df_significant_programs['# programs D0'] = df_significant_programs['significant programs D0'].apply(
        lambda x: str(len(x.split(','))) if x and x != '' else 0
    )

    df_significant_programs['# programs D1'] = df_significant_programs['significant programs sample_D1'].apply(
        lambda x: str(len(x.split(','))) if x and x != '' else  0
    )

    df_significant_programs['# programs D2'] = df_significant_programs['significant programs sample_D2'].apply(
        lambda x: str(len(x.split(','))) if x and x != '' else 0
    )

    df_significant_programs['# programs D3'] = df_significant_programs['significant programs sample_D3'].apply(
        lambda x: str(len(x.split(','))) if x and x != '' else 0
    )

    return df_significant_programs

# Get top correlation terms 
def get_correlation_df(perturbation_path, days=["D0", "sample_D1", "sample_D2", "sample_D3"], top_n=5):
    """
    Create correlation table with top/bottom correlations for each gene across days
    """
    
    correlation_results_all_days = {}

    for day in days:

        correlation_results = []
        perturb_path = f"{perturbation_path}_{day}.txt"
        
        # Compute correlations for this day
        waterfall_correlation = compute_gene_waterfall_cor(perturb_path)
        
        # Use the keys from waterfall_correlation as the gene list
        for gene in waterfall_correlation.keys():
            gene_correlations = waterfall_correlation[gene]
            corr_series = pd.Series(gene_correlations).sort_values(ascending=False)
            
            # Get top 5 positive and bottom 5 negative correlations
            top_positive = corr_series.head(top_n)
            top_negative = corr_series.tail(top_n)
            
            # Format as strings
            top_pos_targets = '; '.join(top_positive.index.tolist())
            top_pos_values = '; '.join([f"{x:.3f}" for x in top_positive.values])
            top_neg_targets = '; '.join(top_negative.index.tolist())
            top_neg_values = '; '.join([f"{x:.3f}" for x in top_negative.values])
            
            correlation_results.append({
                'target_name': gene,
                f'top 5 pos correls targets (program log2fc) {day}': top_pos_targets,
                f'top 5 pos correls (program log2fc) {day}': top_pos_values,
                f'top 5 neg correls targets (program log2fc) {day}': top_neg_targets,
                f'top 5 neg correls (program log2fc) {day}': top_neg_values
            })
    
             # Convert to DataFrame and reorganize by gene
        
        correlation_results_all_days[day] = pd.DataFrame(correlation_results).set_index('target_name') 
            
    
    # Pivot to get all day columns for each gene
    final_df = pd.concat(correlation_results_all_days.values(), axis=1)
    
    return final_df

#-------------- helper methods starget summary--------------


# final function to compile target Summary sheet 
def Compile_Target_Summary_sheet(adata, perturbation_path, Sample = ["D0", "sample_D1","sample_D2","sample_D3"], adj_pval_threshold= 0.05, top_n=5, T= 0.5):

    df_mean_expr_per_day = get_guide_mean_expr_per_day(adata)
    df_guide_days = get_guide_cells_per_days(adata)
    df_significant_program = get_significant_programs_df(perturbation_path ,Sample, adj_pval_threshold)
    df_specificity_program =  get_specificity_program(perturbation_path, Sample, T = T)
    df_correlation = get_correlation_df(perturbation_path ,Sample, top_n)

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
def simple_Summary_cols(df, df_GO, df_Perturbation, df_Program_loading, df_Explained_Variance, Sample = ["D0", "sample_D1","sample_D2","sample_D3" ]):

    # set program #
    k = len(df_GO['program_name'].unique())

    # create GO summary col
    df_GO_enriched = df_GO.loc[df_GO['Adjusted P-value']<=0.05]
    df['Total Enriched GO Terms'] = [df_GO_enriched[df_GO_enriched['program_name']==i].shape[0] for i in range(k)] 

    # create perturbation program summary col
    df_Perturbation_enriched = df_Perturbation.loc[df_Perturbation['adj_pval']<=0.05]
    df_Perturbation_positive = df_Perturbation_enriched.loc[df_Perturbation_enriched['log2FC'] > 0 ]
    df_Perturbation_negative = df_Perturbation_enriched.loc[df_Perturbation_enriched['log2FC'] < 0 ]

    df['Regulators with positive effect'] = [df_Perturbation_positive[df_Perturbation_positive['program_name']==i].shape[0] for i in range(k)] 
    df['Regulators with negative effect'] = [df_Perturbation_negative[df_Perturbation_negative['program_name']==i].shape[0] for i in range(k)] 

    # create motif summary col
    df['Total Enriched Enhancer Motifs'] = [''] * k
    df['Total Enriched Promoter Motifs'] = [''] * k


    # create top gene summary col
    df['top10_loaded_genes'] = [';'.join(df_Program_loading.iloc[i][:10]) for i in range(k)]

    # create explained variance summary col
    df['variance_explained'] = df_Explained_Variance

    # create perturbation gene summary col
    for samp in Sample:
        df_Perturbation_D = df_Perturbation_enriched.loc[df_Perturbation_enriched['Sample'] == samp]
        df_Perturbation_D_program = df_Perturbation_D.loc[df_Perturbation_D['program_name'] == 0]
        df[f'sigfdr0.05_targets_sorted_abslog2fc_{samp}'] = [';'.join(df_Perturbation_D.loc[df_Perturbation_D['program_name'] == i].index.unique()) for i in range(k)]

# make the program info in summary sheet
def get_program_info_Summary_cols(mdata):

    # create cell info col summary
    df_cell = pd.DataFrame(data=mdata['cNMF'].X, index=mdata['cNMF'].obs_names)
    results = []

    # program #
    k = df_cell.shape[1]

    # Loop through all k values
    for i in range(k):
        df_cell_program = pd.DataFrame({
            "expression": df_cell.iloc[:, i],
            "cell_type": mdata['rna'].obs["sample"].values  
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

    top_GO = df_GO.groupby('program_name').apply(
        lambda x: ';'.join(x.sort_values('Adjusted P-value').index[:10])
    )

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
def Compile_Summary_sheet(mdata, df_GO, df_Geneset, df_Perturbation, df_Program_loading, df_Explained_Variance, Sample = ["D0", "sample_D1","sample_D2","sample_D3"]):

    # set program #
    k = len(df_GO['program_name'].unique())

    df = pd.DataFrame({
    'manual_annotation_label': [''] * k,
    'manual_timepoint': [''] * k,
    'Notes': [''] * k,
    'Automatic Timepoint': [''] * k }, index=pd.Index(range(k), name='program_name'))

    simple_Summary_cols(df, df_GO, df_Perturbation, df_Program_loading, df_Explained_Variance,  Sample = ["D0", "sample_D1","sample_D2","sample_D3"])
    df_cell_info_cols = get_program_info_Summary_cols(mdata)
    df_top_terms = get_top_terms_Summary_cols(df_GO, df_Geneset)

    # refill automatic time point
    df_mean = df_cell_info_cols[['Mean program score D0', 'Mean program score sample_D1','Mean program score sample_D2', 'Mean program score sample_D3']]
    
    # Create mapping from column names to 0,1,2,3
    col_mapping = {
        'Mean program score D0': 0,
        'Mean program score sample_D1': 1,
        'Mean program score sample_D2': 2,
        'Mean program score sample_D3': 3
    }

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
    

    return final_merged_df
    
