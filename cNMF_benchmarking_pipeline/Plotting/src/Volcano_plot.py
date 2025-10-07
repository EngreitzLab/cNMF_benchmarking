from bioinfokit import visuz
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import xarray as xr
from pathlib import Path

# plot one volcone plot given file path
def plot_volcano(path, down_thred_log = -0.05, up_thred_log = 0.05, p_value_thred = 0.05, program_num = 0, folder_name = None, file_name = None):

    df = pd.read_csv(path, sep = "\t")
    df_program = df.loc[df["program_name"] == program_num]

    plt.scatter(x=df_program['log2FC'],y=df_program['adj_pval'].apply(lambda x:-np.log10(x)),s=1,label="Not significant")

    # highlight down- or up- regulated genes
    down = df_program[(df_program['log2FC']<=down_thred_log)&(df_program['adj_pval']<=p_value_thred)]
    up = df_program[(df_program['log2FC']>=up_thred_log)&(df_program['adj_pval']<=p_value_thred)]

    plt.scatter(x=down['log2FC'],y=down['adj_pval'].apply(lambda x:-np.log10(x)),s=3,label="Down-regulated",color="blue")
    plt.scatter(x=up['log2FC'],y=up['adj_pval'].apply(lambda x:-np.log10(x)),s=3,label="Up-regulated",color="red")

    for i,r in up.iterrows():
        plt.text(x=r['log2FC'],y=-np.log10(r['adj_pval']),s=r['target_name'])

    for i,r in down.iterrows():
        plt.text(x=r['log2FC'],y=-np.log10(r['adj_pval']),s=r['target_name'])
        
    plt.xlabel("log2FC")
    plt.ylabel("adj_pval")
    plt.axvline(down_thred_log,color="grey",linestyle="--")
    plt.axvline(up_thred_log,color="grey",linestyle="--")
    plt.axhline(-np.log10(p_value_thred),color="grey",linestyle="--")
    plt.legend()
    plt.title(f"volcano plot for program {program_num} on {file_name}")

    if folder_name and file_name:
        fig.savefig(f"{folder_name}/{file_name}.png")

    plt.show()
    plt.close()

# plot volcano plots by date
def plot_all_days_valcano(in_folder_name, in_file_name,  down_thred_log = -0.05, up_thred_log = 0.05, p_value_thred = 0.05, program_num = 0, folder_name = None, file_name = None):
    for samp in ['D0', 'sample_D1', 'sample_D2', 'sample_D3']:
        path = f"{in_folder_name}/{in_file_name}_{samp}_perturbation_association.txt"
        plot_volcano(path, down_thred_log=down_thred_log, up_thred_log=up_thred_log, p_value_thred=p_value_thred, program_num=program_num,folder_name=folder_name, file_name = in_file_name)