# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 14:15:31 2021
This script plots heatmaps. 

@author: Shuyue Xue
"""
import os, matplotlib
import pandas as pd, matplotlib.pyplot as plt, seaborn as sns
from preprocess_raw_data import configure_SLV_dataframe, configure_Bcell_datasets
import dir_conf

# set filenames and their paths
# dir_conf.initialize('SLV')  # test code with SLV data
cd = os.getcwd()
data_set = dir_conf.data_set
network_data_path = os.path.join(dir_conf.results_dir, 'network_data')
communities_data_path = os.path.join(dir_conf.results_dir, 'reactome_analysis')
comm_gene_file_path = os.path.join(cd, network_data_path,
                               'M0.125_Qf0.75_Qc0.995_DN_full_name_communities.csv')
sgnf_communities_file_path = os.path.join(cd, communities_data_path,
                               'M0.125_Qf0.75_Qc0.995_DN_community_reactome.xlsx')


# load the gene expression files
os.chdir(dir_conf.data_path)
if 'SLV' == data_set:
    du, dt = configure_SLV_dataframe()
    treated_file = 'community_signals_h2(treated).xlsx'
    untreated_file = 'community_signals_h1(untreated).xlsx'
if 'Bcell' == data_set:
    drt, dru, dt, du = configure_Bcell_datasets()
    treated_file = 'community_signals_treated.xlsx'
    untreated_file = 'community_signals_untreated.xlsx'
    
# load the community files
comm_genes = pd.read_csv(comm_gene_file_path, header=0, index_col=0)
comm_xls = pd.ExcelFile(sgnf_communities_file_path)
comms_for_plot = comm_xls.sheet_names

# change to the heatmap directory to export the plots
os.chdir(os.path.join(cd, dir_conf.results_dir, 'network_plots'))
if not os.path.exists('heatmaps'):
    os.makedirs('heatmaps')
os.chdir('./heatmaps/')

# group gene signals by communities with significant pathways (export in 2 files)
with pd.ExcelWriter(treated_file, engine='xlsxwriter') \
    as writer:
    for comm in comms_for_plot:
        this_comm_genes = comm_genes[comm].dropna()
    
        if 'SLV' == data_set:
            comm_signals = dt.loc[ this_comm_genes ]
        
        if 'Bcell' == data_set:
            comm_signals_t = dt.loc[ this_comm_genes ]
            comm_signals_rt = drt.loc[ this_comm_genes ]            
            comm_signals = 1/2*(comm_signals_t + comm_signals_rt)

        comm_signals.to_excel(writer, sheet_name=comm)

with pd.ExcelWriter(untreated_file, engine='xlsxwriter') \
    as writer:
    for comm in comms_for_plot:
        this_comm_genes = comm_genes[comm].dropna()
        
        if 'SLV' == data_set:
            comm_signals = du.loc[ this_comm_genes ]
        
        if 'Bcell' == data_set:
            comm_signals_u = du.loc[ this_comm_genes ]
            comm_signals_ru = dru.loc[ this_comm_genes ] 
            comm_signals = 1/2*(comm_signals_u + comm_signals_ru)
        
        comm_signals.to_excel(writer, sheet_name=comm)
        
        
def plot_heatmap_and_time_average(community_signal_normalized_df, raw_df, title):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 13))
    plt.rcParams['font.size'] = 14
    plt.rcParams['font.family'] = 'serif'
    
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('CP',['steelblue','steelblue','steelblue','w','firebrick','firebrick','firebrick'], N=100)
    
    heatmap = sns.clustermap(community_signal_normalized_df,
               metric = 'euclidean',
               method = 'complete',
               yticklabels=False,
               col_cluster=False,
               cmap = cmap,
               vmin=-1, vmax=1,

               cbar_pos=(0.93, 0.06, 0.014, 0.91),
                )
    hm = heatmap.ax_heatmap.get_position()
    heatmap.ax_heatmap.set_position([hm.x0, hm.y0, hm.width*0.86, hm.height*1.25])
    
    row = heatmap.ax_row_dendrogram.get_position()
    heatmap.ax_row_dendrogram.set_position([row.x0, row.y0,
                                            row.width, row.height*1.25])
    
    file_name = title.replace(':','')
    # file_name = file_name.replace(' ','_')+'.svg'
    file_name = file_name.replace(' ','_')+'.jpeg'
    plt.savefig( file_name )
    plt.close()

    return


def normalize_df(raw_signal_df):
    """normalize the df"""
    from sklearn.preprocessing import Normalizer
    # 1. subtract t=0 from all time points
    # normalized_df =raw_signal_df.iloc[:, 1:].sub(df[0], axis=0) # keep t=0
    normalized_df = raw_signal_df.sub(raw_signal_df[0], axis=0)   # norm t=0

    # 2. normalize as vector norm
    normalized_df.iloc[:,:] = \
                    Normalizer(norm='l2').fit_transform(normalized_df)
    return normalized_df


treated_xls = pd.ExcelFile(treated_file)
untreated_xls = pd.ExcelFile(untreated_file)
for comm in treated_xls.sheet_names:
    # print('\niteration %s' %comm)
    df_t = pd.read_excel(treated_xls, index_col=0, sheet_name=comm)
    df_ut = pd.read_excel(untreated_xls, index_col=0, sheet_name=comm)
        
    norm_df_t = normalize_df(df_t)
    norm_df_ut = normalize_df(df_ut)

    diff = norm_df_t - norm_df_ut  
    diff.columns.name = 'Time (hours)'
    title = 'Differentail Signal: ' + comm + ' Heatmap'
    plot_heatmap_and_time_average(diff, diff, title=title)


os.chdir(cd) # Back to the Code directory