# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 15:31:41 2020

Test: 
    * plot the histograms of 1) aggregated the number of FDR <0.05 in all communities in the community-wise results 2) the number of FDR <0.05 in Network-wise results
    * To see if there's a postive correlation between the two
    * Then we will stick with one plot for number of FDR historgram plot

Plot the ration of entites found / entites total as a second parameter to determine the "good" vs. "bad" results

@author: Shuyue Xue
"""
import os, glob, matplotlib
import pandas as pd, numpy as np, matplotlib.pyplot as plt
from matplotlib import transforms

import pyiomica as pio
from pyiomica.enrichmentAnalyses import KEGGAnalysis, GOAnalysis, \
ExportEnrichmentReport, ReactomeAnalysis, ExportReactomeEnrichmentReport

"""helper: count the number of significant pathways by FDR"""
num_significant_pathways = lambda df: np.sum( df['Entities FDR'] <=0.05 )


"""Return: total number of sig_path across all communties + community-wise #sig_paths."""
def get_reactome(dn, filename='react'):
    """
    Does reactome enrichment anaylsis and writes community results into excel sheets.
    Return: 
        total number of sig_path across all communties . 
    """
    if not os.path.exists('reactome_analysis'):
                os.makedirs('reactome_analysis')
    
    # remove the duplicates:
    dn = dn.stack()
    dn.drop_duplicates(keep='first', inplace=True)
    dn = dn.unstack()
    
    # set a counter for the total number of sig_pathwas in all communities
    total_sig_paths = 0      
    community_sig_paths = {}
    max_comm_sig_paths = -1                  # get the community with most significant FDRs
     
    writer  = pd.ExcelWriter('./reactome_analysis/'+ filename + \
                        '_community_reactome.xlsx', engine='xlsxwriter')
    for colm in dn.columns:
        if dn[colm].count() < 8: continue  # Ignore communities whose size <=8 genes 
        rect = ReactomeAnalysis( dn[colm].dropna().values.tolist() )
        if not rect.empty:
                rect.to_excel(writer, sheet_name=colm)
                
                this_comm_sig_paths = num_significant_pathways(rect)
                total_sig_paths += this_comm_sig_paths
                community_sig_paths[colm] = this_comm_sig_paths
                
                if max_comm_sig_paths < this_comm_sig_paths: 
                    max_comm_sig_paths = this_comm_sig_paths
        else:
            print(colm+' column analysis is empty')        
    writer.save()
    writer.close()
    
    non_empty_comm_sig_paths = {k:v for k,v in community_sig_paths.items() if 0 != v}
    
    if max_comm_sig_paths > 0:
        max_comms = [k for k, v in community_sig_paths.items() 
                     if v == max_comm_sig_paths]
        print('Communities with the most significant FDRs are ', max_comms)
    else:
        max_comms = []
        print('None of the %s communities have significant pathways' 
              %len(dn.columns))
    
    return total_sig_paths, community_sig_paths, non_empty_comm_sig_paths, max_comms


def plot_ntwk_enrichment_overview(non_empty_comm_sig_paths, total_sig_paths, 
                                  ax='', figname=''):
    """ 
    Plot the bar chart of the # of siganificant FDR in each non-empty community. 
    """
    os.chdir('./reactome_analysis')
    if not ax:
        fig, ax = plt.subplots()
    plt.bar(non_empty_comm_sig_paths.keys(), non_empty_comm_sig_paths.values(),
            width=0.85, color='forestgreen')
    plt.xlabel('Community')
    plt.ylabel('# of FDRs <=0.05')
    plt.title(figname+ ' significant pathways' \
              +' (total: '+str(total_sig_paths)+')' )
    plt.savefig(figname+'_overview.pdf')
    plt.close(fig)
    os.chdir('../')
    return 


def get_lowest_fdr(file, reactome_dir='reactome_analysis'):
    """
    check FDRs, pathwayID and community number from a reactome result file.
    """
    os.chdir(reactome_dir)
    xls = pd.ExcelFile(file)
    sheets = xls.sheet_names
    df = pd.read_excel(file, index_col=0, header=0, sheet_name=None)
    
    lowest_fdr = 0.05
    pathwayID = ''
    lowest_fdr_comm = '0'
    for comm in sheets:
        if lowest_fdr >= df[comm]['Entities FDR'].min():
            lowest_fdr = df[comm]['Entities FDR'].min()
            pathwayID =  df[comm].index[0]
            lowest_fdr_comm = comm
    os.chdir('../')
    
    return lowest_fdr, pathwayID, lowest_fdr_comm


if __name__ == '__main__':
    os.chdir("./Bcell_results/network_data")
    
    df = pd.read_csv('M0.125_Qf0.85_Qc0.95_PDN_communities.csv', index_col=0)
    get_reactome(df)