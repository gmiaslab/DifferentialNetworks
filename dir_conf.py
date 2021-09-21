# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 15:31:41 2021
This script initializes file names and directory paths 
for plotting heatmaps and the netowrks. 
    
@author: Shuyue Xue
"""
import os

cd = os.getcwd()

def initialize(data):
    """
    Configure data paths and some file names 
    
    parameters
    -------------   
    data: string
        'SLV' or 'Bcell'
    """
    global data_set, data_path, results_dir, network_filename, significant_comms, \
           reactome_filename, reactome_dir, community_filename, community_dir, \
           network_data_dir, network_plot_dir

    data_set = data    
    
    # create a results directory 
    results_dir_name = data_set + '_results'
    if not os.path.exists(results_dir_name): os.makedirs(results_dir_name)
    os.chdir(results_dir_name)    # Working dir in xx_results dir
    for folder in ["network_plots", 'network_data']:
        if not os.path.exists(folder):
            os.makedirs(folder)
    
    # set the community number for individual community plots
    if data_set == 'SLV':
        significant_comms = [0, 1, 2, 4, 8, 9] 
    if data_set == 'Bcell':
        significant_comms = [2, 4, 5, 6, 7, 9, 10, 12, 13]

    # set the file names and dir paths 
    data_path = os.path.join(cd, data_set+'_data')
    results_dir = os.path.join(cd, results_dir_name)
    
    network_filename = data_set+'_final_dn.xlsx'
    network_data_dir = os.path.join(results_dir,'network_data', network_filename)
    network_plot_dir = os.path.join(results_dir,'network_plots')
    
    reactome_filename = 'M0.125_Qf0.75_Qc0.995_DN_community_reactome.xlsx'
    reactome_dir = os.path.join(results_dir, 'reactome_analysis',
                                     reactome_filename )
    
    community_filename = 'M0.125_Qf0.75_Qc0.995_DN_full_name_communities.csv'
    community_dir = os.path.join(results_dir, 'network_data',
                                     community_filename)
    
    os.chdir(cd)