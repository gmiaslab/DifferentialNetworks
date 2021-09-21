import matplotlib
matplotlib.use('agg')
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 4 11:31:53 2020
This script:
    * Produces Bcell DN 
    * Performs reactome analysis 
    * Plots Heatmaps ----> must initialize before import the module 
    * Plots the DN ----> must initialize before import the module 
@author: Shuyue Xue
"""
import os, warnings
import networkx as nx
import preprocess_raw_data as prd, differential_network as dn, \
    reactome_results_evaluation as rre
import dir_conf

dir_conf.initialize('Bcell') # Initialize files and paths for heatmaps/DN plots
warnings.filterwarnings('ignore')

def pre_process(df, NaN_in_sparse_genes):
    # df = df.iloc[:13131]   # FOR CODE TESTING
    df = prd.clean_raw_data(df)
    df = prd.preprocess_data(df, missing_value_threshold=NaN_in_sparse_genes)
    return df

def get_high_fold_change_gene_pool(du, dt, dru, drt, q_fold, fold_dist_ttl='',
                                   dist=False,):
    gene_pool = prd.get_pairwise_fold_change_genes( dt, du, q_fold= q_fold,
            distribution=dist, title=fold_dist_ttl)
    
    gene_pool_repeated = prd.get_pairwise_fold_change_genes(
            drt, dru, q_fold= q_fold,
            distribution=dist, title=fold_dist_ttl+' (R)')
    
    gene_pool_list = [gene_pool, gene_pool,
                      gene_pool_repeated, gene_pool_repeated]
    return gene_pool_list

def make_network(df, target_gene_pool, q):
    df = df.loc[target_gene_pool]
    corr_mtrx = dn.make_correlation_matrix(df, quantile=q, export=False)
    G = dn.make_correlation_network(corr_mtrx)
    return G


def make_differential_networks(G1, G2,
                              g1_min_g2_file_name, g2_min_g1_file_name):

    g1, g2 = dn.make_differential_networks(G1, G2,
                              g1_min_g2_file_name = g1_min_g2_file_name,
                              g2_min_g1_file_name = g2_min_g1_file_name )
    return g1, g2


def pipeline(missing, q_fold, q_edge):
    Code_dir_path = os.getcwd()
    
    # load Bcell dataset: u=untreated, t=treated, ru(ru)= repeated untreated(treated)
    os.chdir(dir_conf.data_path)
    dfs_list = prd.configure_Bcell_datasets()
    os.chdir(dir_conf.results_dir)    # Working in Bcell_results dir
    
    diff_name_list = ['U_min_T', 'T_min_U', 'RU_min_RT', 'RT_min_RU']

    # 1. preprocess
    dfs_list = [pre_process(df, missing) for df in dfs_list]
    du, dt, dru, drt = dfs_list
    
    # 2. get target gene pool
    fold_dist_ttl = 'Missing' + str(missing) + ' Qf'+ str(q_fold)
    gene_pool_list = get_high_fold_change_gene_pool(du, dt, dru, drt, q_fold,
                                                    fold_dist_ttl, dist=False,)
    # 3. build corr networks
    U, T, RU, RT =[make_network(df, gene_pool, q=q_edge)
                    for (df, gene_pool) in zip(dfs_list, gene_pool_list)]
    
    # 4. build DN
    n, p = make_differential_networks(U, T,
                                      diff_name_list[0], diff_name_list[1])
    rn, rp = make_differential_networks(RU, RT,
                                        diff_name_list[2], diff_name_list[3])
    
    # find intersection of the two +DN (note: -DN, fn, will be skipped)
    # p = positive, n = negative in DNs
    fp = dn.get_DN_intersection(p, rp, graph_name='Bcell Final +DN')
    if 0 == len(fp):
        print('No intersection found')
        return 0
    
    # 5. export the final positive Differentail Net data
    fp_df = nx.to_pandas_edgelist(fp)
    fp_df_export = os.path.join('network_data', 
                                dir_conf.network_filename)
    fp_df.to_excel(fp_df_export)

    # 6. pass into the enrichment analysis : if find good FDR, save the community
    # So the communities where there's no good FDR found will not be saved 
    dn_name = 'M' + str(round(missing, 3)) + '_Qf'+ str(q_fold) + '_Qc'+str(q_edge) + '_DN'
    
    # export genes in each community:
    comm_df = dn.export_communities(fp, 'B', filename=dn_name)  
    
    
    # Simple overview of the network info    
    print("Network info:")
    print("Total Nodes = ", len(fp.nodes))
    print("Total Edges = ", len(fp.edges))
    print("Total communities =", len(comm_df.columns) )
    
    # export reactome analysis + count the total sig_pathways in this +DN
    # Note: different parameters give different DN
    total_sig_paths, community_sig_paths, non_empty_comm_sig_paths, max_comms = \
        rre.get_reactome(comm_df, dn_name)

    rre.plot_ntwk_enrichment_overview(non_empty_comm_sig_paths,
                                      total_sig_paths, 
                                  figname=dn_name)
    
    os.chdir(Code_dir_path)  # Back to the Code directory
    
    return comm_df


if __name__ == '__main__':
    # ------------- parameters ---------------------------------------
    NaNs = 1/8          # Missing Values in sparse genes
    q_fold = 0.75       # Quantile in fold change distribution:
    q_edge =  0.995     # Quantile in edge distribution
    #-----------------------------------------------------------------

    ntwk_df = pipeline(NaNs, q_fold, q_edge)
        
    # plot heatmaps
    import heatmaps
    heatmaps
        
    # plot DN network and communities : need Mathematica
    import network_plot_by_mathematica
    network_plot_by_mathematica