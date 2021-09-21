# -*- coding: utf-8 -*-
"""
Created on Fri March 13 13:31:41 2019

This module:
    * Creates the correlation matrix using the correlation cutoff in the
        filtered dataset.
    * Builds the coexpression network.
    * Computes the differential networks under two states.
      and returns both the postive and negtive differnetial networks.
    * Plots differential networks.
    * Exports the differential network structures to excel and JSON (D3).
    * Exports gene sets in each Louvain communities.

@author: Shuyue Xue
"""
import os, gc, time, warnings, json
from sknetwork.clustering import Louvain
import pandas as pd, numpy as np, networkx as nx, matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize
import community
from networkx.readwrite import json_graph
warnings.filterwarnings('ignore')


format_gene_name_SLV = lambda m: m.split('_')[0]
format_gene_name_Bcell = lambda m: m.split(':')[0]

def get_correlation_matrix(df): # Helper function
    """ Calculate correlation matrix and store the upper triangle only"""
    corr_df = df.corr()
    positions_to_fill_nans = pd.np.tril(corr_df, k=0).astype(bool)
    corr_df.mask(positions_to_fill_nans, inplace=True)
    return corr_df


def make_correlation_matrix(df_clean, quantile = 0.99,
                           export=False, file_name='',):
    """
    Calculate the correlation matrix above the corrleaiton cutoff 
    
    parameters
    -------------   
    df_clean : data frame
        The clean data frame results from the preprocess_data() method.
    correlation_cutoff : float
        correlation value results from the calculate_correlation_cutoff().

    
    Returns
    -------------
    corr_df : data frame
    """

    # print('\nStart calculating correlation matrix ....')
    # print('Starting dataframe has the shape: ', df_clean.shape)
    corr_df = get_correlation_matrix(df_clean.T)
    corr_df = corr_df.stack()
    
    # Exclude all 1's from the origianl correlation distribution
    exclude_1s = corr_df[
            corr_df.apply(np.floor) !=1.0000 ]

    # Find the 95% quantile at the distribution of the exclude_1s-distribution
    correlation_at_quantile = exclude_1s.quantile( q=quantile )
    # print('correlation at %s quantile : ' %quantile, correlation_at_quantile)
    # print('End calculation: cutoff ----')

    # Select the correlations above the cutoff
    corr_df = corr_df[ corr_df >= correlation_at_quantile ]

    gc.collect()
    # print('Size of the correlation df [GB] ): ',
    #       corr_df.memory_usage(deep=True)*10**(-9) )

    if export==True:
        corr_df.to_pickle(file_name+'_corr.pkl')
        # corr_df.unstack().to_excel(file_name+'_corr.xlsx')

    # print('End calculation: correlation matrix for the network ----')

    return corr_df


def make_correlation_network(correlation_df):
    """
    Build correlation network.

    Parameters
    -------------
    correlation_df : data frame
        The stacked correlation matrix by the make_corrlation_matrix().
    
    Returns
    -------------
    G : graph
    The gene coexpression network.  
    """
    # print('\nStart building correlation network ....')
    tic = time.time()
    G = correlation_df

    G.index.names = ['n1', 'n2']
    G = G.reset_index()
    G.columns = ['source', 'target', 'weight']
    G = nx.from_pandas_edgelist(G, 'source', 'target', edge_attr=True)
    toc = time.time()
    # print('The corr network has %s Nodes and %s edges'
    #       %(len(G.nodes), len(G.edges)) )
    # print('End building correlation network ----\nTime taken= %s sec'
    #       % round(toc-tic) )

    del correlation_df
    gc.collect()

    return G


def make_differential_networks(g1, g2, export=False,
                         g1_min_g2_file_name='',
                         g2_min_g1_file_name='', 
                         subgroup_method='skLouvain'):
    # print('\nStart making differential networks ....')
    # tic = time.time()

    temp1 = g1.copy()
    temp2 = g2.copy()

    # In g1: leave edges in g1 but not in g2
    g1.remove_edges_from( ed for ed in temp2.edges() )
    isolated_nodes = list(nx.isolates(g1))
    # print('# of single nodes =', len(isolated_nodes))
    g1.remove_nodes_from(isolated_nodes)

    # In g2: leave edges in g2 but not in g1
    g2.remove_edges_from( ed for ed in temp1.edges() )
    isolated_nodes = list(nx.isolates(g2))
    # print('# of single nodes =', len(isolated_nodes))
    g2.remove_nodes_from(isolated_nodes)

    del temp1, temp2
    gc.collect()

    # Set louvain communities
    set_subgroup(g1, method=subgroup_method)
    set_subgroup(g2, method=subgroup_method)
    
    
    # toc = time.time()
    # print('End building correlation network ----\nTime taken= %s sec'
    #       % round(toc-tic) )
    return g1, g2


def get_DN_intersection(p, rp, graph_name='final DN', 
                        subgroup_method='skLouvain'):
    # print('\nStart making intersection DN ....')
    # tic = time.time()

    # find intesection nodes:
    common_nodes = set(p.nodes).intersection(set(rp.nodes))

    # create subgraphs
    temp1 = p.subgraph(common_nodes).copy()
    temp2 = rp.subgraph(common_nodes).copy()
    # print(len(temp1.edges), len(temp2.edges))

    # create the graph that contains only the edges existing in both.
    final_p = nx.intersection(temp1, temp2)
    # isolated_nodes = list(nx.isolates(final_p))
    # print('# of single nodes =', len(isolated_nodes))
    # final_p.remove_nodes_from(isolated_nodes)
    # temp1 = temp1.subgraph(final_p.nodes).copy()
    # temp2 = temp2.subgraph(final_p.nodes).copy()
    if 0 == len( final_p.edges() ):
        print("Intersection DN doesn't exist")
        toc = time.time()
        # print('End building final DN ----\nTime taken= %s sec'
        #   % round(toc-tic) )
        return final_p.edge_subgraph(final_p.edges())

    # clean up isolated nodes
    temp1 = temp1.edge_subgraph( final_p.edges() )
    temp2 = temp2.edge_subgraph( final_p.edges() )
    
    # set edge weight to be the average of the 2 DNs
    d1 = nx.to_pandas_edgelist(temp1)
    d2 = nx.to_pandas_edgelist(temp2)

    # print('Match? ',
    #       ( np.sort(d1['source']) == np.sort(d2['source']) ).all(),
    #       ( np.sort(d1['target']) == np.sort(d2['target']) ).all()
    #       )
    """find the intersection of the 2 +DNs."""

    d1['weight'] = 0.5*(d1['weight']+d2['weight'])
    final_p = nx.from_pandas_edgelist(d1, edge_attr='weight')


    del temp1, temp2
    gc.collect()

    # Set new communities: 
    # set_subgroup(final_p, method='pyLouvain')   # use the python-louvin
    set_subgroup(final_p, method=subgroup_method)

    # toc = time.time()
    print(graph_name+' has %s Nodes and %s edges'
          %(len(final_p.nodes), len(final_p.edges)) )
    # print('End building final DN ----\nTime taken= %s sec'
    #       % round(toc-tic) )

    return final_p


def export_gene_list(gene_list, B_or_S): 
    """Export general gene list and format the gene names 
    as in the Enrichment analysis"""
    
    if B_or_S == 'S':
        gene_names_formatted = gene_list.stack().apply(
                                format_gene_name_SLV).unstack()
    if B_or_S == 'B':
        gene_names_formatted = gene_list.stack().apply(
                                format_gene_name_Bcell).unstack()
    return gene_names_formatted


def export_communities(diff_network, B_or_S, filename=''): 
    """
    Export gene sets in each of the communities.

    Parameters
    -------------
    diff_network : networkx graph 
        One differential network from make_differential_networks().
    B_or_S: str,
        'B' for B cell data, 'S' for SLV dataset, 
        for the gene name formatting.
    filename: str
        file name for saving.  

    Returns
    -------------
    None.
    """
    # print('\nSaving the genes in each communities ....')
    tic = time.time()
    
    # Save into the results path
    if not filename:
        filename = 'communities'
    if not os.path.exists('./network_data'):
        os.makedirs('network_data')

    community_gene_sets = pd.DataFrame()
    nodes_df = pd.DataFrame.from_dict(
            dict(diff_network.nodes(data=True)), orient='index'
            )
    
    # combine communities as sheets in one dataframe
    for i in get_communities(diff_network):
        temp = nodes_df[ nodes_df['community']==i ].index.values
        temp = pd.DataFrame(temp, columns=['C'+str(i)])
        community_gene_sets = pd.concat([community_gene_sets, temp],
                                        axis=1)

    ############# ############# ############# #############
    # Modified 12/4/2020: output the orginial gene names before name formatting
    community_gene_sets.to_csv('./network_data/'+ \
                                  filename + \
                                  '_full_name' + \
                                  '_communities.csv')
    ############# ############# ############# #############
    
    # Format the gene names for Genontology analysis
    if B_or_S == 'S':
        community_gene_sets = community_gene_sets.stack().apply(
                                format_gene_name_SLV).unstack()
    if B_or_S == 'B':
        community_gene_sets = community_gene_sets.stack().apply(
                                format_gene_name_Bcell).unstack()


    community_gene_sets.to_csv('./network_data/'+ \
                                  filename +'_communities.csv')
    # community_gene_sets.to_excel('./network_data/'+ \
    #                              filename +'_communities.xlsx')
    # community_gene_sets.to_pickle('./network_data/'+ \
    #                              filename +'_communities.pkl')

    toc = time.time()
    # print('Finished Saving ----\nTime taken= %s sec' % round(toc-tic) )
    return community_gene_sets


def set_subgroup(graph, method='skLouvain'):
    """
    Set subgroups of the graph either by .

    Parameters
    -------------
    graph : networkx graph 
    method : str,
        skLouvain, pyLouvain, greedy, naive-greedy, girvan-newman, (k-means)?

    Returns
    -------------
    None.
    """
    if '' == method: 
        print('No method selected, we use skLouvain') 
        set_by_sk_louvain(graph)
    elif 'skLouvain' == method: set_by_sk_louvain(graph)
    elif 'pyLouvain' == method: set_by_python_louvain(graph)
    else: print('selection is out of the availabel options....') 
    return 

def set_by_python_louvain(G):    # the python-louvain version 
    partition = community.best_partition(G)
    return nx.set_node_attributes(G, partition, name='community')

def set_by_sk_louvain(graph): # set the community by sknetwork
    adjacency = nx.to_scipy_sparse_matrix(graph)

    """Check this new seed today"""
    seed = np.random.seed(seed=13)                 # Lock the random state in Louvain
    comm = Louvain(random_state=seed).fit_transform(adjacency)

    comm_dict = {node : comm 
                 for (node, comm) in zip (list(graph.nodes), comm)}
    nx.set_node_attributes(graph, comm_dict, name='community')
    return

def gene_id_map( graph_w_id , graph_wo_id):     # reset the gene labels 
    # create the gene id map
    mapping = {node: gene_id for (node, gene_id) in 
               zip( list(graph_w_id.nodes),  list(graph_wo_id.nodes))}
    
    # reset nodes labels in the graph_wo_id
    nx.relabel_nodes(graph_wo_id, mapping, copy=False)    
    return 

def get_communities(G):
    return set( nx.get_node_attributes(G, 'community').values() )


def make_JSON(graph, filename, into_type=''):
    "Making jason files for either cytoscape or D3 plot"
    if into_type == 'cyto': 
        jsonData =  json_graph.cytoscape_data(graph)
        filename = filename,+"_JSON_cyto"
    else:
        jsonData = json_graph.node_link_data(graph)
        filename = filename+"_JSON_D3"
    with open('./JSON_files'+'/'+ filename+".json", 'w') as write_file:
            json.dump(jsonData, write_file, indent=4)
    return jsonData


def plot(graph, legend=True, labels=False, 
         major=True, filename='', save=True):
    """If major=True, plot the major component, else plot the entire graph"""

    # print('Start plotting the network ....')
    tic = time.time()

    plt.ioff()
    plt.clf()
    
    if major:
        graph = major_comp(graph)
        # set_subgroup(graph, method='skLouvain')    # Recalculate communities
        filename = filename+'_major_component'

    random_state = np.random.RandomState(13)
    pos=nx.spring_layout(graph, seed=random_state)
    # pos=nx.spring_layout(graph)

    edges, weights = \
    zip(*nx.get_edge_attributes(graph,'weight').items())
    nx.draw_networkx_edges(graph, pos,
                        edge_color = weights,
                        edge_cmap=plt.cm.PRGn,
                        width = 0.017,
                        alpha=0.1)
    
    communities = list(get_communities(graph))
    # colormap = get_cmap('viridis')
    # colormap = get_cmap('hsv')
    colormap = get_cmap('tab20c')
    
    my_norm = Normalize( vmin=min(communities), vmax= max(communities) )
    
    communities_color = {k: k for k in communities}
    # else:
    #     communities_color = {k: (4*k)%len(communities) for k in communities}
        
    if len(communities) >=20:
        communities_color = {k: (4*k)%len(communities) for k in communities}
    
    # colormap = mpl.colors.ListedColormap(['red', 'green', 'm', \
    #                                       'blue', 'cyan', \
    #                                       'k', 'y', 'purple'],
    #                                      N = len(communities))

    node_color_list = [ communities_color[graph.nodes[node]['community']]
                        for node in graph.nodes ]
    node_color_list = my_norm(node_color_list)
    
    nx.draw_networkx_nodes(graph, pos,
                        node_size=0.68,                      
                        cmap = colormap,
                        # norm = my_norm,
                        node_color = node_color_list,
                        alpha = 0.68)
    
    if labels == True: 
        node_labels_list = dict( zip (list(graph.nodes), graph.nodes) )
        nx.set_node_attributes(graph, node_labels_list, 'label')
        nx.draw_networkx_labels(graph, pos, labels=node_labels_list, font_size = 3)
           
    if legend == True:
        for comm in communities:

            plt.scatter([], [],
                        color=colormap( my_norm( communities_color[comm] ) ),
                        alpha=0.68, 
                        s=3.1,
                        
                        label=str(comm), )
        plt.legend(scatterpoints=1,  labelspacing=0.3,
                   frameon=False,
                   loc="best",
                   title='communities', fontsize=4,)

    plt.axis('off')
    if save == True:
        if not os.path.exists('./network_plots'):
            os.makedirs('network_plots')
        plt.savefig('./network_plots'+'/'+filename+".pdf", format="PDF")
    plt.close()

    toc = time.time()
    # print('Finished plotting ----\nTime taken= %s sec' % round(toc-tic) )
    return graph

def major_comp(G):
    "A helper for the plot() function."
    num_components = nx.number_connected_components(G)
    print('This graph has %s components.'\
          ' Here we plot the major component.' %num_components)
    majr_comp = max(nx.connected_components(G), key=len)
    H = G.subgraph(majr_comp).copy()
    return H



if __name__ == '__main__':
    os.chdir('SLV_results')

    file1 = "H1_corr.pkl"
    file2 = "H2_corr.pkl"
    corr_df1 = pd.read_pickle(file1).iloc[:13131]
    corr_df2 = pd.read_pickle(file2).iloc[:13131]

    for folder in ["network_plots", 'network_data']:
        if not os.path.exists(folder):
            os.makedirs(folder)

    G1 = make_correlation_network(corr_df1)
    G2 = make_correlation_network(corr_df2)

    save_g1_file ='h1_min_h2'
    save_g2_file ='h2_min_h1'
    g1, g2 = make_differential_networks(G1, G2,
                              g1_min_g2_file_name=save_g1_file,
                              g2_min_g1_file_name=save_g2_file )

    plot(g1, filename=save_g1_file)
    plot(g2, filename=save_g2_file)

    make_JSON(g1, filename=save_g1_file)
    make_JSON(g2, filename=save_g2_file)
    export_communities(g1, 'g1')