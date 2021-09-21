import matplotlib
matplotlib.use('agg')
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 15:31:35 2021

@author: Shuyue Xue
"""
import numpy as np, pandas as pd, networkx as nx 
import os, itertools
from operator import itemgetter

from wolframclient.language import wl, wlexpr
from wolframclient.evaluation import WolframLanguageSession

import dir_conf
global top_fdr_thresh, lead_pathway, reactome_results, minimum_comm_size_for_plot, \
       HPCC 

cd = os.getcwd()
# dir_conf.initialize('SLV') # test code with SLV data

top_fdr_thresh = 0.05
lead_pathway = False                   # only plot the pathway with lowest fdr value
minimum_comm_size_for_plot = 8         # prune the small communities from the plot 
HPCC = False

def set_additional_globals():
    """Set up more global varaibles"""
    global reactome_results, genes_in_each_community, \
           major_comms, major_comms_genes, louvain_groups, whole_net, fp
    
    # get the whole reatome_results
    reactome_results = pd.ExcelFile(dir_conf.reactome_dir)
    
    # get louvain groups (and major communities, genes in major communities)
    genes_in_each_community = pd.read_csv(dir_conf.community_dir, header=0, 
                                          index_col=0)
    major_comms = [c for c in genes_in_each_community.columns 
           if len(genes_in_each_community[c].dropna()) >= minimum_comm_size_for_plot]
    
    major_comms_genes = genes_in_each_community[major_comms].stack().values
    louvain_groups = [genes_in_each_community[c].dropna().tolist() 
                     for c in major_comms]
    
    # create the overall network to be plot
    fp = pd.read_excel(dir_conf.network_data_dir, index_col=0)
    
    whole_net = nx.from_pandas_edgelist(fp, edge_attr=True)  # create the graph
    whole_net = whole_net.subgraph(major_comms_genes)   # keep the major comms genes

    # set the community labels to the nodes
    nodes_comm_labels = {node: int( c.replace('C', '') ) 
                   for c in genes_in_each_community.columns  
                   for node in genes_in_each_community[c]}          
    nx.set_node_attributes(whole_net, nodes_comm_labels, name='comm')
    return


def set_up_whole_net_for_mathematica():
    """Helper: set up adjacency matrix for mathematica """
    global adjacency

    adjacency = nx.to_pandas_adjacency(whole_net)
    np.fill_diagonal(adjacency.values, np.Infinity)
    adjacency.replace(0, np.Infinity, inplace=True)
    adjacency = adjacency.values.tolist()
    return


def format_pylist_to_mathematica_legend(labelname_syntax, pylabel_list):
    """ 
    Helper funciton. convert python list to be mathematica list for the legend
    
    Parameters
    -------------
    labelname_syntax : str
        label type syntax in the Wolfram Mathematica .
    pylabel_list : list
        the list of names in python 
    
    Returns
    -------------
    labels : str
        Mathematica syntax 
    """
    labels ='%s -> %s' %(labelname_syntax, pylabel_list) 
    labels = labels.replace("'", "")
    labels = labels.replace("[", "{")
    labels = labels.replace("]", "}")
    
    return labels


def get_hit_genes(community_num):
    """ 
    Helper funciton. Return the hit genes + full pathways name:
    
    Parameters
    -------------
    community_num : int
        the community number.
    
    Returns
    -------------
    gene_list : list
        hit genes (labeled with their gene IDs, not the full isofrom name)
    pathways : list
        a list of significant pathways with FDR values below the fdr_threshold 
    """
    
    pathway_file = pd.read_excel(reactome_results, sheet_name = community_num)
    
    
    top_fdrs = pathway_file[ pathway_file['Entities FDR'] <= top_fdr_thresh ]
    pathways = top_fdrs['Pathway name'].tolist()        # full name of the path 
    # pathways = top_fdrs['Pathway identifier'].values  # short name of the path
    
    # Note: some entries are a whole string, need to split into a list gene IDs 
    genes_in_sign_pathways = top_fdrs['Submitted entities found'].apply(lambda x: x.split(';') )
    
    gene_list = list( set(itertools.chain(*genes_in_sign_pathways)) )
    
    if not pathways:
        print('No significant pathways in community %s' %community_num)
        
    elif lead_pathway:    # To get the lead pathway
        gene_list = genes_in_sign_pathways[0]
        pathways = pathways[0]
    
    return gene_list, pathways


# ------------- Start plotting -------------
set_additional_globals()
os.chdir(dir_conf.network_plot_dir)    # change the directory to save plots 
set_up_whole_net_for_mathematica()

# Start mathematica session
if HPCC:
    """ 
    If use on HPCC : 
    path = '/opt/software/Mathematica/<version>/Executables/MathKernel' 
    where <version> = 11.0.1, 11.3.0, 12.0.0, 12.1.0 
    """
    session = WolframLanguageSession('/opt/software/Mathematica/12.0.0/Executables/MathKernel')

else:
    session = WolframLanguageSession(
'C:\\Program Files (x86)\\Wolfram Research\\Mathematica\\12.0\\WolframKernel.exe')


#%% plot the major components in the same plot
minimum_component_size = 8
with session as wl_session:
     # wl_session.evaluate(wl.Needs("PlotLegends`"))
     g = wl_session.evaluate( wl.WeightedAdjacencyGraph( 
                            list(whole_net.nodes), 
                            adjacency, 
                            ))
     
     # find all the connected components 
     connected_g = wl_session.evaluate( wl.ConnectedComponents(g) )
     
     # Components to be plot: when nodes number >= minimum_component_size 
     components = [c for c in connected_g if len(c) >= minimum_component_size ]
     major_components = wl_session.evaluate( wl.Subgraph(g, components) )
          
     # group nodes by their louvan communities
     flatten = itertools.chain.from_iterable
     all_nodes = list(flatten(components))
     all_nodes = nx.get_node_attributes(whole_net.subgraph(all_nodes), 'comm')

     louvain_groups = { k: list(map(itemgetter(0), v))
           for k, v in itertools.groupby(sorted(all_nodes.items(), 
                          key=itemgetter(1)), itemgetter(1))
          }
    
     major_comm_labels = set(louvain_groups.keys())
     major_comm_labels = ["C%s" %i for i in major_comm_labels]

     legend_labels = format_pylist_to_mathematica_legend('PlotLegends',
                                                            major_comm_labels)
     major_comm_labels = format_pylist_to_mathematica_legend('CommunityLabels',
                                                            major_comm_labels)

     node_groups = list(louvain_groups.values())
     
     opac = wl_session.evaluate(wl.Opacity(0.86))
     plot_DN = wl_session.evaluate( wl.CommunityGraphPlot(
                              major_components, 
                              node_groups,      # grouped by communities
                              # wlexpr(major_comm_labels),
                              # wlexpr('EdgeStyle -> LightGray'),
                              wlexpr(legend_labels),
                              wlexpr('LabelStyle -> Directive[FontSize -> 3]'),
                              wlexpr('PlotStyle -> %r'%opac),
                              wlexpr('LegendSize -> 0.41'),
                              wlexpr('LegendMarkerSize -> 4'),
                             )
        )
    
                            
     # export the DN plots 
     wl_session.evaluate( wl.Export("differential_network.svg", plot_DN, "SVG"))
     
     
     # plot individual significant communities
     for comm in dir_conf.significant_comms:
         # get the node list 
         genes_in_this_comm = genes_in_each_community[
                              genes_in_each_community.columns[comm]
                                                 ].dropna().tolist()
        
         this_comm = wl_session.evaluate( wl.Subgraph(g, genes_in_this_comm) )
        
         # get highlights, labels and pathway names
         gene_labels, paths = get_hit_genes(comm)
        
         highlights = [ gene 
               for label in gene_labels for gene in genes_in_this_comm 
               if label in gene ]
        
         plot_cm = wl_session.evaluate( 
                 wl.HighlightGraph(
                     wl.Graph(this_comm, wlexpr('EdgeStyle -> LightGray'),
                            # wlexpr(hightlight_labels), # Label vertex
                            # wlexpr('VertexLabelStyle ->'+\
                            #        'Directive[Red, Italic, 31]'),
                           wlexpr('PlotStyle -> %r' %opac),
                         ),

                 wl.Style(wl.Subgraph(this_comm, highlights))
                 ),
         )
        
         # export community plots
         wl_session.evaluate( wl.Export("community %s.svg" %comm, plot_cm, "SVG") )

os.chdir(cd)