# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 13:41:31 2022

* Find hubs for each community

@author: Shuyue Xue (xueshuy1@msu.edu)
"""
import os
import networkx as nx, pandas as pd
cd = os.path.dirname(__file__)
os.chdir(cd)


# Setup
dataset = 'Bcell'

if 'SLV' == dataset:
    get_gene_ids = lambda gene_s: gene_s.split('_')[0]
if 'Bcell' == dataset:
    get_gene_ids = lambda gene_s: gene_s.split(':')[0] 

df_p = pd.read_excel('./Results/%s_results/network_data/%s_final_dn.xlsx'
                     %(dataset, dataset), index_col=0)
df_p.source = df_p.source.apply(get_gene_ids)
df_p.target = df_p.target.apply(get_gene_ids)
net_p = nx.from_pandas_edgelist(df_p, edge_attr='weight')
df_comms = pd.read_csv('./Results/%s_results/network_data/M0.125_Qf0.75_Qc0.995_DN_communities.csv'%dataset, 
                        index_col=0)
p_clusters = []
for c in df_comms.columns:
    comm = df_comms[c].dropna().drop_duplicates()
    if len(comm) < 8: 
        df_comms.drop(columns=c, inplace=True)
        continue
    p_clusters.append( comm.tolist() )
net_p = net_p.subgraph(sum(p_clusters, []))
rec_table = pd.ExcelFile('./Results/%s_results/reactome_analysis/M0.125_Qf0.75_Qc0.995_DN_community_reactome.xlsx' %dataset)


# Community Hubs 
def compute_community_degr_cent(net, community_number):
    """Compute community dcs for each gene."""
    nodes = p_clusters[ community_number ]
    subnet = net.subgraph(nodes)
    comm_dc_dict = nx.degree_centrality(subnet)
    comm_dc = sorted(comm_dc_dict.items(), reverse=True,
                          key=lambda x: x[1])
    f = lambda x: (x[0], format(x[1], '.4f'))
    comm_dc = list(map(f, comm_dc))
    return comm_dc


def add_dc(cell):
    mh13 = lambda gene: (gene, comm_degcen[gene][:5])
    cell = list(map(mh13, cell.split(';')))
    cell_sorted = sorted(cell, reverse=True, key=lambda x: x[1])
    cell_sorted = ', '.join([ '%s(%s)' %(g, dc) for (g, dc) in cell_sorted])
    return cell_sorted


writer  = pd.ExcelWriter("%s_Tables.xlsx"%dataset, engine='xlsxwriter')
hubs_dict = {}
for comm_number in range(len(df_comms.columns)):
    comm_degcen = compute_community_degr_cent(net_p, comm_number) 
    hubs_dict[comm_number] = comm_degcen
    comm_degcen = dict(comm_degcen)
    
    this_rect = pd.read_excel(rec_table, sheet_name = comm_number)    

    this_rect['Submitted entities found'] = \
    this_rect['Submitted entities found'].apply(add_dc)
    this_rect.to_excel(writer, sheet_name = df_comms.columns[comm_number], 
                   engine='xlsxwriter')
dc_table = pd.DataFrame.from_dict(hubs_dict, orient='index').transpose()
dc_table.columns = ['C%s'%i for i in dc_table.columns]
dc_table.to_excel(writer, sheet_name = 'Degree Centraliy', 
                  engine='xlsxwriter')
writer.save()
writer.close()