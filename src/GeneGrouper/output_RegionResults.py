'''
Arrange outputs so that they have all metadata combined and are ready to be plotted 
'''


import os
from os.path import join as pjoin
import re
from Bio.Seq import Seq
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
import warnings
import time
import sqlite3 as sql
from collections import defaultdict
import gc
import sys
import scipy.cluster.hierarchy as shc
from gtr_utils import change_ToWorkingDirectory, make_OutputDirectory, merge_ManyFiles, multiprocesssing_Submission, generate_AssemblyId, chunks


def align_SeedGenes(df_sr):
	'''
	Use global_strand to determine the relative ordering of each region
	Return dataframe df_sr with corrected gene_position called new_pos
	'''
	## arrange genes from start to end, normalize order to strand=1
	# first get the strand orientation of the seed_gene
	df_strand = df_sr[df_sr['seed_region_id']==df_sr['locus_tag']]
	df_strand['global_strand'] = df_strand['strand']
	df_strand = df_strand[['seed_region_id','global_strand']]
	df_sr = df_sr.merge(df_strand,on='seed_region_id')
	df_sr['global_strand'] = df_sr['global_strand'].astype(int)

	# for each locus tag in a gene region, assign a gene position based on whether the global strand is 1 or -1
	df_sr['gene_position'] = df_sr[df_sr['global_strand'] == 1].groupby('seed_region_id').cumcount()+1
	df_sr['gene_position_1'] = df_sr[df_sr['global_strand'] == -1].groupby('seed_region_id').cumcount(ascending=False)+1
	df_sr['gene_position'].fillna(df_sr['gene_position_1'], inplace=True)
	df_sr = df_sr.drop(columns='gene_position_1')

	# add column that tags whether a locus_tag is the seed gene for that region
	df_sr['is_seed'] = np.where(df_sr['locus_tag']==df_sr['seed_region_id'], 1, 0)

	# extract just the positions of each region seed 
	max_seed_position = df_sr[df_sr['is_seed']==1][['seed_region_id','gene_position','is_seed']]
	max_seed_position['relative_seed_position'] = max_seed_position['gene_position']
	max_seed_position = max_seed_position.drop(columns=['gene_position','is_seed'])
	# merge with df_sr
	df_sr = df_sr.merge(max_seed_position,on='seed_region_id')
	# re-calculate seed positions relative to the global_anchor_position
	df_sr['new_pos'] = df_sr['gene_position'] - df_sr['relative_seed_position']
	return(df_sr)


def subset_SimilarityMatrix(rep_list):
	'''
	'''
	df_jac = pd.read_csv(pjoin('internal_data','seed_regions_jaccard_dist.csv'))
	df_jac = df_jac.rename(columns={df_jac.columns[0]:'seed_region_id'},inplace=False)
	df_jac = df_jac[df_jac['seed_region_id'].isin(rep_list)]
	df_jac = df_jac.drop(df_jac.columns.difference(rep_list), axis=1)
	return(df_jac)



def melt_LabelDissimilarty(df_cr, df_jac):
	'''
	'''
	df_crm = df_cr[['dbscan_label','full_index_pos','mean_dissimilarity','representative_relative_dissimilarity']]
	df_crm = pd.melt(df_crm,id_vars=['dbscan_label','full_index_pos'])
	df_crm = df_crm.rename(columns={'variable':'measurement', 'value':'dissimilarity'})
	df_crm = df_crm.merge(df_jac, on ='dbscan_label')
	df_crm.to_csv(pjoin('internal_data','rtable_label_dissimilarity.csv'),index=False)


def melt_IdentityCoverage(df_sr):
	'''
	'''
	df_id = df_sr[['region_id','pident','qcovs','dbscan_label','label_representative','full_index_pos','order']].drop_duplicates()
	df_id['representative_region'] = np.where(df_id['label_representative']==df_id['full_index_pos'], 'representative','member')
	df_id = df_id[['region_id','pident','qcovs','dbscan_label','representative_region','order']]
	df_id = pd.melt(df_id, id_vars = ['region_id','dbscan_label','order','representative_region'])
	df_id.to_csv(pjoin('internal_data','rtable_label_identitycoverage.csv'),index=False)


def sort_DataFrameByHierarchicalClustering(df_input):
	'''
	Sort dataframe by hierarchical clustering. Same order as visualization outputs
	df_horder contains the hierarchical ordering output. 
	'''
	df_horder = pd.read_csv(pjoin(os.getcwd(),'internal_data','full_region.csv'))[['dbscan_label','order']].drop_duplicates()
	df_horder['order'] = np.where(df_horder['dbscan_label']==-1,-1,df_horder['order'])
	df_input = df_input.merge(df_horder,on='dbscan_label')
	df_input = df_input.sort_values(by='order', ascending=False)
	return(df_input)

def rename_ColumnsForUser(df):
	'''
	'''
	rename_dict = {'dbscan_label':'group','order':'hierarchical_order','region_id':'member_id','pseudo_check':'is_pseudogene','pident':'identity_to_query','qcovs':'coverage_to_query','evalue':'evalue_to_query',
	'mean_dissimilarity':'pairwise_mean_dissimilarity','representative_relative_dissimilarity':'relative_dissimilarity_to_representative','relative_seed_position':'relative_gene_position'}
	df = df.rename(columns=rename_dict)
	return(df)

## --------------Start script---------------------##

def OutputRegionResults(
	UserInput_main_dir,
	UserInput_output_dir_name,
	):
	'''
	'''

### Troubleshooting inputs ###
# UserInput_main_dir = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/syntenease_project/gtr/testbed/dataset4/test1'
# UserInput_output_dir_name = 'pdua'
# pd.set_option('display.max_columns', None)

###--------------------###


	change_ToWorkingDirectory(directory_name = pjoin(UserInput_main_dir,UserInput_output_dir_name))
	conn = sql.connect('seed_results.db')
	conn2 = sql.connect(pjoin(UserInput_main_dir,'genomes.db'))

	## prepare taxonomy data ##
	df_meta = pd.read_sql_query("SELECT * FROM gb_metadata", conn2)
	df_meta['genus'] = [t.split(' ')[0].lower() for t in df_meta['organism'].tolist()]
	df_meta['species'] = [t.split(' ')[1].lower() for t in df_meta['organism'].tolist()]
	df_meta = df_meta[['assembly_id','genus','species','organism']]

	## add ortholog, mmseqs cluster ids, dbscan labels, and filtered hit results to the extracted seed regions ##
	df_sr = pd.read_sql_query("SELECT * FROM seed_regions", conn) 
	df_cr = pd.read_sql_query("SELECT * FROM dbscan_label_representatives", conn)
	df_or = pd.read_sql_query("SELECT * FROM select_mcl_results", conn)

	# ----  Make comprehensive region & cluster data table ---- #

	# merge dbscan labels with seed regions
	df_cr = df_cr[['seed_region_id','dbscan_label','mean_dissimilarity','representative_relative_dissimilarity','label_representative','full_index_pos']]
	df_sr = df_sr.merge(df_cr, left_on='region_id', right_on='seed_region_id')

	# merge ortholog clusters with seed regions. Requires doing a small recalculation as in cluster_Regions.py
	# this is becasue not all genes in seed regions have a sequence because they are listed as pseudo genes.
	df_sr = df_sr.merge(df_or, left_on = 'locus_tag',  right_on =  'orf_id', how = 'outer')
	df_sr['orf_id'] = df_sr['orf_id'].fillna(df_sr['locus_tag'])
	max_ortho_cluster_id_val = df_sr['ortho_cluster_id'].max(skipna = True)
	max_mmseqs_cluster_id_val = df_sr['mmseqs_cluster_id'].max(skipna = True)
	df_sr['ortho_cluster_id'] = df_sr['ortho_cluster_id'].fillna(max_ortho_cluster_id_val +1)
	df_sr['mmseqs_cluster_id'] = df_sr['mmseqs_cluster_id'].fillna(max_mmseqs_cluster_id_val +1)

	# add taxonomy data
	df_sr = df_sr.merge(df_meta,on='assembly_id')

	## big merge: gene, product names, for all seed_region members
	df_genomes = pd.read_sql_query("SELECT locus_tag, refseq_locus_tag, refseq_product, refseq_gene, cds_start, cds_end FROM gb_features WHERE locus_tag IN {}".format(tuple(df_sr.locus_tag.tolist())), conn2) 
	df_sr = df_sr.merge(df_genomes,on='locus_tag')
	df_sr = align_SeedGenes(df_sr=df_sr)

	## Use hierarchical clustering to order the output
	# filter for columns and rows that are representatives
	rep_list = df_cr[df_cr['full_index_pos']==df_cr['label_representative']]['seed_region_id'].tolist()
	label_list = df_cr[df_cr['full_index_pos']==df_cr['label_representative']]['dbscan_label'].tolist()
	df_jac = subset_SimilarityMatrix(rep_list=rep_list)

	# make dendrogram with column names as labels
	dend = shc.dendrogram(shc.linkage(y= df_jac, method='ward',
		optimal_ordering = True), labels = df_jac.columns)

	# make the dendrogram to a dataframe
	rep_label_map = pd.DataFrame({'representative_id':rep_list,'dbscan_label':label_list})
	df_reorder = pd.DataFrame({
		'representative_id': dend['ivl'],
		'order': range(0,len(dend['ivl']))})
	df_reorder = df_reorder.merge(rep_label_map, on = 'representative_id')
	df_sr = df_sr.merge(df_reorder.drop(columns='representative_id'), on='dbscan_label')


	# ---- Write comprehensive region & cluster data table for internal use ---- #
	df_sr.to_csv(pjoin('internal_data','full_region.csv'),index=False)

	# ---- Make metadata for plotting outputs ---- #
	melt_LabelDissimilarty(df_cr=df_cr, df_jac=df_reorder)
	melt_IdentityCoverage(df_sr=df_sr)
	df_meta.to_csv(pjoin('internal_data','rtable_genome_input_metadata.csv'),index=False)
	df_sr[['assembly_id','region_id','dbscan_label']].drop_duplicates().to_csv(pjoin('internal_data','rtable_summarized_region_clusters.csv'), index=False)

	# ---- Write summary tabular outputs for user ---- #

	## groups labels with dissimilarity, pident, qcovs info
	df_t1 = pd.read_csv(pjoin('internal_data','full_region.csv'))[['dbscan_label','region_id','mean_dissimilarity','pident','qcovs']].drop_duplicates().melt(id_vars=['dbscan_label','region_id'])
	max_result_vals = df_t1.groupby(['dbscan_label','variable']).max().reset_index().drop(columns='region_id')
	max_result_vals['variable'] = 'max_'+max_result_vals['variable']
	max_result_vals = max_result_vals.pivot_table(index='dbscan_label',columns='variable',values='value').reset_index()
	min_result_vals = df_t1.groupby(['dbscan_label','variable']).min().reset_index().drop(columns='region_id')
	min_result_vals['variable'] = 'min_'+min_result_vals['variable']
	min_result_vals = min_result_vals.pivot_table(index='dbscan_label',columns='variable',values='value').reset_index()
	mean_result_vals = df_t1.groupby(['dbscan_label','variable']).mean().reset_index()
	mean_result_vals['variable'] = 'mean_'+mean_result_vals['variable']
	mean_result_vals = mean_result_vals.pivot_table(index='dbscan_label',columns='variable',values='value').reset_index()
	sd_result_vals = df_t1.groupby(['dbscan_label','variable']).std().reset_index()
	sd_result_vals['variable'] = 'sd_'+sd_result_vals['variable']
	sd_result_vals = sd_result_vals.pivot_table(index='dbscan_label',columns='variable',values='value').reset_index()
	group_counts = pd.read_csv(pjoin('internal_data','full_region.csv'))[['dbscan_label','region_id']].drop_duplicates().groupby('dbscan_label').count().rename(columns={'region_id':'member_counts'}).reset_index()
	df_r1 = max_result_vals
	df_r1 = pd.concat([df_r1, min_result_vals.iloc[:,1:],mean_result_vals.iloc[:,1:],sd_result_vals.iloc[:,1:], group_counts.iloc[:,1]], axis = 1)
	df_r1 = df_r1[['dbscan_label','member_counts',
	'min_mean_dissimilarity','mean_mean_dissimilarity','max_mean_dissimilarity','sd_mean_dissimilarity',
	'min_pident','mean_pident','max_pident','sd_pident',
	'min_qcovs','mean_qcovs','max_qcovs','sd_qcovs']]
	df_r1 = sort_DataFrameByHierarchicalClustering(df_input=df_r1)
	df_r1 = rename_ColumnsForUser(df=df_r1)
	df_r1.to_csv('group_statistics_summmary.csv',index=None)


	## gene content for representative groups
	df_t1 = pd.read_csv(pjoin('internal_data','full_region.csv'))[['dbscan_label','region_id','assembly_id',
		'label_representative','full_index_pos','refseq_locus_tag','refseq_gene','refseq_product','gene_position','is_seed']]
	df_t1 = df_t1[df_t1['full_index_pos']==df_t1['label_representative']]
	df_t1 = df_t1[['dbscan_label','refseq_locus_tag','refseq_gene','refseq_product','gene_position','is_seed']]
	df_t1 = df_t1.replace(np.nan, 'NA', regex=True)
	df_t1['gene_info'] = df_t1['refseq_gene'] + ' | ' + df_t1['refseq_locus_tag'] + ' | ' + df_t1['refseq_product']
	df_t1['gene_info'] = np.where(df_t1['is_seed'] ==1, '* '+df_t1['gene_info'], df_t1['gene_info'])
	df_t1 = df_t1[['gene_position','dbscan_label','gene_info']]
	df_t1['dbscan_label'] = 'group_'+df_t1['dbscan_label'].astype(int).astype(str)
	df_t1 = df_t1[['gene_position','dbscan_label','gene_info']].pivot(index='gene_position',columns='dbscan_label',values='gene_info').reset_index()
	# reorder to hierarchical clustering
	df_horder = pd.read_csv(pjoin(os.getcwd(),'internal_data','full_region.csv'))[['dbscan_label','order']].drop_duplicates()
	df_horder['order'] = np.where(df_horder['dbscan_label']==-1,-1,df_horder['order'])
	df_horder = df_horder.sort_values(by='order', ascending=False)
	group_label_list = ['gene_position']
	for c in df_horder.iterrows():
		group_label_list.append('group_'+c[1][0].astype(int).astype(str))
	df_t1 = df_t1[group_label_list]
	df_t1.to_csv('representative_group_member_summary.csv',index=None)


	## group distribution by genus
	df_g1 = pd.read_csv(pjoin('internal_data','full_region.csv'))[['dbscan_label','assembly_id','genus']].drop_duplicates()
	df_g1 = df_g1.groupby(['dbscan_label','genus']).count().reset_index()
	df_g2 = pd.read_csv(pjoin('internal_data','full_region.csv'))[['dbscan_label','seed_region_id','genus']].drop_duplicates()
	df_g2 = df_g2.groupby(['dbscan_label','genus']).count().reset_index()
	df_g1 = df_g1.merge(df_g2, on = ['dbscan_label','genus'])
	df_t2 = df_meta[['assembly_id','genus']].groupby('genus').count().reset_index().rename(columns={'assembly_id':'genus_total'})
	df_g1 = df_g1.merge(df_t2, on='genus')
	df_g1['perc_present'] = ((df_g1['assembly_id']/df_g1['genus_total'])*100).round(2).astype(str)
	df_g1['perc_present'] = np.where(df_g1['assembly_id']<df_g1['seed_region_id'],df_g1['perc_present']+'*',df_g1['perc_present'])
	df_g1 = df_g1[['genus','perc_present','dbscan_label']]
	df_g1 = df_g1.pivot(index='dbscan_label',columns='genus',values='perc_present').reset_index()
	df_g1 = sort_DataFrameByHierarchicalClustering(df_input=df_g1)
	df_g1 = rename_ColumnsForUser(df=df_g1)
	df_g1.to_csv('group_taxa_summary.csv',index=None)


	## long form datasets for user
	df_l1 = pd.read_csv(pjoin('internal_data','full_region.csv'))
	df_l1 = df_l1[df_l1['label_representative']==df_l1['full_index_pos']][['dbscan_label','region_id']].drop_duplicates().rename(columns={'region_id':'group_representative'})
	df_l1 = df_l1.merge(pd.read_csv(pjoin('internal_data','full_region.csv')), on = 'dbscan_label')
	df_l1 = df_l1[['region_id','dbscan_label','group_representative',
		'assembly_id','contig_id','gene_position','orf_id','pseudo_check','refseq_locus_tag','refseq_gene','refseq_product',
		'is_seed','pident','qcovs','evalue',
		'ortho_cluster_id','ortho_cluster_members','mmseqs_cluster_id','mmseqs_cluster_representative',
		'mean_dissimilarity','representative_relative_dissimilarity',
		'organism','genus','species',
		'cds_start','cds_end','strand','relative_seed_position',
		'order']]
	df_l1 = rename_ColumnsForUser(df=df_l1)
	df_l1['group'] = df_l1['group'].astype(int)
	df_l1.to_csv('group_regions.csv',index=None)


	# close SQL connections
	conn.close()
	conn2.close()

	# clean up
	os.system('cp ./ortholog_clusters/merged_aa.faa ./group_region_seqs.faa')
	os.system('rm -R ./ortholog_clusters')
	os.system('rm -R ./regions')

	print('\nAll outputs files written to {}\n'.format(os.getcwd()))


## --------------End script---------------------##









