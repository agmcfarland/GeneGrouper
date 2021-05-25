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
	df_jac = pd.read_csv(pjoin('results','seed_regions_jaccard_dist.csv'))
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
	df_crm.to_csv(pjoin('results','rtable_label_dissimilarity.csv'),index=False)


def melt_IdentityCoverage(df_sr):
	'''
	'''
	df_id = df_sr[['region_id','pident','qcovs','dbscan_label','label_representative','full_index_pos','order']].drop_duplicates()
	df_id['representative_region'] = np.where(df_id['label_representative']==df_id['full_index_pos'], 'representative','member')
	df_id = df_id[['region_id','pident','qcovs','dbscan_label','representative_region','order']]
	df_id = pd.melt(df_id, id_vars = ['region_id','dbscan_label','order','representative_region'])
	df_id.to_csv(pjoin('results','rtable_label_identitycoverage.csv'),index=False)



## --------------Start script---------------------##

def OutputRegionResults(
	UserInput_main_dir,
	UserInput_output_dir_name,
	):
	'''
	'''
	change_ToWorkingDirectory(directory_name = pjoin(UserInput_main_dir,UserInput_output_dir_name))
	conn = sql.connect('seed_results.db')
	conn2 = sql.connect(pjoin(UserInput_main_dir,'genomes.db'))

	## prepare taxonomy data ##
	df_meta = pd.read_sql_query("SELECT * from gb_metadata", conn2)
	df_meta['genus'] = [t.split(' ')[0].lower() for t in df_meta['organism'].tolist()]
	df_meta['species'] = [t.split(' ')[1].lower() for t in df_meta['organism'].tolist()]
	df_meta = df_meta[['assembly_id','genus','species','organism']]

	## add ortholog, mmseqs cluster ids, dbscan labels, and filtered hit results to the extracted seed regions ##
	df_sr = pd.read_sql_query("SELECT * from seed_regions", conn) 
	df_cr = pd.read_sql_query("SELECT * from dbscan_label_representatives", conn)
	df_or = pd.read_sql_query("SELECT * from select_mcl_results", conn)

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


	# ---- Write comprehensive region & cluster data table ---- #
	df_sr.to_csv(pjoin('region_clusters','full_region.csv'),index=False)

	# ---- Make metadata outputs ---- #
	melt_LabelDissimilarty(df_cr=df_cr, df_jac=df_reorder)
	melt_IdentityCoverage(df_sr=df_sr)
	df_meta.to_csv(pjoin('results','rtable_genome_input_metadata.csv'),index=False)
	df_sr[['assembly_id','region_id','dbscan_label']].drop_duplicates().to_csv(pjoin('results','rtable_summarized_region_clusters.csv'), index=False)

	# close SQL connections
	conn.close()
	conn2.close()



## --------------End script---------------------##









