

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
from gtr_utils import change_ToWorkingDirectory, make_OutputDirectory, merge_ManyFiles, multiprocesssing_Submission, generate_AssemblyId, chunks
from output_RegionResults import subset_SimilarityMatrix
import subprocess
from subprocess import DEVNULL, STDOUT, check_call

def convert_ToGGGenesInput(df_m,filter_on):
	'''
	'''
	df_gggenes = pd.DataFrame()
	for clabel in df_m[filter_on].unique().tolist():
		print(str(int(clabel)))
		df_cl = df_m[df_m[filter_on] == clabel]
		df_cl = df_cl.sort_values(by=['new_pos'],ascending=True)

		# align genes to normalized position of the seed
		start_bp = df_cl[df_cl['is_seed']==1]['cds_start1'].item()

		# process downstream list
		cds_list = list()
		for row in df_cl[df_cl['new_pos']>=0].iterrows():
			row_pos1 = len(row[1])-2
			row_pos2 = len(row[1])-1
			cds_list.append(row[1][row_pos1])
			cds_list.append(row[1][row_pos2])

		downstream_xlist = [0]
		while len(cds_list) >1:
			if cds_list[0] == start_bp:
				x1 = cds_list[1] - cds_list[0]
				downstream_xlist.append(x1)
				del cds_list[0]
				continue

			x1 = cds_list[1] - cds_list[0] + x1
			del cds_list[0]
			downstream_xlist.append(x1)

		# process upstream list
		cds_list = list()
		for row in df_cl[df_cl['new_pos']<=0].iterrows():
			row_pos1 = len(row[1])-2
			row_pos2 = len(row[1])-1
			cds_list.append(row[1][row_pos1])
			cds_list.append(row[1][row_pos2])

		cds_list.reverse()
		del cds_list[0]
			
		upstream_xlist = [0]
		while len(cds_list) >1:
			if cds_list[0] == start_bp:
				x1 = cds_list[1] - cds_list[0]
				upstream_xlist.append(x1)
				del cds_list[0]
				continue

			x1 = cds_list[1] - cds_list[0] + x1
			del cds_list[0]
			upstream_xlist.append(x1)
		del upstream_xlist[0]
		upstream_xlist.reverse()

		# map values back to dataframe
		upstream_xlist = np.asarray(upstream_xlist).reshape(int(len(upstream_xlist)/2),2)
		downstream_xlist = np.asarray(downstream_xlist).reshape(int(len(downstream_xlist)/2),2)
		full_xlist = np.concatenate((upstream_xlist,downstream_xlist), axis=0)

		df_cl['norm_start'] = full_xlist[:,0]
		df_cl['norm_end'] = full_xlist[:,1]
		df_gggenes = df_gggenes.append(df_cl[['locus_tag','norm_start','norm_end']])

	return(df_gggenes)

def align_NormPositions(df_m):
	'''
	Align all positions to start at zero by taking the absolute max value for the start or end of the seed gene
	'''
	max_vals = df_m[df_m['is_seed']==1][['norm_start','norm_end']].abs().max(axis=1)
	df_max_vals = df_m[df_m['is_seed']==1][['region_id']]
	df_max_vals['max_val'] = max_vals
	df_m = df_m.merge(df_max_vals, on='region_id')
	df_m['norm_start'] = np.where(df_m['global_strand']== -1, df_m['norm_start']+df_m['max_val'],df_m['norm_start'])
	df_m['norm_end'] = np.where(df_m['global_strand']== -1, df_m['norm_end']+df_m['max_val'],df_m['norm_end'])
	return(df_m)



def plot_RepresentativeRegions(df_sr, df_cr):
	'''
	Produce the tables necessary to plot the representative regions along with associated metadata
	'''

	## format for plotting
	df_m = df_sr[df_sr['label_representative'] == df_sr['full_index_pos']][['region_id','assembly_id','contig_id','locus_tag','cds_start','cds_end','strand','global_strand','pseudo_check','pident','qcovs','evalue',
	'dbscan_label','mean_dissimilarity','representative_relative_dissimilarity','ortho_cluster_id',
	'genus','species','refseq_product','refseq_gene','new_pos','is_seed','order']]

	#reverse gene order for global_strand == -1 regions by conditionalally multplying by -1
	df_m['cds_start1'] = np.where(df_m['global_strand']==-1, df_m['cds_start']*-1, df_m['cds_start'])
	df_m['cds_end1'] = np.where(df_m['global_strand']==-1, df_m['cds_end']*-1, df_m['cds_end'])

	#convert to gggenes input format and save
	df_gggenes = convert_ToGGGenesInput(df_m=df_m,filter_on='dbscan_label')
	df_m = df_m.merge(df_gggenes,on='locus_tag')
	df_m = align_NormPositions(df_m=df_m)

	df_m.to_csv(pjoin('results','rtable_region_representatives.csv'),index=None)




def get_UniqueClustersWithinCluster(df_cr, df_sr, cluster_label_id):
	'''
	'''
	## subset for the cluster of interest
	rep_list = df_cr[df_cr['dbscan_label']==cluster_label_id]['seed_region_id'].tolist()
	df_jac = subset_SimilarityMatrix(rep_list=rep_list)
	df_jac['seed_region_id'] = rep_list
	df_mjac = pd.melt(df_jac, id_vars=['seed_region_id'])

	## group all subclusters with 0 dissimilarity
	store_unique_cwc = defaultdict(list)
	found_regions = list()
	for c in df_mjac['seed_region_id'].tolist():
		if c in found_regions:
			continue
		temp_list = df_mjac[ (df_mjac['seed_region_id'] == c) & (df_mjac['value'] == 0)]['variable'].tolist()
		store_unique_cwc[c] = temp_list
		[found_regions.append(c1) for c1 in temp_list if c1 not in found_regions]
		
	## get counts of each unique region group
	list(store_unique_cwc.keys())
	store_unique_cwc_count = list()
	for l in store_unique_cwc.keys():
		store_unique_cwc_count.append(len(store_unique_cwc[l]))
	df_cwc = pd.DataFrame({
		'seed_region_id':list(store_unique_cwc.keys()),
		'count':store_unique_cwc_count})

	# add in relative dissimilarity data, sort by ascending dissimilarity, add unique identifier, and store
	df_cwc = df_cwc.merge(df_cr[['seed_region_id','dbscan_label','representative_relative_dissimilarity']], on='seed_region_id',how='inner')
	df_cwc = df_cwc.sort_values(by=['representative_relative_dissimilarity','count'], ascending=[True,False])
	df_cwc['cwc_id'] = [i for i in range(len(df_cwc))]
	df_cwc = df_cwc.rename(columns={'seed_region_id':'region_id'},inplace=False)
	df_cwc.to_csv(pjoin('results','rtable_cwc_dis_counts.csv'),index=False)

	# merge data with df_sr dataframe and store
	df_cwc = df_cwc.merge(df_sr[['region_id','assembly_id','contig_id','locus_tag','cds_start','cds_end','global_strand','strand','pseudo_check','pident','qcovs','evalue','ortho_cluster_id',
	'genus','species','refseq_product','refseq_gene','new_pos','is_seed','order']],
	left_on='region_id',right_on='region_id',how='inner')
	# df_cwc.to_csv(pjoin('results','rtable_cwc_extract_regions.csv'),index=False)

	## convert to ggenes format and save
	#reverse gene order for global_strand == -1 regions
	df_cwc['cds_start1'] = np.where(df_cwc['global_strand']==-1, df_cwc['cds_start']*-1, df_cwc['cds_start'])
	df_cwc['cds_end1'] = np.where(df_cwc['global_strand']==-1, df_cwc['cds_end']*-1, df_cwc['cds_end'])
	#
	df_gggenes = convert_ToGGGenesInput(df_m=df_cwc,filter_on='cwc_id')
	df_cwc = df_cwc.merge(df_gggenes,on='locus_tag')
	df_cwc = align_NormPositions(df_m=df_cwc)
	df_cwc.to_csv(pjoin('results','rtable_cwc_regions.csv'),index=None)




#------- start script ------#
def MakeVisualizations(
	UserInput_main_dir,
	UserInput_output_dir_name,
	UserInput_cluster_label_id,
	UserInput_visualization_type,
	UserInput_script_path):
	'''
	'''
	change_ToWorkingDirectory(directory_name = pjoin(UserInput_main_dir,UserInput_output_dir_name))
	conn = sql.connect('seed_results.db')
	df_sr = pd.read_csv(pjoin('region_clusters','full_region.csv'))
	df_cr = pd.read_sql_query("SELECT * from dbscan_label_representatives", conn)

	# if os.path.isfile(pjoin(UserInput_script_path,'visualize_main_region.R')) == True:

	if UserInput_visualization_type == 'main':
		print('Arranging clusters')
		plot_RepresentativeRegions(df_sr=df_sr,df_cr=df_cr)
		os.system('Rscript {}/visualize_main_region.R {} {}'.format(UserInput_script_path,pjoin(os.getcwd(),'results'), pjoin(os.getcwd(),'visualizations') ) )
		# check_call(['Rscript', '{}/visualize_main_region.R'.format(UserInput_script_path), '{}'.format(pjoin(os.getcwd(),'results')), '{}'.format(pjoin(os.getcwd(),'visualizations')) ], stdout=STDOUT, stderr=STDOUT) #DEVNULL
		# check_call(['Rscript', 'visualize_main_region.R', '{}'.format(pjoin(os.getcwd(),'results')), '{}'.format(pjoin(os.getcwd(),'visualizations')) ], stdout=STDOUT, stderr=STDOUT) #DEVNULL
		return


	if UserInput_visualization_type == 'region_cluster':
		get_UniqueClustersWithinCluster(df_cr=df_cr, df_sr=df_sr, cluster_label_id=UserInput_cluster_label_id)
		os.system('Rscript {}/visualize_cluster.R {} {}'.format(UserInput_script_path,pjoin(os.getcwd(),'results'), pjoin(os.getcwd(),'visualizations') ) )
		return

#------- endscript ------#





