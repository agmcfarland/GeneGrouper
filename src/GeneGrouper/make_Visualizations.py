

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

	df_m.to_csv(pjoin('internal_data','rtable_region_representatives.csv'),index=None)




def get_UniqueClustersWithinCluster(df_cr, df_sr, cluster_label_id):
	'''
	'''
#troubleshooting
# os.chdir('/Users/owlex/Downloads/gg_tutorial/example_search/mexb')


# conn = sql.connect('seed_results.db')
# df_sr = pd.read_csv(pjoin('internal_data','full_region.csv'))
# df_cr = pd.read_sql_query("SELECT * from dbscan_label_representatives", conn)
# cluster_label_id = 10

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

	## extract the index of subgroup representatives and gene clusters within them to a table
	df_subind = pd.DataFrame(list(store_unique_cwc.values()), index=store_unique_cwc.keys())
	df_subind['representative_subgroup'] = df_subind.index
	df_subind = pd.melt(df_subind, id_vars=['representative_subgroup'])
	df_subind = df_subind.drop(labels=['variable'],axis=1)
	df_subind = df_subind.rename(columns={'value':'subgroup_member'})
	df_subind = df_subind.dropna()
		
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
	df_cwc.to_csv(pjoin('internal_data','rtable_cwc_dis_counts.csv'),index=False)

	# merge data with df_sr dataframe and store
	df_cwc = df_cwc.merge(df_sr[['region_id','assembly_id','contig_id','locus_tag','cds_start','cds_end','global_strand','strand','pseudo_check','pident','qcovs','evalue','ortho_cluster_id',
	'genus','species','refseq_product','refseq_gene','new_pos','is_seed','order']],
	left_on='region_id',right_on='region_id',how='inner')

	## convert to ggenes format and save
	#reverse gene order for global_strand == -1 regions
	df_cwc['cds_start1'] = np.where(df_cwc['global_strand']==-1, df_cwc['cds_start']*-1, df_cwc['cds_start'])
	df_cwc['cds_end1'] = np.where(df_cwc['global_strand']==-1, df_cwc['cds_end']*-1, df_cwc['cds_end'])
	#
	df_gggenes = convert_ToGGGenesInput(df_m=df_cwc,filter_on='cwc_id')
	df_cwc = df_cwc.merge(df_gggenes,on='locus_tag')
	df_cwc = align_NormPositions(df_m=df_cwc)
	df_cwc.to_csv(pjoin('internal_data','rtable_cwc_regions.csv'),index=None)

	## save a copy of subgroup data in the subgroup folder with the cluster_label_id as suffix. This is for the user to use in downstream applications
	df_subgroup_ref = df_cwc[['region_id','cwc_id']].drop_duplicates()
	df_subgroup_ref = df_subgroup_ref.merge(df_subind, left_on='region_id',right_on='representative_subgroup')
	df_subgroup_ref = df_subgroup_ref.drop(labels='region_id',axis=1)
	# rename so it's easier for the user to know what columns mean
	df_subgroup_ref = df_subgroup_ref.rename(columns={'cwc_id':'subgroup_id'})
	df_cwc = df_cwc.rename(columns={'cwc_id':'subgroup_id'})

	# if the subgroup folder doesn't exist, save it
	make_OutputDirectory(new_directory='subgroups')
	df_subgroup_ref.to_csv(pjoin('subgroups','subgroup_key_g{}.csv'.format(str(cluster_label_id))),index=None)
	df_cwc.to_csv(pjoin('subgroups','subgroup_regions_g{}.csv'.format(str(cluster_label_id))),index=None)



def retrieve_RepresentativeSeeds(tip_label_type, df_cr):
	'''
	Get the amino acid sequence of each representative seed EXCEPT for the group -1 representative
	'''
	# pd.set_option('display.max_columns', None)
	# tip_label_type = 'full'
	# os.chdir('/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/syntenease_project/gtr/testbed/dataset4/test2/mexb')

	# get just the representative sequences into a new df 
	df_rep = df_cr[(df_cr['label_representative'] == df_cr['full_index_pos']) & (df_cr['dbscan_label'] != -1 )][['dbscan_label','seed_region_id']]

	# create tip_label column with specified content
	if tip_label_type == 'full':
		df_rep['tip_label'] = df_rep['seed_region_id'] + '--g' + df_cr['dbscan_label'].astype(int).astype(str)
	if tip_label_type == 'group':
		df_rep['tip_label'] = 'g'+df_rep['dbscan_label'].astype(int).astype(str)

	# write representative sequences to faa file
	with open(pjoin('internal_data','representative_seed_sequences.faa'), 'w') as outfile:
		for record in SeqIO.parse('group_region_seqs.faa', 'fasta'):
			if record.id in df_rep['seed_region_id'].tolist():
				tip_entry = df_rep[df_rep['seed_region_id'] == record.id]['tip_label'].item()
				outfile.write('>{}\n'.format(tip_entry))
				outfile.write('{}\n'.format(str(record.seq)))


def create_RepresentativeSeedTree():
	'''
	Run mafft and fasttree on the representative seed sequences
	'''
	os.system('mafft {} > {}'.format(pjoin('internal_data','representative_seed_sequences.faa'), pjoin('internal_data','representative_seed_sequences.aln')))
	os.system('FastTree {} > {}'.format(pjoin('internal_data','representative_seed_sequences.aln'), pjoin('internal_data','representative_seed_sequences.nwk')))


#------- start script ------#
def MakeVisualizations(
	UserInput_main_dir,
	UserInput_output_dir_name,
	UserInput_cluster_label_id,
	UserInput_visualization_type,
	UserInput_script_path,
	UserInput_image_format,
	UserInput_tip_label_type,
	UserInput_tip_label_size):
	'''
	'''
	change_ToWorkingDirectory(directory_name = pjoin(UserInput_main_dir,UserInput_output_dir_name))
	conn = sql.connect('seed_results.db')
	df_sr = pd.read_csv(pjoin('internal_data','full_region.csv'))
	df_cr = pd.read_sql_query("SELECT * from dbscan_label_representatives", conn)

	# main visualization
	if UserInput_visualization_type == 'main':
		print('Arranging clusters')
		plot_RepresentativeRegions(df_sr=df_sr,df_cr=df_cr)
		print('plotting in R...\n')
		os.system('Rscript {}/visualize_main_region.R {} {} {}'.format(UserInput_script_path, pjoin(os.getcwd(),'internal_data'), pjoin(os.getcwd(),'visualizations'), UserInput_image_format ))
		# check_call(['Rscript', '{}/visualize_main_region.R'.format(UserInput_script_path), '{}'.format(pjoin(os.getcwd(),'internal_data')), '{}'.format(pjoin(os.getcwd(),'visualizations')) ], stdout=STDOUT, stderr=STDOUT) #DEVNULL
		# check_call(['Rscript', 'visualize_main_region.R', '{}'.format(pjoin(os.getcwd(),'internal_data')), '{}'.format(pjoin(os.getcwd(),'visualizations')) ], stdout=STDOUT, stderr=STDOUT) #DEVNULL
		print('\nAll outputs files written to {}\n'.format(os.getcwd()))
		return

	## group visualization
	if UserInput_visualization_type == 'group':
		get_UniqueClustersWithinCluster(df_cr=df_cr, df_sr=df_sr, cluster_label_id=UserInput_cluster_label_id)
		print('plotting in R...\n')
		os.system('Rscript {}/visualize_group.R {} {} {}'.format(UserInput_script_path, pjoin(os.getcwd(),'internal_data'), pjoin(os.getcwd(),'visualizations'), UserInput_image_format ))
		return

	## tree visualization
	if UserInput_visualization_type == 'tree':
		print('Making a phylogenetic tree from representative seed sequences')
		retrieve_RepresentativeSeeds(tip_label_type=UserInput_tip_label_type, df_cr=df_cr)
		create_RepresentativeSeedTree()
		print('plotting in R...\n')
		os.system('Rscript {}/visualize_tree.R {} {} {} {}'.format(UserInput_script_path, pjoin(os.getcwd(),'internal_data'), pjoin(os.getcwd(),'visualizations'), UserInput_image_format,  UserInput_tip_label_size))
		return

#------- endscript ------#




