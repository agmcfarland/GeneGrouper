
import os
from os.path import join as pjoin
import re
from Bio.Seq import Seq
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
import pandas as pd
pd.options.mode.chained_assignment = None
import math
import numpy as np
import warnings
import time
import sqlite3 as sql
from collections import defaultdict
import gc
import sys
from gtr_utils import change_ToWorkingDirectory, make_OutputDirectory, merge_ManyFiles, multiprocesssing_Submission, generate_AssemblyId, chunks
from sklearn.cluster import DBSCAN
from sklearn.metrics.pairwise import pairwise_distances
from sklearn import metrics
import scipy.stats as stats
from scipy.cluster.hierarchy import dendrogram, linkage, cut_tree
import subprocess
from subprocess import DEVNULL, STDOUT, check_call

def cluster_DBSCAN(jac_dist, table_output_suffix, min_group_size):
	'''
	Run dbscan on a precomputed jaccard distance matrix. 
	Run through step-wise increments of epsilon values. Currently hardcoded to be 0-1 with 0.5 step increments.
	Check for case where only a single cluster is found.
	Return the 'best' clustering result. Currently as measured by max calinksi-harabasz score.
	'''
	# cluster and extract cluster label assignments

	if min_group_size == 'default':
		min_group_size = math.floor(np.log(len(jac_dist)))
	else:
		min_group_size = int(min_group_size)

	print('Using a minimum group size of {}'.format(min_group_size))

	cluster_label_assignments = {}
	cluster_stats=[]
	# eps_val_list=[round(n,2) for n in np.arange(0.05,1,0.05)]
	eps_val_list=[round(n,2) for n in np.arange(0.02,1,0.02)]
	for x, eps_val in enumerate(eps_val_list):
		#
		db = DBSCAN(eps=eps_val,min_samples=min_group_size,metric='precomputed').fit(jac_dist)
		core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
		core_samples_mask[db.core_sample_indices_] = True
		labels = db.labels_
		if len(set(labels))>1:
			# Number of clusters in labels, ignoring noise if present.
			n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
			n_noise_ = list(labels).count(-1)
			cluster_stats.append(
				{'eps_val':eps_val,
				'n_clusters':n_clusters_,
				'n_noise':n_noise_,
				'silh_score':round(metrics.silhouette_score(jac_dist,labels),2),
				'cali_score':round(metrics.calinski_harabasz_score(jac_dist, labels),2),
				'dbi_score': round(metrics.davies_bouldin_score(jac_dist,labels),2)})

		else:
			n_clusters_ = 1
			n_noise_ = list(labels).count(-1)
			cluster_stats.append(
				{'eps_val':eps_val,
				'n_clusters':n_clusters_,
				'n_noise':n_noise_,
				'silh_score':0,
				'cali_score':0,
				'dbi_score':0})
			# the actual per seed region cluster labels
		cluster_label_assignments['{}'.format(str(eps_val))] = pd.Series(labels)
		# Run DBSCAN again with the lowest possible input if the lowest input in the tested range has a cluster of 1
		#under construction
		print(eps_val)

	## Save two tables: the clustering statistic results and the labels for each of the clusters
	df_cluster_stats = pd.DataFrame(cluster_stats)
	df_cluster_stats.to_csv(pjoin('internal_data','dbscan_cluster_statistics{}.csv'.format(table_output_suffix)),index=False)

	df_cluster_labels = pd.DataFrame(cluster_label_assignments)
	df_cluster_labels.to_csv(pjoin('internal_data','dbscan_cluster_labels{}.csv'.format(table_output_suffix)),index=False)

	# check for case where only 1 cluster in across all tested eps values. if this is true, then return the lowest tested eps value as the 'best' value
	if df_cluster_stats['cali_score'].sum() == 0:
		best_cluster_result = df_cluster_stats.iloc[0,0]
	else:
		best_cluster_result = df_cluster_stats[df_cluster_stats['cali_score']==max(df_cluster_stats['cali_score'])]['eps_val'].iloc[0]
	df_cr = df_cluster_labels[[str(best_cluster_result)]]
	return(df_cr)



def determine_ClusterRepresentative(df_cr,df_c):
	'''
	Find the most representative seed_region for each cluster. 
	Currently using the seed_region with the lowest mean dissimilarity 
	Uses the index from df_c, which is a seed_regionX<ortho clustering id> presence/absence matrix. This makes subsequent matrix operations fast.
	IMPORTANT: Uses the jac_dist global 
	'''

	df_clusterep = pd.DataFrame()

	df_cr['seed_region_id'] = df_c.index.tolist()
	for cluster_id in sorted((df_cr.iloc[:,0]).unique()):
		# get the index postions of each seed_region_id for the current cluster
		# make a relational set of indices for the subsetted cluster
		seed_region_indices = np.array(df_cr[df_cr.iloc[:,0]==cluster_id].index.tolist())
		seed_region_indices = np.column_stack((seed_region_indices,np.asarray(range(0,len(seed_region_indices)))))

		# filter the full distance matrix for those in cluster
		interdist_array =  np.take(jac_dist,seed_region_indices[:,0],axis=0)
		interdist_array =  np.take(interdist_array,seed_region_indices[:,0],axis=1)

		n_seedregions = len(interdist_array)
		# for each row, calculate the mean dissimilarity. subtract by 1 to correct for self comparison.
		mean_dissimilarity_array = np.asarray([(sum(i)/(n_seedregions-1)) for i in interdist_array])
		mean_dissimilarity_array = np.column_stack((np.asarray(seed_region_indices),mean_dissimilarity_array))
		# sort the list to retrieve entry with lowest mean
		lowest_mean_dis_info = mean_dissimilarity_array[mean_dissimilarity_array[:,2].argsort()][0]
		# retrieve the row with dissimilarites that correspond to the region with the lowest mean dissimilarity. use the reduced index to find.
		# columns for lowest_relative_dis_info goes: full index, subsetted index, mean dissimilarity
		lowest_relative_dis = interdist_array[:,int(lowest_mean_dis_info[1])]
		## put all the pieces together into a single array:
		# add column with the representative seed_region index
		mean_dissimilarity_array = np.column_stack((np.full(shape=n_seedregions,fill_value=lowest_mean_dis_info[0]),mean_dissimilarity_array))
		# add column with cluster_id 
		mean_dissimilarity_array = np.column_stack((np.full(shape=n_seedregions,fill_value=cluster_id),mean_dissimilarity_array))
		# add column with disimilartiy values relative to the region with lowest dissimilarity
		mean_dissimilarity_array = np.column_stack((mean_dissimilarity_array,lowest_relative_dis))
		# add to master pandas dataframe. label representative uses number from full_index_pos
		colnames = ['dbscan_label','label_representative','full_index_pos','relative_index_pos','mean_dissimilarity','representative_relative_dissimilarity']
		df_clusterep = df_clusterep.append(pd.DataFrame(mean_dissimilarity_array,columns=colnames))


	# add seed_region_id to the df_clusterep dataframe. save output
	df_clusterep = df_clusterep.sort_values(by='full_index_pos')
	df_clusterep['seed_region_id'] = df_cr['seed_region_id'].tolist()
	return(df_clusterep)


def update_FinalCluster(df, cluster_dict):
	'''
	'''
	# for indices that have a clustering result, update the dictionary key.
	for c in df[df.iloc[:,0]>-1].iterrows():
		cluster_dict[c[0]] = c[1][0] # update index with cluster value
	return(cluster_dict)



def recluster_UnclusteredRegions(df_lastcluster, max_cluster_val, recluster_round):
	'''
	IMPORTANT: Uses the jac_dist global 
	'''
	# get the jaccard distance matrix index positions for all unclustered regions
	unclustered_positions = df_lastcluster[df_lastcluster.iloc[:,0] == -1]
	unclustered_positions = unclustered_positions.rename_axis('position').reset_index()
	# convert to an array that is used to filter the jaccard distance matrix
	unclustered_index = np.array(unclustered_positions.position.tolist())
	ju = jac_dist[unclustered_index[:, None], unclustered_index]
	df_unclustered = cluster_DBSCAN(jac_dist = ju, table_output_suffix = '_{}'.format(recluster_round)) 
	# add the jaccard distance matrix positions to the clustering assignments
	df_unclustered.index = unclustered_positions.position
	# update clustering assingments to be 1+ the max_cluster_val when clustering assingment is not -1
	df_unclustered.iloc[:,0] = np.where(df_unclustered.iloc[:,0] >= 0, df_unclustered.iloc[:,0] + max_cluster_val + 1, -1)
	#
	return(df_unclustered)



## --------------Start script---------------------##

def ClusterRegions(
	UserInput_main_dir,
	UserInput_output_dir_name,
	UserInput_reclustering_iterations,
	UserInput_min_group_size
	):
	'''
	'''

	# # #### TROUBLESHOOTING INPUTS #####		
	# UserInput_main_dir = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/syntenease_project/gtr/testbed/dataset1/test1'
	# UserInput_output_dir_name = 'pdua'
	# UserInput_reclustering_iterations = 1
	# pd.set_option('display.max_columns', None)
	# # #################################		

	# make jac_dist global so it isn't copied in functions
	global jac_dist

	## read in ortholog relationship information and assign each locus_tag it's seed_region identifier##
	change_ToWorkingDirectory(directory_name = pjoin(UserInput_main_dir,UserInput_output_dir_name))
	conn = sql.connect('seed_results.db')
	df_orthos = pd.read_sql_query("SELECT * from select_mcl_results", conn)
	df_seed_regions = pd.read_sql_query("SELECT * from seed_regions", conn)

	## pseudogenes are present in the seed_regions df but not in the mcl output ##
	## assign all pseudo genes a single cluster value that is 1+ the maximum cluster_id
	df_m = df_orthos[['orf_id','ortho_cluster_id']].merge(df_seed_regions[['locus_tag','region_id']],  left_on = 'orf_id',  right_on =  'locus_tag', how = 'outer')
	df_m['orf_id'] = df_m['orf_id'].fillna(df_m['locus_tag']) #this identifies pseudo genes
	max_ortho_cluster_id_val = df_m['ortho_cluster_id'].max(skipna = True)
	df_m['ortho_cluster_id'] = df_m['ortho_cluster_id'].fillna(max_ortho_cluster_id_val+1) #this adds pseudogene identifier
	df_m['ortho_cluster_id'] = df_m['ortho_cluster_id'].astype(int)

	if len(df_m[df_m['ortho_cluster_id'].isnull()]) != 0:
		print('Some sequences dont have ortholog ids')

	## construct a jaccard distance matrix ##
	## create a presence absence table using region_id and ortho_cluster_id. cast the table so that it is binary ##
	## and then convert it to a jaccard distance ##
	df_m = df_m[['region_id','ortho_cluster_id']]
	df_m['count'] = 1
	df_m = df_m.pivot_table('count',index='region_id',columns='ortho_cluster_id').fillna(0).astype(int)
	#
	ortho_cluster_pa_array = df_m.to_numpy()
	jac_dist = pairwise_distances(ortho_cluster_pa_array, metric = "jaccard")
	## save distance matrix to contain labels and columns
	df_jd = pd.DataFrame(jac_dist,columns=df_m.index.tolist(),index=df_m.index.tolist())
	df_jd.to_csv(pjoin('internal_data','seed_regions_jaccard_dist.csv'))

	print ('Clustering {} regions'.format(len(jac_dist)))
	## Run DBSCAN on jaccard distance using different epsilon values ##
	## Returns info on clustering quality, statistics, and for the "best" cluster, dissimilarities ##
	df_cr = cluster_DBSCAN(jac_dist = jac_dist, table_output_suffix = '_0', min_group_size = UserInput_min_group_size) 


	# intialize a dictionary with all -1 clustering assignments for all cluster regions indices.
	final_cluster = {i:-1 for i in df_cr.index.tolist()}
	# add clustering results from first round
	final_cluster = update_FinalCluster(df = df_cr, cluster_dict = final_cluster)

	# reclustering procedure START
	if UserInput_reclustering_iterations > 0:
		for r in range(UserInput_reclustering_iterations):
			
			# check that -1 exist
			if sum(value == -1 for value in final_cluster.values()) <= 1:
				print('No more unclustered regions to cluster')
				break

			# Recluster
			print('Clusering round {}'.format(r + 1))
			print('Clustering {} regions'.format(sum(value == -1 for value in final_cluster.values())))
			df_cr =  recluster_UnclusteredRegions(df_lastcluster = df_cr, max_cluster_val = max(final_cluster.values()), recluster_round = r + 1)
			final_cluster = update_FinalCluster(df = df_cr, cluster_dict = final_cluster)

			# exit loop once recluster limit reached
			if r + 1 == UserInput_reclustering_iterations:
				break
	# reclustering procedure END

	# convert final_cluster dictionary to pandas df
	df_final_cluster = pd.DataFrame(final_cluster.items(), columns = ['index','clustering'])
	df_final_cluster = df_final_cluster.set_index('index')

	# find cluster representatives
	df_clusterep = determine_ClusterRepresentative(df_cr=df_final_cluster,df_c=df_m)

	## Store representative cluster labels in a new table in the seed_results sql database ##
	## also write to internal_data ##
	conn = sql.connect('seed_results.db')
	df_clusterep.to_sql(name='dbscan_label_representatives',con=conn, if_exists='replace',index_label='seed_region_id',index=False)
	df_clusterep.to_csv(pjoin('internal_data','dbscan_label_representatives.csv'),index=False)

	## remove large dataframes from memory ##
	del [[df_jd,df_clusterep,ortho_cluster_pa_array,jac_dist]]
	gc.collect()

	## close sql connection
	conn.close()


## --------------End script---------------------##

















