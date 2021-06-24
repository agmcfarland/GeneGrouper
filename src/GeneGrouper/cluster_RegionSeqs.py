
import os
from os.path import join as pjoin
import re
import multiprocessing as mp
from multiprocessing import Pool
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
import subprocess
from subprocess import DEVNULL, STDOUT, check_call

def count_OrfSeqs(input_fasta_file):
	'''
	returns the total number of sequences in a fasta. starts at 1. Used to create an index for the large merged_aa fasta file.
	'''
	number_of_seqs = 0
	for record in SeqIO.parse(input_fasta_file,'fasta'):
		number_of_seqs += 1
	return(number_of_seqs)


def ava_Blast(query,db_name,output_filename):
	'''
	'''
	# os.system('blastp -query {} -db {} -max_target_seqs 100 -evalue 1e-6 -outfmt "10 qseqid sseqid mismatch positive gaps ppos pident qcovs evalue bitscore qframe sframe sstart send slen qstart qend qlen" -num_threads 1  -out {}'.format(query,db_name,output_filename))
	check_call(['blastp', '-query', query, '-db', db_name, '-max_target_seqs', '100', '-evalue', '1e-6', '-outfmt', "10 qseqid sseqid mismatch positive gaps ppos pident qcovs evalue bitscore qframe sframe sstart send slen qstart qend qlen", '-num_threads', '1', '-out', output_filename], stdout=DEVNULL, stderr=STDOUT)


def chunk_BlastInput(chunksize,fasta_to_chunk):
	'''
	chunk up the merged_aa file to form several query fasta files
	The input fasta file can either be a unclustered or preclustered orf fasta
	'''
	fasta_to_chunk_input = '{}_DB_clu_rep.fasta'.format(fasta_to_chunk)
	total_aa_seqs = count_OrfSeqs(input_fasta_file = fasta_to_chunk_input)


	## populate a list with the number of chunks and the remaining chunk value.
	## calculate the cumulative value for each row, this forms an index that is equivalent to the number from enumerate
	## note: -1 is necessary on the final input to match enumerate counter
	aa_chunk_splitsize = int(total_aa_seqs/chunksize)

	aa_chunk_size_index = [aa_chunk_splitsize for x in range(chunksize)]
	aa_chunk_size_index.append(total_aa_seqs-sum(aa_chunk_size_index)-1)
	if sum(aa_chunk_size_index) != total_aa_seqs-1:
		print('BAD MESSAGE: chunking error')
	aa_chunk_size_index = pd.Series(aa_chunk_size_index).cumsum().tolist()

	## chunk_tracker is updated and temp_sequences is cleared everytime the chunk size is reached and a file with the data is written
	chunk_filesuffix = 'aa_chunk_'
	chunk_tracker = 0
	temp_sequences = []
	for e, record in enumerate(SeqIO.parse(fasta_to_chunk_input,'fasta')):
		if e == aa_chunk_size_index[0]:
			temp_sequences.append('>'+record.id+'\n')
			temp_sequences.append(str(record.seq)+'\n')
			with open('{}{}.txt'.format(chunk_filesuffix, chunk_tracker),'w') as outfile:
				[outfile.write(i) for i in temp_sequences]
			temp_sequences = []
			chunk_tracker += 1
			aa_chunk_size_index = aa_chunk_size_index[1:]
		else:
			temp_sequences.append('>'+record.id+'\n')
			temp_sequences.append(str(record.seq)+'\n')


def blast_ChunkedSeqs(chunksize, fasta_to_chunk_prefix,working_directory,UserInput_processes):
	'''
	Returns tabular blast results for a chunked input blast search
	'''
	chunk_BlastInput(chunksize=chunksize,fasta_to_chunk=fasta_to_chunk_prefix)
	query_chunks = [[f,fasta_to_chunk_prefix,f.replace('.txt','.csv')] for f in os.listdir(working_directory) if f.startswith('aa_chunk_') if f.endswith('.txt')]
	multiprocesssing_Submission(function_name=ava_Blast,submission_list=query_chunks,processor_n=UserInput_processes,output_type='file')
	merged_ava_output = '{}_avablast'.format(fasta_to_chunk_prefix)
	merge_ManyFiles(input_filepath=working_directory,output_filepath=working_directory,wildcard_search='aa_chunk_*.csv',output_filename=merged_ava_output+'.csv')


def cluster_MCLOrtho(df_mcl, similarity_scoring_method, inflation_values_list):
	'''
	Convert a tabular blast result into a .abc formatted dataframe that is then clustered using mcl at different inflation values. save a csv of cluster relationship at each inflation value
	'''
	# mcl_file_prefix='merged_aa'
	df_mcl = df_mcl[['qseqid','sseqid',similarity_scoring_method]]
	df_mcl.to_csv('merged_aa.abc',sep='\t',index=None,header=None)
	#
	# os.system("mcxload -abc {}.abc --stream-mirror --stream-neg-log10  -o seq_{}.mci -write-tab seq_{}.tab -stream-tf 'ceil(200)'".format(mcl_file_prefix,mcl_file_prefix,mcl_file_prefix))
	check_call( ["mcxload", "-abc", "merged_aa.abc", "--stream-mirror", "--stream-neg-log10",  "-o", "seq_merged_aa.mci", "-write-tab", "seq_merged_aa.tab", "-stream-tf", 'ceil(200)'], stdout=DEVNULL, stderr=STDOUT)
	#
	for inval in inflation_values_list:
		## find clusters for each inflation value
		i_val = inval.replace('.','')
		# os.system('mcl seq_{}.mci -I {}'.format(mcl_file_prefix,inval))
		check_call( ['mcl', 'seq_merged_aa.mci', '-I', str(inval)], stdout=DEVNULL, stderr=STDOUT)

		# os.system('mcxdump -icl out.seq_{}.mci.I{} -tabr seq_{}.tab -o dump.seq_{}.mci.I{}'.format(mcl_file_prefix,i_val,mcl_file_prefix,mcl_file_prefix,i_val))
		check_call( ['mcxdump', '-icl', 'out.seq_merged_aa.mci.I{}'.format(i_val), '-tabr', 'seq_merged_aa.tab', '-o', 'dump.seq_merged_aa.mci.I{}'.format(i_val) ], stdout=DEVNULL, stderr=STDOUT)
		## extract clusters into csv format
		## read txt file which has members of a cluster per line
		keeps = []
		with open('dump.seq_merged_aa.mci.I{}'.format(i_val),'r') as infile:
			for line in infile:
				keeps.append(line.replace('\n','').split('\t'))
		## place in dataframe
		df_clong = pd.DataFrame()
		for e,i in enumerate(keeps):
			cl_id = [e for x in range(len(i))]
			z = pd.DataFrame({'ortho_cluster_id':cl_id,'ortho_cluster_members':i})
			df_clong = pd.concat([df_clong,z],axis = 0)
		df_clong.to_csv('mcl_i{}.csv'.format(i_val),index=None)	



def add_UnclusteredOrfs(df_unclustered,df_infval_temp,inf_val):
	'''
	Not all ORFs may end up being found during AVA blast and therefore not part of the mcl results. Adding them back in here. 
	Each orf added in is naturally its own cluster representative
	'''
	new_cluster_range_lower = df_infval_temp['ortho_cluster_id'].max()+1
	new_cluster_range_upper = len(df_unclustered) + new_cluster_range_lower
	new_clusters = [i for i in range(new_cluster_range_lower,new_cluster_range_upper)]
	df_unclustered['ortho_cluster_id']=new_clusters
	df_unclustered['i_val'] = inf_val
	df_unclustered = df_unclustered[['ortho_cluster_id','ortho_cluster_members','i_val']]
	return(df_unclustered)



def assign_ClusterMembership(df_unclustered, inflation_values_list,):
	'''
	Each orf should have a cluster assignment, whether through MCL or manually by adding in unclustered orfs using add_UnclusteredOrfs
	Output is a combined mcl cluster assignments at the different inflation values in a long df format
	'''
	df_infval = pd.DataFrame()
	for inf_val in inflation_values_list:
		i_val = inf_val.replace('.','')
		df_infval_temp = pd.read_csv('mcl_i{}.csv'.format(i_val))
		df_infval_temp['i_val'] = inf_val
		if inf_val == inflation_values_list[0]:
			## the first entry to populate df_infval needs special attention
			df_unclustered = add_UnclusteredOrfs(df_unclustered=df_unclustered,df_infval_temp=df_infval_temp,inf_val=inf_val)
			df_infval_temp = pd.concat([df_infval_temp,df_unclustered])
			df_infval = df_infval_temp
			continue
		df_unclustered = add_UnclusteredOrfs(df_unclustered=df_unclustered,df_infval_temp=df_infval_temp,inf_val=inf_val)
		## do two merges: first to get the clustered and unclustered orfs into the same df, and then to add the temp orf cluster df to the master orf cluster df
		df_infval_temp = pd.concat([df_infval_temp,df_unclustered])
		df_infval = pd.concat([df_infval,df_infval_temp])
	df_infval['ortho_cluster_id'] = df_infval['ortho_cluster_id'].astype(int)
	return(df_infval)


def filter_MCLOorthos(df, UserInput_ortho_select_method):
	'''
	User selects the method to choose which MCL inflation value is chosen
	Returns a float 
		max_div > Returns a float that is t lowest inflation value with the highest ortholog group diversity.
		lowest_i > retruns a float of the lowest i_val used
	'''
	if UserInput_ortho_select_method == 'max_div':
		ortho_diversity_array = pd.DataFrame(df[['i_val','ortho_cluster_id']].groupby('i_val')['ortho_cluster_id'].nunique()).reset_index()
		max_ortho_div_val = ortho_diversity_array['ortho_cluster_id'].max()
		ortho_diversity_array['i_val'] = ortho_diversity_array['i_val'].astype(float)
		selected_ortho_mcl_i_val = ortho_diversity_array[ortho_diversity_array['ortho_cluster_id']==max_ortho_div_val].sort_values('i_val',ascending=True).nsmallest(1,'i_val')['i_val'].item()
	if UserInput_ortho_select_method == 'lowest_i':
		selected_ortho_mcl_i_val = 1.2
	return(selected_ortho_mcl_i_val)



## --------------Start script---------------------##

def ClusterRegionSeqs(
	UserInput_main_dir,
	UserInput_output_dir_name,
	UserInput_faa_chunk_size,
	UserInput_blast_filter_params,
	UserInput_inflation_values,
	UserInput_similarity_scoring_method,
	UserInput_linclust_settings,
	UserInput_ortho_select_method,
	UserInput_processes
	):


	## make an mmseqs database ##
	change_ToWorkingDirectory(directory_name=pjoin(UserInput_main_dir,UserInput_output_dir_name,'ortholog_clusters'))
	make_OutputDirectory(new_directory='tmp_region_orfs')
	os.system('mmseqs createdb merged_aa.faa merged_aa_DB -v 2')
	os.system('mmseqs linclust merged_aa_DB merged_aa_DB_clu tmp_region_orfs {} -v 2'.format(UserInput_linclust_settings))
	os.system('mmseqs createsubdb merged_aa_DB_clu merged_aa_DB merged_aa_DB_clu_rep -v 2')
	os.system('mmseqs convert2fasta merged_aa_DB_clu_rep merged_aa_DB_clu_rep.fasta -v 2')
	os.system('mmseqs createtsv merged_aa_DB merged_aa_DB merged_aa_DB_clu merged_aa_resultsDB_clu.tsv -v 2')

	## run an all-vs-all blast on the reduced amino acid sequence set ##
	# os.system('makeblastdb -in merged_aa_DB_clu_rep.fasta -out merged_aa -dbtype prot -title "ava_db" -parse_seqids')
	check_call(['makeblastdb', '-in', 'merged_aa_DB_clu_rep.fasta', '-out', 'merged_aa', '-dbtype' ,'prot', '-title', "ava_db", '-parse_seqids'], stdout=DEVNULL, stderr=STDOUT)
	blast_ChunkedSeqs(chunksize=UserInput_faa_chunk_size, fasta_to_chunk_prefix='merged_aa',working_directory=os.getcwd(),UserInput_processes=UserInput_processes)
	df_ava = pd.read_csv('merged_aa_avablast.csv',names=['qseqid','sseqid','mismatch', 'positive','gaps', 'ppos','pident','qcovs','evalue','bitscore','qframe','sframe','sstart', 'send', 'slen', 'qstart', 'qend', 'qlen'])
	df_ava = df_ava[ (df_ava['pident'] >= UserInput_blast_filter_params[0]) & (df_ava['qcovs'] >= UserInput_blast_filter_params[1] ) ]

	## before inputting into MCL, need to make sure that all representative sequences have are present in the ava hits table ##
	## This section creates df_interseciton_orfs that contains sequences not picked up by blast and will therefore not be used for MCL ##
	## Those sequences that did not were not a part of the avablast results will be added in after MCL clustering  as their own unique ortholog cluster ##
	df_mmseqsclust = pd.read_csv('merged_aa_resultsDB_clu.tsv',sep='\t',header=None,names=['mmseqs_cluster_representative','orf_id'])
	all_orfs = set(df_mmseqsclust['mmseqs_cluster_representative'].unique().tolist()) #each locus_tag should have at least one mmseqs_cluster_representative
	avablast_orfs = set(df_ava['qseqid'].unique().tolist())
	intersection_orfs = list(all_orfs.difference(avablast_orfs))
	df_intersection_orfs = pd.DataFrame({'ortho_cluster_members':intersection_orfs})
	print('{} sequences not in AVA blast output'.format(len(df_intersection_orfs)))


	## Cluster ava blast results using MCL ##
	cluster_MCLOrtho(df_mcl=df_ava, similarity_scoring_method=UserInput_similarity_scoring_method, inflation_values_list=UserInput_inflation_values)

	## Assign cluster membership: Use MCL results for all sequences that were present in the AVA results ##
	## add missing sequences (if they exist), as their own cluster ID to the results ##
	## Store the results in df_infval, which is the clsuter assignment for each representative sequence, for each inflation value tested
	## df_infval is a long-form dataframe
	df_infval = assign_ClusterMembership(df_unclustered = df_intersection_orfs, inflation_values_list=UserInput_inflation_values)

	## Generate final output final containing mmseqs cluster representative, ortholog cluster representative, over different inflation values, for ALL sequences in merged_aa.faa ##
	df_mmseqsclust['mmseqs_cluster_id'] = df_mmseqsclust.groupby('mmseqs_cluster_representative').ngroup()
	df_infval = df_infval.merge(df_mmseqsclust, left_on='ortho_cluster_members', right_on='mmseqs_cluster_representative')


	## Store full results: All inflation value relationships
	change_ToWorkingDirectory(directory_name = pjoin(UserInput_main_dir,UserInput_output_dir_name))
	conn = sql.connect('seed_results.db')
	df_infval.to_sql(name='full_mcl_results',con=conn, if_exists='replace',index_label='locus_tag',index=False)
	df_infval.to_csv(pjoin('internal_data','full_mcl_results.csv'),index=False)

	## Store the "best" results: currently the max diversity
	selected_i_val = filter_MCLOorthos(df=df_infval, UserInput_ortho_select_method=UserInput_ortho_select_method) 
	df_infval = df_infval[df_infval['i_val'] == str(selected_i_val)]
	df_infval = df_infval.groupby('orf_id').first().reset_index()
	df_infval.to_sql(name='select_mcl_results',con=conn, if_exists='replace',index_label='locus_tag',index=False)
	df_infval.to_csv(pjoin('internal_data','select_mcl_results.csv'),index=False)


	## remove df_meta and df_feat dataframes from memory ##
	del [[df_infval,df_mmseqsclust]]
	gc.collect()

	## close sql connection
	conn.close()

## --------------End script---------------------##














