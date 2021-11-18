
import os
from os.path import join as pjoin
import re
import multiprocessing as mp
from multiprocessing import Pool
from Bio.Seq import Seq
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
import time
import sqlite3 as sql
from collections import defaultdict
import gc
import sys
from gtr_utils import change_ToWorkingDirectory, make_OutputDirectory, merge_ManyFiles, multiprocesssing_Submission, generate_AssemblyId, chunks



def blast_Seed(assembly_id, query_file, blast_db_path, blast_hits_path):
	'''
	'''
	full_db_path = pjoin(blast_db_path, assembly_id)
	full_hits_path = pjoin(blast_hits_path, assembly_id+'.csv')
	os.system('blastp -query {} -db {} -max_target_seqs 100 -evalue 1e-6 -outfmt "10 qseqid sseqid mismatch positive gaps ppos pident qcovs evalue bitscore qframe sframe sstart send slen qstart qend qlen" -num_threads 1 -out {}'.format(query_file, full_db_path, full_hits_path))

## write seed region gene content to file ##

def test_RegionOverlap(x1,x2,y1,y2):
	'''
	Asks: is the largest of the smallest distances less than or equal to the smallest of the largest distances?
	Returns True or False
	'''
	return( max(x1,y1) <= min(x2,y2) )



def filter_BlastHits(df, identity_cutoff, coverage_cutoff):
	'''
	Remove hits that do not meet identity and coverage cutoffs
	'''
	df = df[ (df['pident'] >= identity_cutoff ) & (df['qcovs'] >= coverage_cutoff) ]
	return(df)



def extract_SeedRegions(assembly_id, upstream_search_length, downstream_search_length, identity_cutoff, coverage_cutoff, hits_cutoff):
	'''
	Check for overlaps in seed hits. When this occurs, pick the best hit for overlap.
	Extract x basepairs upstream and downstream of a seed blast hit.
	Write to sql database under table 'seed_regions'
	'''

#--------- Troubleshooting --------#
# pd.set_option('display.max_columns', None)
# upstream_search_length, downstream_search_length, identity_cutoff, coverage_cutoff, hits_cutoff = 10000,10000,20,80, 1

# # UserInput_main_dir = '/projects/b1042/HartmannLab/alex/GeneGrouper_test/testbed/dataset1/test1'#'/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/syntenease_project/gtr/testbed/dataset1/test1'
# UserInput_main_dir = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/syntenease_project/gtr/testbed/dataset3/test1'
# UserInput_output_dir_name = pjoin(UserInput_main_dir,'pdua')
# os.chdir(UserInput_main_dir)
# conn = sql.connect(pjoin(UserInput_main_dir,'genomes.db')) # genome database holds all base info
# conn2 = sql.connect(pjoin(UserInput_output_dir_name,'seed_results.db')) # seed_results database holds all seed search specific info. 


# assembly_id = 'GCF_009800085'  #GCF_009800085_02009     dbscan_label 3     ,global_strand -1 has more
# assembly_id = 'GCF_000583895' #GCF_000583895_02788     dbscan_label 13    ,global_strand 1  has fewer

# assembly_id = 'GCF_001077835'

# assembly_id = 'GCF_000251025'

#--------- Troubleshooting --------#

	try:
		df_hits = pd.read_csv(pjoin('temp_blast_results',assembly_id+'.csv'),names=['qseqid','sseqid','mismatch', 'positive','gaps', 'ppos','pident','qcovs','evalue','bitscore','qframe','sframe','sstart', 'send', 'slen', 'qstart', 'qend', 'qlen'])
	except:
		return(pd.DataFrame())
		pass

	## filter hit results, exit function if not hits pass cutoffs
	df_hits = filter_BlastHits(df=df_hits, identity_cutoff=identity_cutoff, coverage_cutoff=coverage_cutoff)
	if len(df_hits) == 0:
		return(pd.DataFrame())
		pass


	## read in whole genome info
	df_g = pd.read_sql_query("SELECT contig_id, locus_tag, cds_start, cds_end, strand, pseudo_check from gb_features WHERE assembly_id = '{}' ".format(assembly_id), sql.connect('genomes.db'))


	## merge genome level data with blast hit results. sort by identity and evalue so that best hits are at the top of the table
	df_hits = df_g.merge(df_hits[['sseqid','pident','qcovs','evalue']],left_on='locus_tag',right_on='sseqid',how='inner')


	# if there are multiple hits that are identical, keep the one with the highest evalue
	df_hits = df_hits.groupby('locus_tag').first().reset_index()

	df_hits = df_hits.sort_values(['pident','evalue'],ascending=[False,False])

	# keep only n total hits (default is all hits)
	if type(hits_cutoff) == int:
		df_hits = df_hits.iloc[:hits_cutoff, :]


	## define upstream and downstream search lengths. The distance inputs for upstream/start and downstream/end are switched when the seed gene is in the -1 frame
	## The usl and dsl are not strand specific now. They are relative to the genomes' upstream and downstream positions
	df_hits['usl'] = np.where(df_hits['strand']==1, df_hits['cds_start'] - upstream_search_length, df_hits['cds_start'] - downstream_search_length)
	df_hits['dsl'] = np.where(df_hits['strand']==1, df_hits['cds_end'] + downstream_search_length, df_hits['cds_end'] + upstream_search_length)

	## for each contig, test for overlapping region ranges. Keep the region with the best hit.
	## append subsets to new df called df_hits_parsed
	df_hits_parsed = pd.DataFrame()
	for cid in df_hits['contig_id'].unique().tolist():
		#Looping over contigs to prevent instances where they overlap but are not on the same contig.
		df_hits_contig = df_hits[df_hits['contig_id']==cid]
		for hid in df_hits_contig['locus_tag']:
			try:
				# create range to compare all other ranges to			
				hid_x1, hid_x2 = df_hits_contig[df_hits_contig['contig_id']==cid][df_hits_contig['locus_tag']==hid]['usl'].item(), df_hits_contig[df_hits_contig['contig_id']==cid][df_hits_contig['locus_tag']==hid]['dsl'].item()
				# test for overlap with all other ranges on the same contig
				df_hits_contig['overlap'] = df_hits_contig[df_hits_contig['contig_id']==cid].apply(lambda x: test_RegionOverlap(x1=hid_x1, x2=hid_x2, y1=x['usl'], y2=x['dsl']), axis=1)
				keep_index  = (df_hits_contig['overlap']==True).idxmax()
				df_hits_contig['overlap_representative'] = np.where(df_hits_contig.index==keep_index,True,False)
				# remove rows that are not repesentatives but do overlap. keep everything else.
				df_hits_contig = df_hits_contig[ (df_hits_contig['overlap']==False) | (df_hits_contig['overlap_representative'] == True) | (df_hits_contig['overlap']).isnull() == True ]
			except:
				continue
		df_hits_parsed = df_hits_parsed.append(df_hits_contig)



	## build a dataframe that contains genes upsteam and downstream of boundaries for each seed blast hit ##
	df_rkeep = pd.DataFrame()
	for h in df_hits_parsed.iterrows():


		# get coordinates for upstream and downstream gene boundaries

		# ----- attempt start ----- #
		# cds_start for strand == 1 
		# cds_end for strand == -1
		if h[1][4] == 1:
			cds_search_list = df_g[df_g['contig_id'] == h[1][1]]['cds_start'].tolist()
			usl_locus_tag_bp = min(cds_search_list, key = lambda x: abs(x-h[1][10]))

			cds_search_list = df_g[df_g['contig_id'] == h[1][1]]['cds_end'].tolist()
			dsl_locus_tag_bp = min(cds_search_list, key = lambda x: abs(x-h[1][11]))

		if h[1][4] == -1:
			cds_search_list = df_g[df_g['contig_id'] == h[1][1]]['cds_start'].tolist()
			usl_locus_tag_bp = min(cds_search_list, key = lambda x: abs(x-h[1][10]))

			cds_search_list = df_g[df_g['contig_id'] == h[1][1]]['cds_end'].tolist()
			dsl_locus_tag_bp = min(cds_search_list, key = lambda x: abs(x-h[1][11]))



		# subset the main genome dataframe to contain only ORFs that are within the designated bounds
		if h[1][4] == 1:	
			df_r = df_g[ df_g['contig_id'] == h[1][1] ][ ( df_g['cds_end'] >= usl_locus_tag_bp) & (df_g['cds_start'] <= dsl_locus_tag_bp ) ]
		if h[1][4] == -1:	
			df_r = df_g[ df_g['contig_id'] == h[1][1] ][ ( df_g['cds_end'] >= usl_locus_tag_bp) & (df_g['cds_start'] <= dsl_locus_tag_bp ) ]
		# ----- attempt end ----- #

		# append the subsetted dataframe to a new dataframe that will contain all region extractions
		df_r['region_id'] = h[1][0]
		df_rkeep = df_rkeep.append(df_r,ignore_index=True)

	## add blast data to the extracted regions, remove unneeded columns, and return the dataframe
	df_rkeep['assembly_id'] = assembly_id
	df_rkeep = df_rkeep[['region_id','assembly_id','contig_id','locus_tag','strand','pseudo_check']]
	df_rkeep = df_rkeep.merge(df_hits_parsed[['sseqid','pident','qcovs','evalue']], left_on='region_id',right_on='sseqid')
	df_rkeep = df_rkeep.drop(columns='sseqid')


	return(df_rkeep)


def write_RegionSeqsToFile(assembly_id, output_dir_name):
	'''
	Write a .faa containing all sequences that were found to be within a defined seed region for a given assembly_id
	
	'''
	df_s =  pd.read_sql_query("SELECT locus_tag from seed_regions WHERE assembly_id = '{}'".format(assembly_id), conn2)
	if len(df_s) == 0:
		return()

	seqs_to_write = df_s['locus_tag'].tolist()
	with open(pjoin(output_dir_name,'regions',assembly_id+'.faa'),'w') as outfile:
		for record in SeqIO.parse(pjoin('assemblies',assembly_id+'.faa'),'fasta'):
			if record.id in seqs_to_write:
				outfile.write('>{}\n'.format(record.id))
				outfile.write('{}\n'.format(str(record.seq)))




def BuildQueryFiles(
	UserInput_main_dir,
	UserInput_genome_inputs_dir,
	UserInput_output_dir_name,
	UserInput_query_filename_path,
	UserInput_upstream_search_length,
	UserInput_downstream_search_length,
	UserInput_identity_threshold,
	UserInput_coverage_threshold,
	UserInput_hitcount_threshold,
	UserInput_processes,
	):
	'''
	'''

#--------- Troubleshooting --------#
# pd.set_option('display.max_columns', None)
# mp.set_start_method("fork")
# # UserInput_main_dir = '/projects/b1042/HartmannLab/alex/GeneGrouper_test/testbed/dataset1/test1'#'/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/syntenease_project/gtr/testbed/dataset1/test1'
# UserInput_main_dir = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/syntenease_project/gtr/testbed/dataset3/test1'
# UserInput_output_dir_name = pjoin(UserInput_main_dir,'pdua')
# UserInput_genome_inputs_dir = '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/syntenease_project/gtr/testbed/dataset3/core'
# UserInput_query_filename_path = pjoin(UserInput_genome_inputs_dir,'pdua.txt')
# UserInput_upstream_search_length = 2000
# UserInput_downstream_search_length = 18000
# UserInput_identity_threshold = 15
# UserInput_coverage_threshold = 70
# UserInput_hitcount_threshold = 5
# UserInput_processes = 8
# os.chdir(UserInput_main_dir)
# # conn = sql.connect(pjoin(UserInput_main_dir,'genomes.db')) # genome database holds all base info
# # conn2 = sql.connect(pjoin(UserInput_output_dir_name,'seed_results.db')) # seed_results database holds all seed search specific info. 
#--------- Troubleshooting --------#


	make_OutputDirectory(new_directory=UserInput_main_dir)
	change_ToWorkingDirectory(directory_name=UserInput_main_dir)
	make_OutputDirectory(new_directory=UserInput_output_dir_name)
	make_OutputDirectory(new_directory=pjoin(UserInput_output_dir_name,'regions'))
	make_OutputDirectory(new_directory=pjoin(UserInput_output_dir_name,'ortholog_clusters'))
	make_OutputDirectory(new_directory=pjoin(UserInput_output_dir_name,'internal_data'))
	make_OutputDirectory(new_directory=pjoin(UserInput_output_dir_name,'subgroups'))
	make_OutputDirectory(new_directory=pjoin(UserInput_output_dir_name,'visualizations'))
	make_OutputDirectory(new_directory='temp_blast_results')

	global conn, conn2

	conn = sql.connect(pjoin(UserInput_main_dir,'genomes.db')) # genome database holds all base info
	conn2 = sql.connect(pjoin(UserInput_output_dir_name,'seed_results.db')) # seed_results database holds all seed search specific info. 


	##---- Blast seed against database ----- ##
	r_params = [[generate_AssemblyId(input_gbff_file=f), UserInput_query_filename_path, 'blast_database', 'temp_blast_results'] for f in os.listdir(UserInput_genome_inputs_dir) if f.endswith('.gbff')]
	print('Blasting seed against {} genomes'.format(len(r_params)))
	start_t = time.time()
	with Pool(processes=UserInput_processes) as p:
		p.starmap(blast_Seed, r_params[:])
	print((time.time()-start_t)/60)


	##---- identify and extract seed regions ----- ##
	r_params = [[generate_AssemblyId(input_gbff_file=f), UserInput_upstream_search_length, UserInput_downstream_search_length,UserInput_identity_threshold,UserInput_coverage_threshold, UserInput_hitcount_threshold] for f in os.listdir(UserInput_genome_inputs_dir) if f.endswith('.gbff')]
	print('Identifying and extracting seed regions for {} genomes'.format(len(r_params)))
	start_t = time.time()
	df_regions = pd.DataFrame()
	with Pool(processes=UserInput_processes) as p:
		df_regions = pd.concat(p.starmap(extract_SeedRegions, r_params[:]))
	print((time.time()-start_t)/60)
	print('mem usage:', df_regions.memory_usage(deep=True).sum()/1028/1028)

	if len(df_regions) == 0:
		raise ValueError('No genes matched your query!')
		
	df_regions.to_sql(name='seed_regions',con=conn2, if_exists='replace',index_label='locus_tag',index=False)


	##---- index seed_regions from seed_results databse on assembly_id to speed up searches ----- ##
	c2 = conn2.cursor()
	c2.execute('CREATE INDEX assembly_id_index on seed_regions(assembly_id)') # this is key to speeding up search
	conn2.commit()


	##---- write seed region gene content to sql database ----- ##
	r_params = [[generate_AssemblyId(input_gbff_file=f), UserInput_output_dir_name] for f in os.listdir(UserInput_genome_inputs_dir) if f.endswith('.gbff')]
	print('Writing seed region sequences for for {} genomes'.format(len(r_params)))
	start_t = time.time()
	df_regions = pd.DataFrame()
	with Pool(processes=UserInput_processes) as p:
		p.starmap(write_RegionSeqsToFile, r_params[:])
	print((time.time()-start_t)/60)

	## merge sequence files and write to ortholog_clusters df
	merge_ManyFiles(input_filepath=pjoin(UserInput_output_dir_name,'regions'),output_filepath=pjoin(UserInput_output_dir_name,'ortholog_clusters'),wildcard_search='*.faa',output_filename='merged_aa.faa')

	## delete ./temp_blast_results and <seed_region>/regions
	os.system('rm -R ./temp_blast_results')
	#### os.system('rm -R ./{}/regions'.format(UserInput_output_dir_name))


	conn.close()
	conn2.close()

	del conn
	del conn2


## ---- END OF SCRIPT ---- ##
