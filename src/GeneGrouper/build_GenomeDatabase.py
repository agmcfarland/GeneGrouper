
	
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
from gtr_utils import change_ToWorkingDirectory, make_OutputDirectory, merge_ManyFiles, multiprocesssing_Submission, generate_AssemblyId, chunks
import subprocess
from subprocess import DEVNULL, STDOUT, check_call

def extract_GenbankMetadata(input_filename, gbff_raw_dir):
	'''
	'''
	assembly_id = generate_AssemblyId(input_gbff_file=input_filename)
	for x, record in enumerate(SeqIO.parse(pjoin(gbff_raw_dir,input_filename),'gb')):
		if x == 0:
			df_tax = pd.DataFrame({'assembly_id':[assembly_id], 'filename':[input_filename]})
			try:
				df_tax['organism'] = record.annotations['organism']
			except:
				df_tax['organism'] = 'NA'
		break	
	return(df_tax)



def extract_GenbankFeatures(input_filename, gbff_raw_dir, assembly_faa_dir):
	'''
	Large function to process  that extracts all relevant genbank data
	Extracts CDS features
	Writes all CDS to ./assembly folder
	'''
	assembly_id = generate_AssemblyId(input_gbff_file=input_filename)
	store_rows = []
	with open(pjoin(assembly_faa_dir, assembly_id+'.faa'),'w') as outfile:
		n = 0
		for x, record in enumerate(SeqIO.parse(pjoin(gbff_raw_dir,input_filename),'gb')):
			## extract releavant CDS features ##
			for f in record.features:
				if f.type=='CDS':
					try:
						refseq_locus_tag = f.qualifiers['locus_tag'][0]
					except:
						refseq_locus_tag = 'NA'
					try:
						refseq_product = f.qualifiers['product'][0]
					except:
						refseq_product = 'NA'
					try:
						refseq_gene = f.qualifiers['gene'][0]
					except:
						refseq_gene = 'NA'
					try: 
						reseq_translation = f.qualifiers['translation'][0]
					except:
						reseq_translation = 'NA'
					##
					if reseq_translation == 'NA':
						pseudo_check = 'p'
					else:
						pseudo_check = 'g'
					#internal assembly id, internal contig id, internal locus tag, refseq locus tag, cds start, cds end, cds strand, refseq product, refseq gene
					locus_tag_string = str(n).zfill(5) # generates a buffer of up to 5 zeroes so all locus_tag have same length
					locus_tag = assembly_id+'_'+locus_tag_string
					# check that f.location gives the cds coordinates and not the whole contig coordinates
					if len(f.location.parts) > 1:
						start_cds = f.location.parts[1].start
						end_cds = f.location.parts[1].end
					else:
						start_cds = f.location.start.position
						end_cds = f.location.end.position
					# store in list
					store_rows.append(
						[
						assembly_id,
						assembly_id+'_'+str(x), #contig_id
						locus_tag,
						refseq_locus_tag,
						start_cds,
						end_cds,
						f.location.strand,
						refseq_product,
						refseq_gene,
						pseudo_check
						])
					n+=1
					if reseq_translation == 'NA':
						continue
					## Write CDS to file ##
					outfile.write('>{} {}\n'.format(locus_tag, refseq_locus_tag))
					outfile.write('{}\n'.format(reseq_translation))
		## return genbank features
		df_gbf = pd.DataFrame(store_rows,columns=['assembly_id','contig_id','locus_tag','refseq_locus_tag','cds_start','cds_end','strand','refseq_product','refseq_gene','pseudo_check'])
		return(df_gbf)


def make_BlastDatabase(assembly_id, assembly_faa_dir, blast_db_path):
	'''
	'''
	full_faa_path = pjoin(assembly_faa_dir,assembly_id+'.faa')
	full_db_path = pjoin(blast_db_path,assembly_id)
	check_call(['makeblastdb', '-in', full_faa_path, '-out', full_db_path, '-dbtype', 'prot', '-title', '"{}_db"'.format(assembly_id), '-parse_seqids'], stdout=DEVNULL, stderr=STDOUT)



def update_GenomesDatabase():
	'''
	'''
	pass


def check_ForGenomes(genome_inputs_dir):
	'''
	'''

	pass


def build_GenomesDatabase(genomes_list):
	'''
	'''
	pass






## ---- START OF SCRIPT ---- ##
def BuildGenomeDatabase(UserInput_main_dir, UserInput_genome_inputs_dir, UserInput_processes):
	'''
	'''

	# make directory, with subdirectories to store data. and database for genomic data
	make_OutputDirectory(new_directory=UserInput_main_dir)
	change_ToWorkingDirectory(directory_name=UserInput_main_dir)
	make_OutputDirectory(new_directory='blast_database')
	make_OutputDirectory(new_directory='assemblies')
	conn = sql.connect(pjoin(UserInput_main_dir,'genomes.db')) # genome database holds all base info


	##---- Extract organism info ----- ##
	if os.path.isfile('genomes.db'):
		try:
			# if df_meta exists, then check for which files to be processed are not present in df_meta, if any
			df_meta = pd.read_sql_query("SELECT assembly_id from gb_metadata",conn)
			meta_exists_list = df_meta.assembly_id.tolist()
			r_params = [[f, UserInput_genome_inputs_dir] for f in os.listdir(UserInput_genome_inputs_dir) if f.endswith('.gbff') if generate_AssemblyId(input_gbff_file=f) not in meta_exists_list]
		except:
			# if there is an error reading in df_meta, then it doesn't exist. Therefore, process all input files 
			r_params = [[f, UserInput_genome_inputs_dir] for f in os.listdir(UserInput_genome_inputs_dir) if f.endswith('.gbff')]



	if len(r_params) != 0:
		chunk_r_params_input = chunks(l=r_params, n=10)
		print('Extracting metadata for {} genomes'.format(len(r_params)))
		start_t = time.time()
		df_meta = pd.DataFrame()
		with Pool(processes=UserInput_processes) as p:
			df_meta = pd.concat(p.starmap(extract_GenbankMetadata, r_params[:]))
		df_meta.to_sql(name='gb_metadata',con=conn, if_exists='append',index_label='assembly_id',index=False)
		print((time.time()-start_t)/60)

	else:
		print('All files already have metadata extracted')



	##---- Extract gene features ----- ##
	if os.path.isfile('genomes.db'):
		try:
			# if features have already been pre-processed, then check which files, if any, need to be processed
			df_feat = pd.read_sql_query("SELECT DISTINCT assembly_id from gb_features",conn)
			feat_exist_list = df_feat.assembly_id.tolist()
			r_params = [[f, UserInput_genome_inputs_dir, 'assemblies'] for f in os.listdir(UserInput_genome_inputs_dir) if f.endswith('.gbff') if generate_AssemblyId(input_gbff_file=f) not in meta_exists_list]
		except:
			# if there is an error read in df_feat, then it doesn't exist and all files need to be processed
			r_params = [[f, UserInput_genome_inputs_dir, 'assemblies'] for f in os.listdir(UserInput_genome_inputs_dir) if f.endswith('.gbff')]


	update_sql_index = False
	if len(r_params) != 0:

		print('Extracting features for {} genomes'.format(len(r_params)))
		start_t = time.time()
		chunk_r_params_input = chunks(l=r_params, n=100) # Use a chunking strategy to reduce the memory consumption caused by having large pandas dataframes
		for chunk_r in chunk_r_params_input:
			df_feat = pd.DataFrame()
			with Pool(processes=UserInput_processes) as p:
				df_feat = pd.concat(p.starmap(extract_GenbankFeatures, chunk_r)) 
			df_feat.to_sql(name='gb_features',con=conn, if_exists='append',index_label='locus_tag',index=False)
			print('mem usage:', df_feat.memory_usage(deep=True).sum()/1028/1028)
		print((time.time()-start_t)/60)
		update_sql_index = True

	else:
		print('All files already have features extracted')


	## remove df_meta and df_feat dataframes from memory ##
	del [[df_meta,df_feat]]
	gc.collect()
	df_meta=pd.DataFrame()
	df_feat=pd.DataFrame()


	## index gb_features from genome database on assembly_id to speed up searches ##
	if update_sql_index == True:
		try:
			c = conn.cursor()
			c.execute('CREATE INDEX assembly_id_index on gb_features(assembly_id)')
			conn.commit()
			print('Making assembly_id into index for gb_features')
		except:
			print('Did not update gb_features index')
			pass


	##---- Write blast database files ----- ##
	#
	#check which genomes already have blast databases constructed by comparing blast file to files in df_feat
	df_feat = pd.read_sql_query("SELECT DISTINCT assembly_id from gb_features",conn) 
	feat_exist_list = df_feat.assembly_id.tolist()
	blastdb_files = [x[:-4] for x in os.listdir('blast_database') if x.endswith('.pog')]
	files_to_blastdb = [b for b in feat_exist_list if b not in blastdb_files ]


	if len(files_to_blastdb) != 0:

		print('Making blast database for {} genomes'.format(len(r_params)))
		r_params = [[b,'assemblies','blast_database'] for b in files_to_blastdb]
		start_t = time.time()
		with Pool(processes=UserInput_processes) as p:
			p.starmap(make_BlastDatabase, r_params[:])
		print((time.time()-start_t)/60)
	else:
		print('All files already have blast databases constructed')


	conn.close()

	del conn

## ---- END OF SCRIPT ---- ##


