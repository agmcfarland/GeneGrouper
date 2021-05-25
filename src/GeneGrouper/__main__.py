#!/usr/bin/env python
	
import os
from os.path import join as pjoin
import warnings
import time
import argparse
import sys
import multiprocessing as mp
from multiprocessing import Pool


mp.set_start_method("fork")
global script_path

script_path = os.path.split(os.path.realpath(os.path.abspath(__file__)))[0]
sys.path.append(script_path)

from gtr_utils import change_ToWorkingDirectory, make_OutputDirectory, merge_ManyFiles, multiprocesssing_Submission, generate_AssemblyId, chunks
from build_GenomeDatabase import *
from build_QueryFiles import *
from cluster_RegionSeqs import *
from cluster_Regions import *
from output_RegionResults import *
from make_Visualizations import *



def run_BuildGenomeDatabase(args):
	'''
	'''
	print('run_BuildGenomeDatabase\n')

	MainDir = os.path.abspath(args.project_directory)
	GenomeInputsDir = os.path.abspath(args.genomes_directory)

	start_t = time.time()
	BuildGenomeDatabase(
	UserInput_genome_inputs_dir=GenomeInputsDir, 
	UserInput_main_dir=MainDir, 
	UserInput_processes=args.threads
	)
	print('\nTotal time (m): {}'.format((time.time()-start_t)/60))

def run_RegionSearch(args):
	'''
	'''
	print('run_RegionSearch\n')


	MainDir = os.path.abspath(args.project_directory)
	GenomeInputsDir = os.path.abspath(args.genomes_directory)
	SeedFilenamePath = os.path.abspath(args.seed_file)

	if os.path.isdir(MainDir) == False:
		raise ValueError('main directory not found')
	if os.path.isdir(GenomeInputsDir) == False:
		raise ValueError('genomes directory not found')
	if os.path.isfile(SeedFilenamePath) == False:
		raise ValueError('seed file not found')

	start_t = time.time()
	BuildQueryFiles(
		UserInput_main_dir = MainDir,#args.project_directory, 
		UserInput_genome_inputs_dir= GenomeInputsDir,#args.genomes_directory,
		UserInput_output_dir_name = args.search_name,
		UserInput_seed_filename_path = SeedFilenamePath,
		UserInput_upstream_search_length = args.upstream_search,
		UserInput_downstream_search_length = args.downstream_search,
		UserInput_identity_threshold = args.seed_identity,
		UserInput_coverage_threshold = args.seed_coverage,
		UserInput_hitcount_threshold = args.seed_hits_kept,
		UserInput_processes = args.threads
		)

	ClusterRegionSeqs(
		UserInput_main_dir = MainDir,#args.project_directory, 
		UserInput_output_dir_name = args.search_name,
		UserInput_faa_chunk_size = 10,
		UserInput_blast_filter_params = [50,95], #identity, qcoverage
		UserInput_inflation_values = ['1.2','1.4','2.0','4.0','6.0'],
		UserInput_similarity_scoring_method = 'evalue',
		UserInput_linclust_settings = '--min-seq-id 0.95',
		UserInput_ortho_select_method = 'max_div',
		UserInput_processes = args.threads
		)

	ClusterRegions(
		UserInput_main_dir = MainDir,#args.project_directory, 
		UserInput_output_dir_name = args.search_name,
		UserInput_reclustering_iterations = args.recluster_iterations
		)

	OutputRegionResults(
		UserInput_main_dir = MainDir,#args.project_directory, 
		UserInput_output_dir_name = args.search_name,
		)

	print('\nTotal time (m): {}'.format((time.time()-start_t)/60))


def run_MakeVisualizations(args):
	print('make_Visualizations\n')


	MainDir = os.path.abspath(args.project_directory)
	if os.path.isdir(MainDir) == False:
		raise ValueError('main directory not found')

	if args.visual_type == 'main':
		print('MESSAGE: Generating main visualization')
	if args.visual_type == 'region_cluster':
		print('MESSAGE: Generating region cluster visualization for cluster label {}'.format(args.cluster_label))

	start_t = time.time()
	MakeVisualizations(
		UserInput_main_dir = MainDir,
		UserInput_output_dir_name = args.search_name,
		UserInput_cluster_label_id = args.cluster_label,
		UserInput_visualization_type = args.visual_type,
		UserInput_script_path = pjoin(script_path,'Rscripts')
		)

	print('\nTotal time (m): {}'.format((time.time()-start_t)/60))




def main(args = None):
	'''
	'''
	if args is None:
		args = sys.argv[1:]

	# create parser
	parser = argparse.ArgumentParser(prog = 'GeneToRegions')

	# add common args to parser
	parser.add_argument('-d', '--project_directory', type=str, default=os.getcwd(), help = 'main directory to contain the base files used for region searching and clustering')
	parser.add_argument('-n', '--search_name', type=str, default = 'region_search', help = 'name of the directory to contain search-specific results')
	parser.add_argument('-g', '--genomes_directory', type=str, default=pjoin(os.getcwd(),'genomes'), help = 'directory containing genbank-file format genomes with the suffix .gbff')
	parser.add_argument('-t', '--threads', type=int, default= mp.cpu_count(), help = 'number of threads to use')

	# make subparser
	subparsers = parser.add_subparsers(title='subcommands', description ='valid subcommands', help = 'sub-command help', dest='command')

	# build_database subparser 
	parser_bd = subparsers.add_parser('build_database', help ='convert a set of genomes into a useuable format for GeneToRegions ')
	parser_bd.set_defaults(func=run_BuildGenomeDatabase)

	# find_regions args
	parser_fr = subparsers.add_parser('find_regions', help ='find regions given a translated gene and a set of genomes')
	parser_fr.add_argument('-f', '--seed_file', type=str, required = True, help = 'provide the absolute path to a fasta file containing a translated gene sequence')
	parser_fr.add_argument('-us', '--upstream_search', type=int, default=10000, help = 'upstream search length in basepairs')
	parser_fr.add_argument('-ds', '--downstream_search', type=int, default=10000, help = 'downstream search length in basepairs')
	parser_fr.add_argument('-i', '--seed_identity', type=int, default=50, help = 'identity cutoff for initial blast search')
	parser_fr.add_argument('-c', '--seed_coverage', type=int, default=90, help = 'coverage cutoff for initial blast search')
	parser_fr.add_argument('-hk', '--seed_hits_kept', type = int, default = None, help = 'number of blast hits to keep')
	parser_fr.add_argument('-re', '--recluster_iterations', type = int, default = 0, help = 'number of region re-clustering attempts after the initial clustering')
	parser_fr.set_defaults(func=run_RegionSearch)

	# visualize_regions args
	parser_mv = subparsers.add_parser('visualize', help ='visualize region clusters')
	parser_mv.add_argument('-vt', '--visual_type', type=str, choices = ['main','region_cluster'], default = 'main')
	parser_mv.add_argument('-clab', '--cluster_label', type=int, default=-1)
	parser_mv.set_defaults(func=run_MakeVisualizations)

	args = parser.parse_args()
	# print(args)
	args.func(args)


## -------------- Run ---------------------##

if __name__ == '__main__':
	sys.exit(main())


