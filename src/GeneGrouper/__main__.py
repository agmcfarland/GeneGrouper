#!/usr/bin/env python

import warnings
warnings.simplefilter(action='ignore', category=UserWarning)
	
import os
from os.path import join as pjoin
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
	print('run_BuildGenomeDatabase with {} threads\n'.format(args.threads))
	MainDir = os.path.abspath(args.project_directory)
	GenomeInputsDir = os.path.abspath(args.genomes_directory)

	if os.path.isdir(MainDir) == False:
		try:
			os.mkdir(MainDir)
		except:
			raise ValueError('main directory not found')
	if os.path.isdir(GenomeInputsDir) == False:
		raise ValueError('genomes directory not found')

	start_t = time.time()
	BuildGenomeDatabase(
	UserInput_genome_inputs_dir=GenomeInputsDir, 
	UserInput_main_dir=MainDir, 
	UserInput_processes=args.threads
	)
	print('\nFinished building database\nTotal time (m): {}'.format((time.time()-start_t)/60))

def run_RegionSearch(args):
	'''
	'''
	print('run_RegionSearch with {} threads\n'.format(args.threads))


	MainDir = os.path.abspath(args.project_directory)
	GenomeInputsDir = os.path.abspath(args.genomes_directory)
	QueryFilenamePath = os.path.abspath(args.query_file)

	if os.path.isdir(MainDir) == False:
		raise ValueError('main directory not found')
	if os.path.isdir(GenomeInputsDir) == False:
		raise ValueError('genomes directory not found')
	if os.path.isfile(QueryFilenamePath) == False:
		raise ValueError('seed file not found')
	if os.path.isdir(pjoin(MainDir,args.search_name)) == True:
		if args.force == True:
			os.system('rm -R {}'.format(pjoin(MainDir,args.search_name)))
		else:
			raise ValueError('search directory already exists')
			print('error!!!')

	start_t = time.time()
	BuildQueryFiles(
		UserInput_main_dir = MainDir,#args.project_directory, 
		UserInput_genome_inputs_dir= GenomeInputsDir,#args.genomes_directory,
		UserInput_output_dir_name = args.search_name,
		UserInput_query_filename_path = QueryFilenamePath,
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
		UserInput_reclustering_iterations = args.recluster_iterations,
		UserInput_min_group_size = args.min_group_size
		)

	OutputRegionResults(
		UserInput_main_dir = MainDir,#args.project_directory, 
		UserInput_output_dir_name = args.search_name,
		)

	print('\nFinished searching for gene clusters\nTotal time (m): {}'.format((time.time()-start_t)/60))


def run_MakeVisualizations(args):
	print('make_Visualizations\n')


	MainDir = os.path.abspath(args.project_directory)
	if os.path.isdir(MainDir) == False:
		raise ValueError('main directory not found')

	if args.visual_type == 'main':
		print('Generating main visualization')
	if args.visual_type == 'region_cluster':
		print('Generating region cluster visualization for cluster label {}'.format(args.cluster_label))

	start_t = time.time()
	MakeVisualizations(
		UserInput_main_dir = MainDir,
		UserInput_output_dir_name = args.search_name,
		UserInput_cluster_label_id = args.group_label,
		UserInput_visualization_type = args.visual_type,
		UserInput_script_path = pjoin(script_path,'Rscripts'),
		UserInput_image_format = args.image_format,
		UserInput_tip_label_type = args.tip_label_type,
		UserInput_tip_label_size = args.tip_label_size
		)

	print('\nFinished making visualization\nTotal time (m): {}'.format((time.time()-start_t)/60))



def main(args = None):
	'''
	'''
	if args is None:
		args = sys.argv[1:]

	# create parser
	parser = argparse.ArgumentParser(prog = 'GeneGrouper')

	# add common args to parser
	parser.add_argument('-d', '--project_directory', type=str, default=os.getcwd(), help = 'Main directory to contain the base files used for region searching and clustering. Default current directory.', metavar='')
	parser.add_argument('-n', '--search_name', type=str, default = 'region_search', help = 'Name of the directory to contain search-specific results. Default region_search', metavar='')
	parser.add_argument('-g', '--genomes_directory', type=str, default=pjoin(os.getcwd(),'genomes'), help = 'Directory containing genbank-file format genomes with the suffix .gbff. Default ./genomes.', metavar='')
	parser.add_argument('-t', '--threads', type=int, default= mp.cpu_count(), help = 'Number of threads to use. Default all threads.', metavar='')

	# make subparser
	subparsers = parser.add_subparsers(title='subcommands', description ='Valid subcommands', help = 'sub-command help', dest='command')

	# build_database subparser 
	parser_bd = subparsers.add_parser('build_database', help ='Convert a set of genomes into a useable format for GeneGrouper')
	parser_bd.set_defaults(func=run_BuildGenomeDatabase)

	# find_regions args
	parser_fr = subparsers.add_parser('find_regions', help ='Find regions given a translated gene and a set of genomes')
	parser_fr.add_argument('-f', '--query_file', type=str, required = True, help = 'Provide the absolute path to a fasta file containing a translated gene sequence.', metavar='')
	parser_fr.add_argument('-us', '--upstream_search', type=int, default=10000, help = 'Upstream search length in basepairs. Default=10000', metavar='')
	parser_fr.add_argument('-ds', '--downstream_search', type=int, default=10000, help = 'Downstream search length in basepairs. Default=10000', metavar='')
	parser_fr.add_argument('-i', '--seed_identity', type=int, default=60, help = 'Identity cutoff for initial blast search. Default=60', metavar='')
	parser_fr.add_argument('-c', '--seed_coverage', type=int, default=90, help = 'Coverage cutoff for initial blast search. Default=90', metavar='')
	parser_fr.add_argument('-hk', '--seed_hits_kept', type = int, default = None, help = 'Number of blast hits to keep. Default=None', metavar='')
	parser_fr.add_argument('--min_group_size', default = 'default', help = 'The minimum number of gene regions to constitute a group. Default=ln(jaccard distance length)', metavar='')
	parser_fr.add_argument('-re', '--recluster_iterations', type = int, default = 0, help = 'Number of iterations to run DBSCAN on remaining ungrouped regions . Default=0', metavar='')
	parser_fr.add_argument('--force', action='store_true', help = 'Flag to overwrite search name directory.')
	parser_fr.set_defaults(func=run_RegionSearch)

	# visualize_regions args
	parser_mv = subparsers.add_parser('visualize', help ='Visualize GeneGrouper outputs. Three visualization options are provided. Check the --visual_type help description.')
	parser_mv.add_argument('--visual_type', type=str, choices = ['main','group','tree'], default = 'main', help="Choices: [main, group, tree]. Use main for main visualizations. Use group to inspect specific group. Use tree for a phylogenetic tree of representative seed sequencess. Default=main", metavar='')
	parser_mv.add_argument('--group_label', type=int, default=-1, help='The integer identifier of the group you wish to inspect. Default=-1',metavar='')
	parser_mv.add_argument('--image_format', type=str, choices = ['png','svg'], default='png', help='Choices: [png,svg]. Output image format. Use svg if you want to edit the images. Default=png.',metavar='')
	parser_mv.add_argument('--tip_label_type', type=str, choices = ['full','group'], default='full', help='Choices: [full, group]. Use full to include the sequence ID followed by group ID. Use group to only have the group ID. Default=full',metavar='')
	parser_mv.add_argument('--tip_label_size', type=int, default=2, help='Specify the tip label size in the output image. Default=2',metavar='')

	parser_mv.set_defaults(func=run_MakeVisualizations)

	args = parser.parse_args()
	args.func(args)


## -------------- Run ---------------------##

if __name__ == '__main__':
	sys.exit(main())

#os.system('conda list --export > list.txt')


