

import os
from os.path import join as pjoin
import multiprocessing as mp
from multiprocessing import Pool
import re

## -------------- Multi-use functions---------------------##

def change_ToWorkingDirectory(directory_name):
	'''
	General function for changing directory to OUTPUT_FILES_DIR
	'''
	os.chdir(directory_name)
	print('Changing working directory to \n{}\n'.format(directory_name))


def make_OutputDirectory(new_directory):
	'''
	General function for making new output dir
	'''
	if os.path.exists(new_directory) == False:
		print('Making {} directory'.format(new_directory))
		os.mkdir(new_directory)
	else:
		print('{} directory already exists.'.format(new_directory))

def merge_ManyFiles(input_filepath,output_filepath,wildcard_search,output_filename):
	'''
	merge thousands of files. flexible for input/output directory 
	'''
	print('Merging all files with pattern {} and outputting to filename {}'.format(wildcard_search,output_filename))
	input_wildcard_search_path = pjoin(input_filepath,wildcard_search)
	output_filename_path = pjoin(output_filepath,output_filename)
	system_command = 'find {} -type f -exec cat {{}} \; > {}'.format(input_wildcard_search_path,output_filename_path)
	os.system(system_command)

def multiprocesssing_Submission(function_name,submission_list,processor_n,output_type):
	'''
	General function for multiprocessing
	'''
	pool = Pool(processes = processor_n)
	parallel_output = pool.starmap(function_name,submission_list)
	pool.close()
	if output_type == 'array':
		return(parallel_output)

def generate_AssemblyId(input_gbff_file):
	'''
	Take in a gbff filename and split it up to the second underscore. This becomes the assembly id
	'''
	assembly_id_long = re.split('_|\.',input_gbff_file)
	assembly_id = assembly_id_long[0]+'_'+assembly_id_long[1]
	return(assembly_id)


def chunks(l, n):
	'''
	Split a list into evenly sized chunks
	'''
	n = max(1, n)
	return ([l[i:i+n] for i in range(0, len(l), n)])