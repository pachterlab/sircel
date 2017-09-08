"""
Akshay Tambe
Pachter and Doudna groups
"""


from sircel import IO_utils, Plot_utils, Split_reads
import argparse
import os
import sys
import json
import time
import shutil
import subprocess
from pathlib import Path
import gc

from scipy.sparse import coo_matrix
from sklearn.preprocessing import normalize
from sklearn.metrics.pairwise import pairwise_distances
from scipy.spatial.distance import *

def run_all(args):
	print('\nInspecting and pre-processing input data')
	
	if(args['output_dir'][-1] == '/'):
		args['output_dir'] = args['output_dir'][0:-1]
	if not os.path.exists(args['output_dir']):
		os.makedirs(args['output_dir'])
	if not os.path.exists(args['output_dir'] + '/plots'):
		os.makedirs(args['output_dir'] + '/plots')
	with (Path(__file__).parent / 'params.json').open() as r:
		kallisto = json.load(r)['kallisto']
	assert kallisto

	split_args = {}
	check_pipeline_input(args, kallisto)
	
	if args['dropseq']:
		args['barcode_start']	= 0
		args['barcode_end'] 		= 12
		args['umi_start']			= 12
		args['umi_end'] 			= 20	
		if args['kmer_size'] == None:
			args['kmer_size'] = 8
	
		print('Unzipping and indexing files (temporary)')
		reads_unzipped, reads_offsets = \
			IO_utils.unzip(args['reads'].split(','))
		barcodes_unzipped, barcodes_offsets = \
			IO_utils.unzip(args['barcodes'].split(','))
		args['reads'] = reads_unzipped
		args['barcodes'] = barcodes_unzipped
		args['reads_offsets'] = reads_offsets
		args['barcodes_offsets'] = barcodes_offsets
		
	elif args['10xgenomics']:
		args['barcode_start']	= 0
		args['barcode_end'] 		= 26
		args['umi_start'] 		= 26
		args['umi_end'] 			= 34
		if args['kmer_size'] == None:
			args['kmer_size'] = 16		
	
		print('Unzipping and indexing files (temporary)')
		reads_unzipped, reads_offsets = \
			IO_utils.unzip(args['reads'].split(','))
		barcodes_unzipped, barcodes_offsets = IO_utils.merge_barcodefiles_10x(
			args['barcodes'].split(','),
			args['umis'].split(','))
				
		args['reads'] = reads_unzipped
		args['barcodes'] = barcodes_unzipped	
		args['reads_offsets'] = reads_offsets
		args['barcodes_offsets'] = barcodes_offsets	

	else:
		if args['kmer_size'] == None:
			args['kmer_size'] = 8
	
		print('Unzipping and indexing files (temporary)')
		reads_unzipped, reads_offsets = \
			IO_utils.unzip(args['reads'].split(',')[0])
		barcodes_unzipped, barcodes_offsets = \
			IO_utils.unzip(args['barcodes'].split(',')[0])
		args['reads'] = reads_unzipped
		args['barcodes'] = barcodes_unzipped	
		args['reads_offsets'] = reads_offsets
		args['barcodes_offsets'] = barcodes_offsets
	
	assert len(reads_offsets) == len(barcodes_offsets), \
		'Reads file and barcodes file do not have the same number of reads\n%i, %i' % \
		(len(reads_offsets), len(barcodes_offsets))
	
	check_split_input(args)
	reads_nuc_content = \
		IO_utils.get_nuc_content(reads_unzipped, len(reads_offsets))
	barcodes_nuc_content = \
		IO_utils.get_nuc_content(barcodes_unzipped, len(barcodes_offsets))
		
	nuc_content_plt = \
		Plot_utils.plot_nuc_content(
			args['output_dir'], 
			reads_nuc_content,
			barcodes_nuc_content)
	print('\n')
	
	output_files, elapsed_time = Split_reads.run_all(args)
	output_files['args'] = args
	print('Split reads. Time elapsed %s seconds' % elapsed_time)
	
	output_files['nuc_content'] = nuc_content_plt
	
	del(reads_offsets)
	del(barcodes_offsets)
	_ = gc.collect()
	
	#print(args['kallisto_idx'])
	if(args['kallisto_idx'] != None):
		print('Running kallisto')
		kallisto_dir = '%s/kallisto_outputs' % args['output_dir']
		if not os.path.exists(kallisto_dir):
			os.makedirs(kallisto_dir)
		
		output_files['kallisto'] = run_kallisto(
			args,
			kallisto,
			kallisto_dir,
			output_files)
		print('Getting transcript compatibility counts')
		output_files['tcc'] = write_transcript_compatability_counts(
			args,
			output_files,
			kallisto_dir)
				
	print('Removing temp files')
	os.unlink(reads_unzipped)
	os.unlink(barcodes_unzipped)
	
	output_files['run_outputs'] = \
		'%s/run_outputs.json' % args['output_dir']
	with open(output_files['run_outputs'], 'w') as writer:
		writer.write(json.dumps(output_files, indent=3))
	print('Done.')
	return output_files
	
def run_kallisto(args, kallisto, kallisto_dir, output_files):
	kallisto_start_time = time.time()
	
	kallisto_cmd = [
		kallisto, 'pseudo',
		'-b',			output_files['split']['batch'],
		'-i', 		args['kallisto_idx'],
		'-o',			kallisto_dir,
		'-t',			str(args['threads']),
		'--umi']
	kallisto = subprocess.check_call(kallisto_cmd)
	kallisto_output = {
		'tsv' :				'%s/matrix.tsv' % kallisto_dir,
		'run_info'	:		'%s/run_info.json' % kallisto_dir,
		'equiv_classes' :	'%s/matrix.ec' % kallisto_dir,
		'cells' : 			'%s/matrix.cells' % kallisto_dir}
	current_time = time.time()
	elapsed_time = current_time - kallisto_start_time
	print('Time elapsed %0.2f' % elapsed_time)
	
	return kallisto_output

def write_transcript_compatability_counts(args, input_files, kallisto_dir):
	"""
	Modifed from Vasilis Ntranos scRNA-seq-TCC-prep
		https://github.com/pachterlab/scRNA-Seq-TCC-prep
	"""	
	
	import numpy as np
	import pickle
	
	for fname in input_files['kallisto'].values():
		assert(os.path.exists(fname)), \
			print('kallisto output file not found: %s' % fname)
		
	output_files = {
		'tcc_coo':		'%s/tcc_matrix.dat' % kallisto_dir,
		'tcc_csr':		'%s/tcc_matrix_csr.dat' % kallisto_dir,
		'l1_dist' : 	'%s/pairwise_l1_distance.npy' % kallisto_dir,
		'nonzero_eq' :	'%s/nonzero_equiv_classes.npy' % kallisto_dir,
		'tcc_norm_t' :	'%s/tcc_normalized_transposed.npy' % kallisto_dir}
	
	print('\tLoading kallisto matrices')
	#matrix.ec file
	equiv_class_fname =  input_files['kallisto']['equiv_classes']
	tsv_fname = input_files['kallisto']['tsv']
	
	tsv = np.genfromtxt(tsv_fname, delimiter='\t' , dtype=int)
	tsv_rows, tsv_cols, tsv_data = tsv.T
	nonzero_equiv_classes = np.unique(tsv_rows)
	
	map_rows = {
		val:ind for ind,val in enumerate(nonzero_equiv_classes)}
	map_cols = {
		val:ind for ind,val in enumerate(np.unique(tsv_cols))}
	
	print('\tPreparing TCC matrix')
	#use scipy coo_matrix for sparse data
	tcc_matrix = coo_matrix(
		(tsv_data.astype(float),
		( [map_rows[r] for r in tsv_rows], [map_cols[c] for c in tsv_cols])))	
	tcc_matrix_csr = tcc_matrix.tocsr()
	#save coo_matrix and csr_matrix using pickle (?)
	with open(output_files['tcc_coo'], 'wb') as outf:
		pickle.dump(tcc_matrix, outf, pickle.HIGHEST_PROTOCOL)
	with open(output_files['tcc_csr'], 'wb') as outf:
		pickle.dump(tcc_matrix_csr, outf, pickle.HIGHEST_PROTOCOL)
	
		
	print('\tNormalizing TCCs')
	tcc_normalized = normalize(tcc_matrix_csr, norm='l1', axis=0) 
	tcc_normalized_transposed = tcc_normalized.transpose()
	
	print('\tComputing L1 distance')
	l1_distance = pairwise_distances(
		tcc_normalized_transposed,
		metric=get_l1_distance,
		n_jobs=args['threads'])
	
	print('\tSaving outputs as .npy matrices')
	np.save(output_files['l1_dist'], l1_distance)
	np.save(output_files['nonzero_eq'], nonzero_equiv_classes)
	np.save(output_files['tcc_norm_t'], tcc_normalized_transposed)
	return output_files	

def get_l1_distance(p,q):
    return cityblock(p,q).sum()

def check_pipeline_input(args, kallisto):
	
	assert not (args['10xgenomics'] and args['dropseq']), \
		'10xgenomics and dropseq options are mutually exclusive'
	assert Path(kallisto).is_file() or shutil.which(kallisto), \
		'Cannot find kallisto executable %s' % kallisto
	assert os.path.exists(args['reads']), \
		'Cannot find reads file %s' % args['reads']	
	assert os.path.exists(args['barcodes']), \
		'Cannot find barcodes file %s' % args['barcodes']
	
	if args['kallisto_idx'] is not None:
		assert(os.path.exists(args['kallisto_idx'])), \
			'Cannot find kallisto index %s' % \
			args['kallisto_idx']
	if args['10xgenomics']:
		assert(os.path.exists(args['umis'])), \
			'Cannot find reads file %s' % args['umis']

def check_split_input(args):
	"""
	More here
	"""
	assert args['barcode_start'] < args['barcode_end'], \
		'Barcode end position %i is <= start position %i' % \
		(args['barcode_end'], args['barcode_start'])	
	assert args['umi_start'] <= args['umi_end'], \
		'UMI end position %i is less than start position %i' % \
		(args['umi_end'], args['umi_start'])	
	assert os.path.exists(args['reads']), \
		'Cannot find reads file %s' % args['reads']	
	assert os.path.exists(args['barcodes']), \
		'Cannot find barcodes file %s' % args['barcodes']
	assert args['kmer_size'] > 0, \
		'Kmer size must be positve. %i' % args['kmer_size']

def get_args(args=None):
	if args is None:
		args = sys.argv[1:]

	parser = argparse.ArgumentParser(
		description = 'This script splits reads for dropseq data',
		formatter_class = argparse.ArgumentDefaultsHelpFormatter)
		
	parser.add_argument('--dropseq',
		help='Use barcode / umi positions from Macosko et al Dropseq dataset',
		required=False,
		action='store_true')
	parser.add_argument('--10xgenomics',
		help='Use barcode / umi positions from 10xGenomics v2 chemistry data',
		required=False,
		action='store_true')
	
	parser.add_argument('--barcodes', 
		type=str, 
		help='Cell barcodes file name or comma separarted list. File I1 for 10x genomics', 
		required=True)
	parser.add_argument('--umis', 
		type=str,
		help='Umis (R1) file name or comma separarted list. Only required for 10x genomics', 
		required=False)
	parser.add_argument('--reads', 
		type=str, 
		help='RNAseq reads file name or comma separarted list. File R2 for 10x genomics', 
		required=True)
	parser.add_argument('--output_dir', 
		type=str, 
		help='Directory where outputs are written', 
		required=True)
	
	parser.add_argument('--kmer_size', 
		type=int, 
		help='Size of kmers for making barcode De Bruijn graph.', 
		default=None)
	parser.add_argument('--breadth', 
		type=int, 
		help='How many nodes search.', 
		default=10000)
	parser.add_argument('--depth', 
		type=int, 
		help='How many paths to search for at each node.', 
		default=10)
	parser.add_argument('--threads', 
		type=int, 
		help='Number of threads to use.', 
		default=32)
	parser.add_argument('--kallisto_idx',
		type=str,
		help='Path to kallisto index. Optional.' + \
			'If not provided, script terminates after splitting barcodes',
		default=None)
	parser.add_argument('--min_dist',
		type=int,
		help='Minimum Hamming distance between error-corrected barcodes.',
		default=1)
	parser.add_argument('--num_cells',
		type=int,
		help='Estimated number of cells.',
		default=None)
	
	parser.add_argument(
		'--barcode_start',
		type=int,
		default=None)
	parser.add_argument(
		'--barcode_end',
		type=int,
		default=None)
	parser.add_argument(
		'--umi_start',
		type=int,
		default=0)
	parser.add_argument(
		'--umi_end',
		type=int,
		default=0)

	return vars(parser.parse_args(args))

