"""
Akshay Tambe
Pachter and Doudna groups

Evaluate_cell_errors.py
	Read run_outputs.json, get list of *barcodes.fastq.gz
	For each barcode file
		get the consensus / error-corrected barcode seq (from fname)
		for each read
			compute levenshtein and hamming distance to consensus
			record error freq vs position as numpy array
"""

import sys
from sircel import IO_utils
from sircel import Plot_utils
from Levenshtein import distance, hamming
from multiprocessing import Pool
from collections import Counter
from itertools import repeat

def run_all():
	run_outputs = sys.argv[1]
	split_files = read_run_outputs(run_outputs, 'split')
	split_args = read_run_outputs(run_outputs, 'args')
	
	del split_files['batch']
	del split_files['cell_unassigned']
	distances = get_cell_error_rate(split_files, split_args)
		#list of tuples: (cell_name (str), lev_dist (counter), ham_dist (counter) )
	bc_len = split_args['barcode_end'] - split_args['barcode_start']
	edit_distances_plt = Plot_utils.plot_cell_distance_hmap( \
		split_args['output_dir'], distances, bc_len)

def read_run_outputs(run_outputs, field = None):
	import json
	with open(run_outputs) as r:
		if field != None:
			return json.load(r)[field]
		else:
			return json.load(r)

def get_cell_error_rate(split_files, split_args):
	pool = Pool(processes = 8)
	cell_names = list(split_files.keys())
	
	distances = pool.map(
		get_single_cell_error_rate,
		zip(cell_names, repeat(split_files), repeat(split_args)))
		#list of tuples (lev_dist, ham_dist), which are counters
		#contain distributions of levenshtein and hamming distances
	
	ret = []
	for cell_name, (lev_dist, ham_dist, dist_diff) in zip(cell_names, distances):
		ret.append((cell_name, lev_dist, ham_dist, dist_diff))
	return ret

def get_consensus_seq(cell_name):
	assert cell_name.startswith('cell_'), \
		'Cell name does not follow expected format: %s' % cell_name
	return cell_name[5:]

def get_single_cell_error_rate(params):
	(cell_name, split_files, split_args) = params
	
	consensus = get_consensus_seq(cell_name)
	bc_fq = split_files[cell_name]['barcodes']
	
	lev_dist = Counter()
	ham_dist = Counter()
	dist_diff = Counter()
	for read in read_fastq_gz(bc_fq):
		seq = read[1].decode('utf-8').strip()[ \
			split_args['barcode_start']: split_args['barcode_end']]
		lev = distance(consensus, seq)
		ham = hamming(consensus, seq)
		
		lev_dist.update(str(lev))
		ham_dist.update(str(ham))
		dist_diff.update(str(lev - ham))
	return (lev_dist, ham_dist, dist_diff)

def read_fastq_gz(fq_file):
	import gzip as gz
	with gz.open(fq_file, 'rb') as r:
		for lines in IO_utils.grouper(r, 4):
			yield lines


if __name__ == '__main__':
	run_all()	
		