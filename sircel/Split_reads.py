"""
Akshay Tambe
Pachter and Doudna groups


Split reads for dropseq data
1. Index kmers
		Produce a dict kmer_index 
		kmer_index[kmer] -> list of read line numbers that contain this kmer
2. Find cyclic paths
		pick a popular kmer
		get all reads that contain the kmer
		make subgraph from that subset of reads
		get best path(s) starting at the starting kmer
3. Threshold paths 
		Histogram of path weights has a minimum
4. Assign reads
		For each read, find the path that it shares the most kmers with
"""

import os
import sys
import time
import json
import gc
import numpy as np

from collections import Counter, namedtuple
from itertools import repeat, chain
from multiprocessing import Pool
from Levenshtein import distance, hamming
from scipy import signal
 
from sircel.utils import IO_utils, Plot_utils, Logger
from sircel.utils.Graph_utils import Edge, Graph, Path

np.random.seed(0)

args = {}
output_files = {}
output_dir = ''

def run_all(cmdline_args):
	print('Splitting reads by barcodes')
	
	global args
	global output_files
	global output_dir
	
	args = cmdline_args
	output_dir = args['output_dir']
	output_files['log'] = '%s/run_log.txt' % output_dir
	
	Logger.start(output_files['log'])
	start_time = time.time()
	
	reads_unzipped = args['reads']
	barcodes_unzipped = args['barcodes']
	print('Building kmer index')
	kmer_index, kmer_counts, subsamp_pearson = get_kmer_index(barcodes_unzipped)
	output_files['subsamp_pearson_plot'] = subsamp_pearson
	print('\t%i unique kmers indexed' % len(kmer_counts.items()))
	
	print('Finding cyclic paths in the barcode de Briujn graph')
	cyclic_paths = find_paths(
		(kmer_index,
		kmer_counts,
		barcodes_unzipped, 
		reads_unzipped,
		output_dir))	
	print('\t%i cyclic paths found' % len(cyclic_paths))
	output_files['all_paths'] = IO_utils.save_paths_text(
		output_dir, cyclic_paths, prefix='all')
	
	print('Thresholding paths')
	(top_paths, fit_out) = threshold_paths(
		output_dir, cyclic_paths, args['num_cells'])
	output_files.update(fit_out)
	consensus_bcs = set([tup[0] for tup in top_paths])
	
	print('Assigning reads by kmer compatability')
	reads_assigned_db, reads_assigned_pipe = assign_all_reads(
		(consensus_bcs,
		reads_unzipped, 
		barcodes_unzipped))
	
	print('Splitting reads by cell')
	output_files['split'], reads_per_cell = write_split_fastqs(
		(consensus_bcs,
		reads_assigned_db, 
		reads_assigned_pipe,
		output_dir,
		reads_unzipped,
		barcodes_unzipped))
	
	#update paths list
	top_paths = update_paths_list(top_paths, reads_per_cell)
	output_files['thresholded_paths'] = IO_utils.save_paths_text(
		output_dir, top_paths, prefix='threshold')
	
	current_time = time.time()
	elapsed_time = current_time - start_time
	
	Logger.stop()
	return(output_files, elapsed_time)
	
def get_kmer_index(barcodes_unzipped):
	"""
	Args:
		barcodes_unzipped (str): filename for unzipped barcodes fq
	
	Returns
		kmer_idx (dict): map of kmer to list of line offsets for reads 
			that contain that kmer
		kmer_counts (dict): map of kmer to absolute counts
	
	This method returns a kmer index and counts dict for a random
	subset of the dataset. The size of the subset attempts to be the
	minimal number of reads whose kmer spectrum is representative
	of the data
	
	General approach:
		initialize:
			get a random chunk of reads based on line offsets
			compute kmer counts
		loop:
			get a new chunk of reads and combine with prevoius chunks
			compute kmer counts for the new chunk
			compare kmer counts with previous iteration
		terminate when:
			pearsonR >= some cutoff value
	
	"""
	PEARSONR_CUTOFF = 0.999
	MIN_ITERS = 10
	BUFFER_SIZE = 10000
	
	length = args['barcode_end'] - args['barcode_start']
	pool = Pool(processes = args['threads'])
	
	read_count = 0
	kmer_idx = {}
	counts_corr_coefs = []
	num_reads = []	
	
	bc_file = open(barcodes_unzipped, 'rb')
	read_chunks_iter = IO_utils.get_read_chunks(
			bc_file,
			random = True,
			BUFFER_SIZE = BUFFER_SIZE)
	chunk_num = 0
	while True:
		try:
			reads_chunk = next(read_chunks_iter)
			chunk_num += 1
		except StopIteration:
			break
								
		read_count += len(reads_chunk)
		num_reads.append(read_count)
		chunk_kmer_indices = pool.map(
			index_read,
			reads_chunk)
			#chunk_kmer_indices is a list of dicts
		old_kmer_counts = get_kmer_counts(kmer_idx)
			#kmer counts before updating with chunk_kmer_indexes
		
		for element in chunk_kmer_indices:
			for(key, read_offsets) in element.items():
				#read_offsets: [offset1, offset2, offset3 ...]
				if key not in kmer_idx:
					kmer_idx[key] = []
				kmer_idx[key] = kmer_idx[key] + read_offsets
		
		del(chunk_kmer_indices)
		_ = gc.collect()
		
		new_kmer_counts = get_kmer_counts(kmer_idx)
		#check kmer count correlation
		counts_corr_coef = get_kmer_count_correlation(
			old_kmer_counts, new_kmer_counts)
		counts_corr_coefs.append(counts_corr_coef)
		print('\t%i reads indexed. Running pearsonr is %f' % \
			(read_count, counts_corr_coef))
		
		if(len(counts_corr_coefs) >= MIN_ITERS) and \
			(counts_corr_coef > PEARSONR_CUTOFF):
			break
		
	bc_file.close()
	pool.close()
	
	return (kmer_idx, 
		new_kmer_counts,
		Plot_utils.plot_kmer_subsamp_pearson(
			output_dir,
			counts_corr_coefs,
			num_reads))

def index_read(params):
	"""
	Args
		params (tuple):
			barcodes_data (str): sequence of read_1 (barcode)
			barcodes_offset (int): line offset for this read
	Returns
		kmer_index (dict): 
	"""
	(barcodes_data, barcodes_offset) = params
	
	kmer_index = {}
	read_kmers = IO_utils.get_cyclic_kmers(
		barcodes_data, 
		args['kmer_size'],
		args['barcode_start'], 
		args['barcode_end'])
	for(kmer, _) in read_kmers:
		if(kmer not in kmer_index.keys()):
			kmer_index[kmer] = []
		kmer_index[kmer].append(barcodes_offset)
	return kmer_index

def get_kmer_counts(kmer_idx):
	"""
	
	
	"""
	kmer_counts = {}
	for kmer, offsets in kmer_idx.items():
		kmer_counts[kmer] = len(offsets)
	return kmer_counts

def get_kmer_count_correlation(kmer_counts_a, kmer_counts_b):
	from scipy.stats.stats import pearsonr
	
	common_elements = kmer_counts_a.keys() & kmer_counts_b.keys()
	if(len(common_elements) <= 1):
		return 0
	
	x = []
	y = []
	for element in common_elements:
		x.append(kmer_counts_a[element])
		y.append(kmer_counts_b[element])
	
	corr_coef, pval = pearsonr(x, y)
	return corr_coef
	
def find_paths(params, starting_kmers = None):
	(	kmer_index,
		kmer_counts,
		barcodes_unzipped, 
		reads_unzipped,
		output_dir) = params
	barcode_length = args['barcode_end'] - args['barcode_start']
	kmers_sorted = [tup[0] for tup in sorted(
		list(kmer_counts.items()),
		key = lambda tup: tup[1],
		reverse = True)]
	
	if(starting_kmers == None):
		starting_kmers = []
		for kmer in kmers_sorted:
			if(kmer[0] == '$'):
				starting_kmers.append((kmer, kmer_index[kmer]))
			if(len(starting_kmers) >= args['breadth']):
				break	
	else:
		starting_kmers_tmp = []
		for kmer in starting_kmers:
			starting_kmers_tmp.append(kmer, kmer_index[kmer])
		starting_kmers = starting_kmers_tmp
		
	
	pool = Pool(processes = args['threads'])
	paths = []
	for group in IO_utils.grouper(
		starting_kmers, args['threads']):
		
		kmers_group = [tup[0] for tup in group]
		offsets_group = [tup[1] for tup in group]
		paths_group = pool.map(find_path_from_kmer, zip(
			kmers_group,
			offsets_group,
			repeat(barcodes_unzipped),
			repeat(barcode_length)))
		paths += [item for sublist in paths_group for item in sublist]
	pool.close()
	return paths

def find_path_from_kmer(params):
	(	starting_kmer,
		offsets, 
		barcodes_unzipped,
		barcode_length) = params
	#1. build subgraph
	subgraph = build_subgraph(offsets, barcodes_unzipped)
	#2. find paths
	node = starting_kmer[0:-1]
	neighbor = starting_kmer[1:]
	paths = []
	paths_iter = subgraph.find_all_cyclic_paths(
			node, neighbor, barcode_length + 1)
	counter = 1
	
	while(True):
		try:
			path = next(paths_iter)
		except StopIteration:
			break
		if(not path.is_cycle()):
			break
		seq = path.get_sequence_circular()
		weight = path.get_cycle_weight()
		nodes = [edge.get_sequence() for edge in path.edges]
		paths.append((seq, weight, counter))
		if(counter > args['depth']):
			break
		counter += 1
	return merge_paths(paths)

def build_subgraph(reads_in_subgraph, barcodes_unzipped):
	bc_file = open(barcodes_unzipped, 'rb')
	barcodes_iter = IO_utils.read_fastq_random(
		bc_file, offsets = reads_in_subgraph)
	subgraph_kmer_counts = Counter()
	while(True):
		try:
			barcode_data, _ = next(barcodes_iter)
		except StopIteration:
			break	
		read_kmers = IO_utils.get_cyclic_kmers(
			barcode_data, 
			int(args['kmer_size']),
			int(args['barcode_start']), 
			int(args['barcode_end']))		
		for (kmer, _ ) in read_kmers:
			subgraph_kmer_counts[kmer] += 1
	bc_file.close()
	
	edges = []
	for(kmer, count) in subgraph_kmer_counts.items():
		edge = Edge(kmer[0:-1], kmer[1:], count)
		edges.append(edge)
	subgraph = Graph(edges)
	return subgraph

def threshold_paths(output_dir, paths, num_cells):
	LOCAL_WINDOW_LEN = 50
	MIN_CAPACITY = 0
	MIN_WEIGHT = 10
	
	threshold_out = {
		'slopes' : '%s/slopes.txt' % output_dir,
	}
	unique_paths = {}
	for tup in paths:
		key = tup[0]
		if(tup[1] > MIN_CAPACITY):
			if(key not in unique_paths):
				unique_paths[key] = tup
			else:
				#note- because of the subgraph method, 
				#	the same path in two different subgraphs might have different weight
				#keep the instance of the path with the higher capacity
				old_capacity = unique_paths[key][1]
				current_capacity = tup[1]
				if(current_capacity > old_capacity):
					unique_paths[key] = tup
					#keep only unique paths with a capacity higher than some threshold value
	unique_paths_sorted = sorted(
		unique_paths.values(),
		key = lambda tup: tup[1],
		reverse = True)
	
	path_weights = [tup[1] for tup in unique_paths_sorted if tup[1] >= MIN_WEIGHT]
	for i in range(2 * LOCAL_WINDOW_LEN):
		path_weights.append(MIN_WEIGHT)
	
	grad = [-1 * i for i in \
				local_lin_fit(np.log10(path_weights),
				window_len=LOCAL_WINDOW_LEN)]   
	second_grad = local_lin_fit(grad, window_len = LOCAL_WINDOW_LEN)
	
	lmax = get_lmax(second_grad, LOCAL_WINDOW_LEN)
	threshold = get_threshold((
		grad, 
		second_grad,
		lmax, 
		num_cells, 
		unique_paths_sorted,
		LOCAL_WINDOW_LEN))
	top_paths = unique_paths_sorted[0:threshold]

	print('\t%i paths remain after thresholding' % len(top_paths))
	threshold_out['paths_threshold_plot'] = Plot_utils.plot_path_threshold(
		(output_dir, 
			path_weights,
			grad,
			second_grad,
			lmax,
			threshold,
			LOCAL_WINDOW_LEN))
	
	return top_paths, threshold_out

def get_lmax(second_grad, LOCAL_WINDOW_LEN):
	#finds zeros in
	lmax = []
	for i in range(len(second_grad) - 1):
		if(second_grad[i] > 0 and second_grad[i + 1] <= 0):
			lmax.append(int(i + LOCAL_WINDOW_LEN))
	return lmax

def get_threshold(params):
	(grad, 
		second_grad, 
		lmax, 
		num_cells, 
		unique_paths_sorted, 
		LOCAL_WINDOW_LEN) = params
	
	#if there is no lmax, return the number of paths
	if len(lmax) == 0:
		return len(unique_paths_sorted)
	
	#if there is an expected number of cells, 
	#filter local maxima to those nearest to the lmax
	if num_cells != None:
		MAX_DISTANCE = 150
		lmax_thresholded = []
		for i in lmax:
			if np.fabs(i - num_cells)  <= MAX_DISTANCE:
				lmax_thresholded.append(i)
		lmax = lmax_thresholded
		
	change_coords = lambda i : int(i - LOCAL_WINDOW_LEN / 2)
	#return the local max with highest value (steepest inflection)
	try:
		threshold = lmax[-1]
	except IndexError:
		return len(unique_paths_sorted)
	for i in lmax:		
		if(grad[change_coords(i)] > grad[change_coords(threshold)]):
			threshold = i
	return min(threshold, len(unique_paths_sorted))

def local_lin_fit(y, window_len=10):
	from scipy.optimize import curve_fit
	
	num_windows = len(y) - window_len
	slopes = []
	for window_start in range(0, num_windows):
		window_x = range(window_start, window_start + window_len)
		window_y = y[window_start : window_start + window_len]
		coeff, var_matrix = curve_fit(
			linear,
			window_x,
			window_y,
			p0 = [window_y[-1] - window_y[0], window_y[0]])
		(slope, intercept) = coeff
		slopes.append(slope)
	return slopes

def linear(x, *p):
	(slope, intercept) = p
	return slope * x + intercept
	
def merge_paths(paths, MIN_DIST = 1):	
	paths_sorted = sorted(paths, key = lambda tup: tup[1])
	num_paths = len(paths)
	
	paths_merged = {tup[0] : tup for tup in paths_sorted}
	get_seq = lambda tup: tup[0]
	for (i, path) in enumerate(paths_sorted):
		for j in range(i+1, num_paths):
			ham_dist = hamming(get_seq(paths[i]), get_seq(paths[j]))
			if(ham_dist <= MIN_DIST):
				bad_path = min([paths[i], paths[j]], key = lambda tup: tup[1])
				if(get_seq(bad_path) in paths_merged.keys()):
					del(paths_merged[get_seq(bad_path)])
	return list(paths_merged.values())

def assign_all_reads(params):
	(	consensus_bcs,
		reads_unzipped, 
		barcodes_unzipped) = params
	
	BUFFER_SIZE = 100000
	MAX_KMER_SIZE = args['barcode_end'] - args['barcode_start']
	MIN_KMER_SIZE = 6
	
	reads_assigned_db, reads_assigned_pipe = IO_utils.initialize_redis_pipeline()
	pool = Pool(processes = args['threads'])
	
	#print('\tMapping kmers to consensus barcodes')
	#kmer_map = map_kmers_to_bcs(consensus_bcs, MIN_KMER_SIZE, MAX_KMER_SIZE)
	
	print('\tAssigning reads to consensus barcodes')
	read_count = 0
	num_unassigned = 0
	reads_f = open(reads_unzipped, 'rb')
	barcodes_f = open(barcodes_unzipped, 'rb')
		
	for reads_chunk, barcodes_chunk in zip(
		IO_utils.get_read_chunks(
			reads_f,
			random = False,
			BUFFER_SIZE = BUFFER_SIZE),
		IO_utils.get_read_chunks(
			barcodes_f,
			random = False,
			BUFFER_SIZE = BUFFER_SIZE)):
		read_count += len(reads_chunk)
		
		#if args['split_levenshtein']:
		assignments = pool.map(assign_read_levenshtein,
			zip(
				repeat(consensus_bcs),
				reads_chunk,
				barcodes_chunk))		
		"""
		else:
			#this is a pipeline for reviwer expts only
			#works quite poorly, see simulation results
			assignments = pool.map(assign_read_kmers, 
				zip(
				repeat(kmer_map),
				repeat(MIN_KMER_SIZE),
				repeat(MAX_KMER_SIZE),
				reads_chunk,
				barcodes_chunk))
		"""
			
		for (assignment, offset1, offset2) in assignments:
			if(assignment == 'unassigned'):
				num_unassigned += 1
			#reads_assigned[assignment].append((offset1, offset2))
			reads_assigned_pipe.rpush(
			 	assignment.encode('utf-8'), 
			 	offset1, offset2)
				
		reads_assigned_pipe.execute()
		print('\tProcessed %i reads' % read_count)
	
	reads_f.close()
	barcodes_f.close()
	pool.close()
	
	print('\t%i reads could not be assigned' % num_unassigned)
	#return pickle_files
	return reads_assigned_db, reads_assigned_pipe

def initialize_reads_assigned(consensus_bcs):
	reads_assigned = {}
		#key / value map of: [cell name] :-> list of line offsets
	for bc in consensus_bcs:
		reads_assigned[bc] = []
	reads_assigned['unassigned'] = []
	return reads_assigned

def map_kmers_to_bcs(consensus_bcs, MIN_KMER_SIZE, MAX_KMER_SIZE):
	
	kmer_map = {}
	for kmer_size in range(MAX_KMER_SIZE, MIN_KMER_SIZE, -1):
		kmer_map_ = \
			map_kmers_to_bcs_fixed_k(consensus_bcs, kmer_size)
		kmer_map = dict(list(kmer_map_.items()) + list(kmer_map.items()))		
	return kmer_map

def map_kmers_to_bcs_fixed_k(consensus_bcs, kmer_size):
	kmers_to_paths = {}
	for cell_barcode in consensus_bcs:
		kmers = IO_utils.get_cyclic_kmers(
			['na', cell_barcode, 'na', cell_barcode],
			kmer_size,
			0,
			len(cell_barcode),
			indel=True)
		for (kmer, _) in kmers:
			if(kmer not in kmers_to_paths.keys()):
				kmers_to_paths[kmer] = []
			kmers_to_paths[kmer].append(cell_barcode)
	return kmers_to_paths
	
def assign_read_kmers(params):
	"""
	Assigns a single read to a cell barcode by kmer compatibility
	args (tuple)
		kmers_to_paths: dict of kmer -> list of paths that contain it
		min_kmer_size
		max_kmer_size
		read: list of fastq entry lines
	"""
	(kmer_map,
		min_kmer_size,
		max_kmer_size,
		(reads_data, reads_offset),
		(barcodes_data, barcodes_offset)) = params
	
	for kmer_size in range(max_kmer_size, min_kmer_size, -1):
		read_kmers = IO_utils.get_cyclic_kmers(
			barcodes_data, 
			kmer_size,
			args['barcode_start'], 
			args['barcode_end'],
			indel = True)
		bcs, is_assigned, is_unique = get_most_common_bc(
			kmer_map, read_kmers)
		if is_assigned and is_unique:
			return (bcs[0], reads_offset, barcodes_offset)
		#outherwise decrement kmer size and try again
	return ('unassigned', reads_offset, barcodes_offset)

def get_most_common_bc(kmer_map, read_kmers):
	compatable_bcs = {}
	for (kmer, _) in read_kmers:
		bcs = kmer_map.get(kmer, None)
		if bcs != None:
			increment = 1.0 / len(bcs)
			for bc in bcs:
				if bc not in compatable_bcs:
					compatable_bcs[bc] = 0
				compatable_bcs[bc] += increment
	most_common = None
	highest_count = 0
	for bc, count in compatable_bcs.items():
		if count > highest_count:
			highest_count = count
			most_common = [bc]
		elif count == highest_count and count > 0:
			most_common.append(bc)
	
	if most_common == None:
		return None, False, False
	elif len(most_common) == 1:
		return most_common, True, True
	else:
		return most_common, True, False

def assign_read_levenshtein(params):	
	(consensus_bcs,
		(reads_data, reads_offset),
		(barcodes_data, barcodes_offset)) = params
	
	obs_bc = reads_data[1].strip()[ \
		args['barcode_start']: args['barcode_end']]
	#first check for perfect match	
	if obs_bc in consensus_bcs:
		return (obs_bc, reads_offset, barcodes_offset)
	
	#otherwise minimize levenshtein distance
	min_lev_dist = None
	assignment = []
	for consensus_bc in consensus_bcs:
		lev_dist = distance(obs_bc, consensus_bc)
		if min_lev_dist == None or lev_dist < min_lev_dist:
			min_lev_dist = lev_dist
			assignment = [consensus_bc]
		#in the case of a tie,
		elif lev_dist == min_lev_dist:
			assignment.append(consensus_bc)
	#return the best unique assignment
	if len(assignment) == 1:
		return (assignment[0], reads_offset, barcodes_offset)
	#or don't assign read (in the case of a tie)
	return ('unassigned', reads_offset, barcodes_offset)

def update_paths_list(top_paths, reads_per_cell):
	updated_paths = []
	for (seq, capacity, depth) in top_paths:
		num_reads = reads_per_cell.get(seq, [None])
		updated_paths.append((seq, capacity, depth, num_reads))
	return updated_paths

def write_split_fastqs(params):
	import gzip
	(	consensus_bcs,
		reads_assigned_db,
		reads_assigned_pipe,
		output_dir,
		reads_unzipped,
		barcodes_unzipped) = params
	
	split_dir = '%s/reads_split' % output_dir
	if not os.path.exists(split_dir):
		os.makedirs(split_dir)
	output_files = {'batch' : '%s/batch.txt' % (split_dir)}
	batch_file = open(output_files['batch'], 'w')
	
	reads_per_cell = {}
	consensus_bcs.add('unassigned')
	
	for cell in consensus_bcs:
		
		try:
			cell_offsets = IO_utils.get_from_db(reads_assigned_pipe, [cell])
		except IndexError:
			pass
		
		#cell_offsets = IO_utils.read_from_pickle(reads_assigned_pickled, cell)
		cell_name = 'cell_%s' % cell
		
		#initialie all readers and writers
		output_files[cell_name] = {
			'reads' : '%s/%s_reads.fastq.gz' % (split_dir, cell_name),
			'barcodes' : '%s/%s_barcodes.fastq.gz' % (split_dir, cell_name),
			'umi' : '%s/%s.umi.txt' % (split_dir, cell_name)}
		batch_file.write('%s\t%s\t%s\n' % \
			(cell_name, 
			output_files[cell_name]['umi'], 
			output_files[cell_name]['reads']))
		reads_writer = gzip.open(output_files[cell_name]['reads'], 'wb')
		barcodes_writer = gzip.open(output_files[cell_name]['barcodes'], 'wb')
		umi_writer = open(output_files[cell_name]['umi'], 'wb')
		reads_f = open(reads_unzipped, 'rb')
		barcodes_f = open(barcodes_unzipped, 'rb')
		
		reads_iter = IO_utils.read_fastq_random(
			reads_f, 
			offsets = 
				[cell_offsets[i] for i in range(len(cell_offsets)) if i % 2 == 0])
		barcodes_iter = IO_utils.read_fastq_random(
			barcodes_f,
			offsets = 
				[cell_offsets[i] for i in range(len(cell_offsets)) if i % 2 == 1])
		reads_in_cell = 0
		while(True):
			try:
				reads_data, _ = next(reads_iter)
				barcodes_data, _ = next(barcodes_iter)
				reads_in_cell += 1
			except StopIteration:
				break
			
			reads_data[0] += ' %s' % cell_name.replace('_', ':')
			reads_data[0] = reads_data[0].replace(' ', '_')
			barcodes_data[0] += ' %s' % cell_name.replace('_', ':')	
			barcodes_data[0] = barcodes_data[0].replace(' ', '_')
					
			umi = barcodes_data[1][
				int(args['umi_start']): int(args['umi_end'])]
			reads_writer.write(
				('\n'.join(reads_data) + '\n').encode('utf-8'))
			barcodes_writer.write(
				('\n'.join(barcodes_data) + '\n').encode('utf-8'))
			umi_writer.write((umi + '\n').encode('utf-8'))
		
		reads_writer.close()
		umi_writer.close()
		barcodes_writer.close()
		reads_f.close()
		barcodes_f.close()
		
		print('\tWrote %i reads to file:\t%s' % \
			(reads_in_cell, cell_name))
		reads_per_cell[cell] = reads_in_cell
	batch_file.close()
	return output_files, reads_per_cell

def get_args():
	import argparse
	
	parser = argparse.ArgumentParser(
		description = 'This script splits reads for dropseq data')
	parser.add_argument('--barcodes', 
		type=str, 
		help='Barcodes file name (unzipped)', 
		required=True)
	parser.add_argument('--barcodes_offset',
		type=list,
		help='Fq entry line offsets for barcoddes file',
		required=True)
	parser.add_argument('--reads', 
		type=str, 
		help='RNAseq reads file name (unzipped)', 
		required=True)
	parser.add_argument('--reads_offset',
		type=list,
		help='Fq entry line offsets for reads file',
		required=True)
		
	parser.add_argument('--output_dir', 
		type=str, 
		help='Directory where outputs are written', 
		required=True)
	parser.add_argument('--barcode_start', 
		type=int, 
		help='Start position of barcode.', 
		default=0)
	parser.add_argument('--barcode_end', 
		type=int, 
		help='End position of barcode.', 
		default=12)
	parser.add_argument('--umi_start', 
		type=int, 
		help='Start position of UMI.', 
		default=12)
	parser.add_argument('--umi_end', 
		type=int, 
		help='End position of UMI.', 
		default=20)
	parser.add_argument('--kmer_size', 
		type=int, 
		help='Size of kmers for making barcode De Bruijn graph.', 
		default=7)
	parser.add_argument('--depth', 
		type=int, 
		help='Fraction of edge weight at starting node to assign to path.', 
		default=10)
	parser.add_argument('--breadth', 
		type=int, 
		help='How many nodes search.', 
		default=10000)
	parser.add_argument('--threads', 
		type=int, 
		help='Number of threads to use.', 
		default=32)
	parser.add_argument('--num_cells',
		type=int,
		help='Estimated number of cells.',
		default=None)
	
	#only for reviewer expts. never actually use this!
	parser.add_argument('--split_levenshtein',
		type = bool,
		help = argparse.SUPPRESS,
		default = False)
	
	return vars(parser.parse_known_args())

if __name__ == '__main__':
	cmdline_args = get_args()	
	output_files, elapsed_time = run_all(cmdline_args)
	print('Done. Time elapsed: %f seconds' % elapsed_time)



