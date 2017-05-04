"""
Akshay Tambe
Pachter and Doudna groups


Split reads for dropseq data
1. Index kmers
		Produce a dict kmer_index 
		kmer_index[kmer] -> list of read line numbers that contain this kmer
	Use redis database to store
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

import IO_utils
import argparse
import os
import sys
import time
import json
import collections
import itertools
import gc
import numpy as np
np.random.seed(0)


from Graph_utils import Edge, Graph, Path
from multiprocessing import Pool 


args = {}

def run_all(cmdline_args):
	global args
	args = cmdline_args
	output_files = {}
	output_dir = args['output_dir']
	output_files['log'] = '%s/run_log.txt' % output_dir
	sys.stdout = IO_utils.Logger(output_files['log'])
	
	print('Running Dropseq_subgraphs\nArgs:%s' % \
		json.dumps(args, indent = 5))
	
	start_time = time.time()
	reads_unzipped = args['reads']
	barcodes_unzipped = args['barcodes']
	print('Initializing redis-db server')
	kmer_idx_db, kmer_idx_pipe = IO_utils.initialize_redis_pipeline(db=0)
	print('Indexing reads by circularized kmers')
	kmer_idx_pipe = get_kmer_index_db(
		(kmer_idx_pipe, barcodes_unzipped, reads_unzipped))
	kmer_counts = get_kmer_counts(kmer_idx_db, kmer_idx_pipe)
	print('\t%i kmers indexed' % len(kmer_counts.items()))
	
	print('Finding cyclic paths')
	cyclic_paths = find_paths(
		(kmer_idx_pipe,
		kmer_counts,
		barcodes_unzipped, 
		reads_unzipped,
		output_dir))	
	print('\t%i cyclic paths found' % len(cyclic_paths))
	output_files['all_paths'] = IO_utils.save_paths_text(
		output_dir, cyclic_paths, prefix='all')		
	"""
	print('Merging similar paths by Hamming distance')
	merged_paths = merge_paths(cyclic_paths)
	output_files['merged_paths'] = IO_utils.save_paths_text(
		output_dir, merged_paths, prefix='merged')	
	print('\t%i paths remaining after merging' % len(merged_paths))
	"""
			
	print('Thresholding paths')
	(threshold, top_paths, fit_out) = threshold_paths(output_dir, cyclic_paths)
	output_files.update(fit_out)
	print('\tThreshold is %i' % threshold)
	print('\t%i paths have weight higher than the threshold' % len(top_paths))
	output_files['thresholded_paths'] = IO_utils.save_paths_text(
		output_dir, top_paths, prefix='threshold')
	kmer_idx_db.flushall()
	
	print('Assigning reads')
	reads_assigned_db, reads_assigned_pipe = assign_all_reads(
		top_paths, reads_unzipped, barcodes_unzipped)
	
	print('Splitting reads by cell')
	output_files['split'] = write_split_fastqs(
		(reads_assigned_db,
		reads_assigned_pipe,
		output_dir,
		reads_unzipped,
		barcodes_unzipped))
	
	print('Flushing assignments index')
	reads_assigned_db.flushall()
	
	current_time = time.time()
	elapsed_time = current_time - start_time
	return(output_files, elapsed_time)
	
def get_kmer_index_db(params):
	
	(	kmer_idx_pipe,
		barcodes_unzipped, 
		reads_unzipped) = params
	length = args['barcode_end'] - args['barcode_start']
	pool = Pool(processes = args['threads'])
	
	read_count = 0
	for reads_chunk in IO_utils.get_read_chunks(
		barcodes_unzipped,
		reads_unzipped,
		lines = None,
		random_subset = args['index_depth']):
				
		read_count += len(reads_chunk)
		chunk_kmer_indices = pool.map(
			index_read,
			reads_chunk)
		#chunk_kmer_indices is a list of dicts
		for element in chunk_kmer_indices:
			for(key, offsets) in element.items():
				#value is a list: [offset1, offset2, offset3 ...]
				value = ','.join(['%i' % i for i in offsets]) + ','
				kmer_idx_pipe.append(key.encode('utf-8'), value.encode('utf-8'))
		pipe_out = kmer_idx_pipe.execute()#slow. chunks are large...
		del(chunk_kmer_indices); _ = gc.collect()
		print('\t%i reads indexed' % read_count)
		
		"""
		#snapshot the database 
		try:
			pipe_out = kmer_idx_pipe.save()
		except redis.exceptions.ResponseError:
			pass
		"""
	return kmer_idx_pipe

def index_read(params):
	(	reads_data,
		reads_offset,
		barcodes_data, 
		barcodes_offset) = params
		
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

def get_kmer_counts(kmer_idx_db, kmer_idx_pipe):
	kmer_counts = {}
	for key in kmer_idx_db.keys():
		entries = IO_utils.get_from_db(kmer_idx_pipe,[key])[0]
		kmer_counts[key.decode('utf-8')] = len(entries)
	return kmer_counts

def find_paths(params):
	(	kmer_idx_pipe,
		kmer_counts,
		barcodes_unzipped, 
		reads_unzipped,
		output_dir) = params
	barcode_length = args['barcode_end'] - args['barcode_start']
	kmers_sorted = [tup[0] for tup in sorted(
		list(kmer_counts.items()),
		key = lambda tup: tup[1],
		reverse = True)]
	
	starting_kmers = []
	for kmer in kmers_sorted:
		#if(kmer[0] == '$'):
		starting_kmers.append(kmer)
		if(len(starting_kmers) >= args['breadth']):
			break	
	
	pool = Pool(processes = args['threads'])
	paths = []
	for kmers_group in IO_utils.grouper(starting_kmers, args['threads']):
		offsets_group = IO_utils.get_from_db(kmer_idx_pipe, kmers_group)
		paths_group = pool.map(find_path_from_kmer, zip(
				kmers_group,
				offsets_group,
				itertools.repeat(barcodes_unzipped),
				itertools.repeat(barcode_length)))
		paths += [item for sublist in paths_group for item in sublist]
	
	
	
	plot_cycles_multi(
		get_cycles_multi(paths), output_dir)
	
	#keep only unique paths
	unique_paths = {}
	for lst in paths:
		key = lst[0]
		if(key not in unique_paths):
			unique_paths[key] = lst
	return list(unique_paths.values())

def get_cycles_multi(paths):
	paths_multi = {}
	for path in paths:
		(seq, weight, depth, nodes) = path
		start_node = nodes[0]
		if(seq not in paths_multi.keys()):
			paths_multi[seq] = []
		paths_multi[seq].append((start_node, weight))
	return paths_multi

def plot_cycles_multi(paths_multi, output_dir):
	#plots a scatter plot of mean vs variance over mean for capacity
	
	import matplotlib as mpl
	mpl.use('Agg')
	from matplotlib import pyplot as plt
	from scipy import signal
	
	fig, ax = plt.subplots(
		nrows = 1, 
		ncols = 1,
		figsize = (4,4))
	
	mean_capacity = []
	std_capacity = []
	for lst in paths_multi.values():
		mean_capacity.append(np.mean([tup[1] for tup in lst]))
		std_capacity.append(np.std([tup[1] for tup in lst]))
	
	ax.scatter(mean_capacity, std_capacity**2)
	
	fig.savefig('%s/mean_variance_paths.pdf' % output_dir)


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
		paths.append( (seq, weight, counter, nodes) )
		if(counter > args['depth']):
			break
		counter += 1
	#merge similar paths by hamming distance
	merged_paths = merge_paths(paths)
	return merged_paths

def build_subgraph(reads_in_subgraph, barcodes_unzipped):
	barcodes_iter = IO_utils.read_fastq_random(
		barcodes_unzipped, reads_in_subgraph)
	subgraph_kmer_counts = collections.Counter()
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
	edges = []
	for(kmer, count) in subgraph_kmer_counts.items():
		edge = Edge(kmer[0:-1], kmer[1:], count)
		edges.append(edge)
	subgraph = Graph(edges)
	return subgraph

def merge_paths(paths):
	paths_sorted = sorted(paths, key = lambda tup: tup[1])
	num_paths = len(paths)
	
	get_seq = lambda paths, i: paths[i][0]
	paths_merged = {tup[0] : tup for tup in paths_sorted}
	
	for (i, path) in enumerate(paths_sorted):
		for j in range(i+1, num_paths):
			hamming = hamming_distance(get_seq(paths, i), get_seq(paths, j))
			if(hamming <= args['min_dist']):
				bad_path = min([paths[i], paths[j]], key = lambda tup: tup[1])
				if(bad_path[0] in paths_merged.keys()):
					del(paths_merged[bad_path[0]])
	return list(paths_merged.values())

def hamming_distance(seq1, seq2):
	hamming = 0
	for (i,j) in zip(seq1, seq2):
		if(i != j):
			hamming += 1
	return hamming

def threshold_paths(output_dir, paths):
	WINDOW = [200, 1000]
	LOCAL_WINDOW_LEN = 25

	import matplotlib as mpl
	mpl.use('Agg')
	from matplotlib import pyplot as plt
	from scipy import signal

	threshold_out = {
		'slopes' : '%s/slopes.txt' % output_dir,
		'paths_plot' : '%s/paths_plotted.pdf' % output_dir}

	fig, ax = plt.subplots(
		nrows = 1, 
		ncols = 2,
		figsize = (8,4))
		
	x = range(0, len(paths))
	y = sorted( [tup[1] for tup in paths], reverse=True)
	ax[0].step(x,y, label='Cum dist')

	slopes, x = local_lin_fit(
		np.log10(y), window_len=LOCAL_WINDOW_LEN)
	ax[1].scatter(
		x, slopes, color='r', alpha=0.2, s=2, label='Local gradient')
	savgol = signal.savgol_filter(slopes, 251, 4)
	ax[1].step(x, savgol, label='Savgol filter')

	threshold = int(np.argmax(
		savgol[WINDOW[0]:WINDOW[1]]) + WINDOW[0] + x[0])
	ax[1].axvline(threshold, color='k', ls='--', label='Threshold')
	ax[1].legend(loc=1)
	ax[0].axvline(threshold, color='k', ls='--', label='Threshold')
	ax[0].legend(loc=1)
	ax[0].set_yscale('log')
	
	paths_sorted = sorted(
		paths, key = lambda tup: tup[1], reverse = True)
	top_paths = paths_sorted[0:threshold]
	fig.savefig(threshold_out['paths_plot'])
	
	return threshold, top_paths, threshold_out
	
def local_lin_fit(y, window_len=10):
	from scipy.optimize import curve_fit
	num_windows = len(y) - window_len
	slopes = []
	x = []
	for window_start in range(0, num_windows):
		window_x = range(window_start, window_start + window_len)
		window_y = y[window_start : window_start + window_len]
		coeff, var_matrix = curve_fit(
			linear,
			window_x,
			window_y,
			p0=[window_y[-1] - window_y[0], window_y[0]])
		(slope, intercept) = coeff
		slopes.append(-slope)
		x.append(window_start + window_len / 2)
	return slopes, x

def linear(x, *p):
	(slope, intercept) = p
	return slope*x + intercept
		
def assign_all_reads(top_paths, reads_unzipped, barcodes_unzipped):
	MIN_KMER_SIZE = 4
	MAX_KMER_SIZE = args['barcode_end'] - args['barcode_start']
	
	#initialize vars
	reads_assigned_db, reads_assigned_pipe = IO_utils.initialize_redis_pipeline(db=1)
	kmers_to_paths = {}
	
	print('\tGetting kmers in paths')
	for path in top_paths:
		cell_barcode = path[0]
		for kmer_size in range(MIN_KMER_SIZE, MAX_KMER_SIZE):
			kmers = IO_utils.get_cyclic_kmers(
				['na', cell_barcode, 'na', cell_barcode],
				kmer_size,
				0,
				len(cell_barcode),
				indel=False)
			for (kmer, _) in kmers:
				if(kmer not in kmers_to_paths.keys()):
					kmers_to_paths[kmer] = []
				kmers_to_paths[kmer].append(cell_barcode)

	print('\tAssigning reads to paths')
	pool = Pool(processes = args['threads'])	
	read_count = 0
	num_unassigned = 0
	for reads_chunk in IO_utils.get_read_chunks(
		barcodes_unzipped,
		reads_unzipped,
		lines=None):
		
		read_count += len(reads_chunk)
		assignments = pool.map(assign_read, 
			zip(itertools.repeat(kmers_to_paths),
			itertools.repeat(MIN_KMER_SIZE),
			itertools.repeat(MAX_KMER_SIZE),
			reads_chunk))
		for (assignment, offset1, offset2) in assignments:
			if(assignment == 'unassigned'):
				num_unassigned += 1
			reads_assigned_pipe.append(
				assignment.encode('utf-8'), 
				('%i,%i,' % (offset1, offset2)).encode('utf-8'))
		reads_assigned_pipe.execute()
		print('\t%i reads assigned' % read_count)
	print('%i reads could not be assigned' % num_unassigned)
	return(reads_assigned_db, reads_assigned_pipe)
	
def assign_read(params):
	"""
	Assigns a single read to a cell barcode by kmer compatibility
	args (tuple)
		kemers_to_paths: dict of kmer -> list of paths that contain it
		min_kmer_size
		max_kmer_size
		read: list of fastq entry lines
	returns
	
	"""
	(kmers_to_paths,
		min_kmer_size,
		max_kmer_size,
		read) = params
	(reads_data,
		reads_offset,
		barcodes_data, 
		barcodes_offset) = read
	
	read_assignment = collections.Counter()
	for kmer_size in range(max_kmer_size, min_kmer_size, -1):
		read_kmers = IO_utils.get_cyclic_kmers(
			barcodes_data, 
			kmer_size,
			args['barcode_start'], 
			args['barcode_end'])
	
		for (kmer, _ ) in read_kmers:
			paths_with_kmer = kmers_to_paths.get(kmer, [])
			for path in paths_with_kmer:
				read_assignment[path] += 1
		most_common = read_assignment.most_common(1)
		if(len(most_common) == 1):
			assignment = most_common[0][0]
			return (assignment, reads_offset, barcodes_offset)
		#else decremenet kmer size
	return ('unassigned', reads_offset, barcodes_offset)

def write_split_fastqs(params):
	import gzip
	(	reads_assigned_db,
		reads_assigned_pipe,
		output_dir,
		reads_unzipped,
		barcodes_unzipped) = params
	
	split_dir = '%s/reads_split' % output_dir
	if not os.path.exists(split_dir):
		os.makedirs(split_dir)
	output_files = {'batch' : '%s/batch.txt' % (split_dir)}
	batch_file = open(output_files['batch'], 'w')
		
	for cell in reads_assigned_db.keys():
		cell_name = 'cell_%s' % cell.decode('utf-8')
		print('\tWorking on cell %s' % cell_name)
		
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
		
		cell_offsets = IO_utils.get_from_db(reads_assigned_pipe, [cell])[0]
		assert len(cell_offsets) % 2 == 0, \
			'Cell offsets must contain an even number of entries'
		reads_iter = IO_utils.read_fastq_random(
			reads_unzipped, 
			[cell_offsets[i] for i in range(len(cell_offsets)) if i % 2 == 0])
		barcodes_iter = IO_utils.read_fastq_random(
			barcodes_unzipped,
			[cell_offsets[i] for i in range(len(cell_offsets)) if i % 2 == 1])
		
		while(True):
			try:
				reads_data, _ = next(reads_iter)
				barcodes_data, _ =  next(barcodes_iter)
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
	batch_file.close()
	return output_files


def get_args():
	parser = argparse.ArgumentParser(
		description = 'This script splits reads for dropseq data')
	parser.add_argument('--barcodes', 
		type=str, 
		help='Barcodes file name (unzipped)', 
		required=True)
	parser.add_argument('--reads', 
		type=str, 
		help='RNAseq reads file name (unzipped)', 
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
		default=3)
	parser.add_argument('--breadth', 
		type=int, 
		help='How many nodes search.', 
		default=2000)
	parser.add_argument('--threads', 
		type=int, 
		help='Number of threads to use.', 
		default=32)
	parser.add_argument('--min_dist', 
		type=int, 
		help='Minimum Hamming distance between barcodes.', 
		default=3)
	parser.add_argument('--index_depth',
		type=float,
		help='Fraction of reads to build kmer index from',
		default=0.1)
	
	return vars(parser.parse_args())

if __name__ == '__main__':
	cmdline_args = get_args()	
	output_files, elapsed_time = run_all(cmdline_args)
	print('Done. Time elapsed: %f seconds' % elapsed_time)



