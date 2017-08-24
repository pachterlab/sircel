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
import os
import sys
import time
import json
import collections
import itertools
import gc
import numpy as np
from multiprocessing import Pool 
from sircel import IO_utils
from sircel.Graph_utils import Edge, Graph, Path

import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt

np.random.seed(0)
args = {}
output_files = {}
output_dir = ''

def run_all(cmdline_args):
	print('Running Split_reads')
	
	global args
	global output_files
	global output_dir
	
	args = cmdline_args
	output_dir = args['output_dir']
	output_files['log'] = '%s/run_log.txt' % output_dir
	sys.stdout = IO_utils.Logger(output_files['log'])
	
	start_time = time.time()
	
	reads_unzipped = args['reads']
	barcodes_unzipped = args['barcodes']
	reads_offsets = args['reads_offsets']
	barcodes_offsets = args['barcodes_offsets']
	
	print('Indexing reads by circularized kmers')
	kmer_index, kmer_counts = get_kmer_index(barcodes_unzipped, list(barcodes_offsets))
		#list(barcodes_offsets) makes a copy of barcodes_offsets. 
		#this allows use of the pop() method without modifying the old list
	print('\t%i unique kmers indexed' % len(kmer_counts.items()))
	
	print('Finding cyclic paths')
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
	(threshold, top_paths, fit_out) = threshold_paths(output_dir, cyclic_paths)
	output_files.update(fit_out)
	output_files['thresholded_paths'] = IO_utils.save_paths_text(
		output_dir, top_paths, prefix='threshold')
	
	print('Assigning reads')
	reads_assigned_db, reads_assigned_pipe = assign_all_reads(
		(top_paths,
		reads_unzipped, 
		reads_offsets, 
		barcodes_unzipped, 
		barcodes_offsets))
	
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
	
def get_kmer_index(barcodes_unzipped, barcodes_offsets):
	PEARSONR_CUTOFF = 0.999
	MIN_ITERS = 5
	
	length = args['barcode_end'] - args['barcode_start']
	pool = Pool(processes = args['threads'])
	
	read_count = 0
	kmer_idx = {}
	counts_corr_coefs = []
	num_reads = []
	
	for (chunk_num, reads_chunk) in enumerate(
		IO_utils.get_read_chunks(
			barcodes_unzipped,
			barcodes_offsets,
			random = True,
			BUFFER_SIZE = 10000)):
				
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
		
		del(chunk_kmer_indices); _ = gc.collect()
		
		new_kmer_counts = get_kmer_counts(kmer_idx)
		#check kmer count correlation
		counts_corr_coef = get_kmer_count_correlation(
			old_kmer_counts, new_kmer_counts)
		counts_corr_coefs.append(counts_corr_coef)
		print('\t%i reads indexed. %f running pearsonr' % \
				(read_count, counts_corr_coef))
		
		if(len(counts_corr_coefs) >= MIN_ITERS):
			if(all(corr >= PEARSONR_CUTOFF for corr in counts_corr_coefs[-MIN_ITERS:])):
				break
		
	plot_kmer_subsamp_pearson(counts_corr_coefs, num_reads)
	return kmer_idx, new_kmer_counts

def index_read(params):
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
	
def plot_kmer_subsamp_pearson(counts_corr_coefs, num_reads):
	fig, ax = plt.subplots(
		nrows = 1, 
		ncols = 1,
		figsize = (4,4))
	
	ax.plot(num_reads, counts_corr_coefs)
	ax.set_xlabel('Number of reads indexed')
	ax.set_ylabel('Pearson R')
	ax.set_title('Correlation between relative kmer counts \nas number of reads are incrementally increased')
	ax.grid()
	ax.set_ylim([0,1])
	fig.savefig('%s/indexed_reads_correlation.png' % output_dir, dpi = 300)
	
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
			itertools.repeat(barcodes_unzipped),
			itertools.repeat(barcode_length)))
		paths += [item for sublist in paths_group for item in sublist]
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
		paths.append( (seq, weight, counter, nodes) )
		if(counter > args['depth']):
			break
		counter += 1
	
	return merge_paths(paths)

def build_subgraph(reads_in_subgraph, barcodes_unzipped):
	barcodes_iter = IO_utils.read_fastq(
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

def threshold_paths(output_dir, paths):
	from scipy import signal
	
	NUM_CELLS = [100, 5000]
	LOCAL_WINDOW_LEN = 50
	MIN_CAPACITY = 10
	
	threshold_out = {
		'slopes' : '%s/slopes.txt' % output_dir,
		'paths_plot' : '%s/paths_plotted.png' % output_dir,
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
	unique_paths_sorted = sorted(unique_paths.values(), key = lambda tup: tup[1], reverse = True)
	
	fig, ax = plt.subplots(
		nrows = 1, 
		ncols = 3,
		figsize = (12,4))

	y = [tup[1] for tup in unique_paths_sorted]
	x = range(len(y))
	ax[0].step(x, y, label='Path weights')
	ax[0].set_yscale('log')

	grad = [-1 * i for i in \
				local_lin_fit(np.log10(y), window_len=LOCAL_WINDOW_LEN)]   
	x = [LOCAL_WINDOW_LEN / 2 + i for i in range(len(grad))]
	ax[1].plot(x, grad, color='b', label='First derivative (smoothed)')
	second_grad = local_lin_fit(grad, window_len = LOCAL_WINDOW_LEN)
	x = [LOCAL_WINDOW_LEN + i for i in range(len(second_grad))]
	ax[2].plot(x, second_grad, color = 'b', label = 'Second derivative (smoothed)')

	lmax = get_lmax(second_grad, LOCAL_WINDOW_LEN)
	threshold = get_threshold(grad, lmax)
	print('\tThreshold is %i' % threshold)
	top_paths = unique_paths_sorted[0:threshold]
	
	for a in ax:
		for m in lmax:
			a.axvline(m, color = 'gray', alpha = 0.25)
		a.axvline(threshold, color='k', label='Threshold')
		a.legend(loc = 1, fontsize = 6)
		a.grid()
		#a.set_xscale('log')
		a.set_xlim([0, 1000])
		a.set_xlabel('Path (sorted)')
	
	plt.tight_layout()
	fig.savefig(threshold_out['paths_plot'], dpi = 300)
	print('\tMerging similar paths by Hamming distance')
	top_paths = merge_paths(top_paths)#merges by hamming distance
	print('\t%i paths remaining after merging' % len(top_paths))
	return len(top_paths), top_paths, threshold_out

def get_lmax(second_grad, LOCAL_WINDOW_LEN):
	#finds zeros in
	lmax = []
	for i in range(len(second_grad) - 1):
		if(second_grad[i] > 0 and second_grad[i + 1] <= 0):
			lmax.append(int(i + LOCAL_WINDOW_LEN))
	return lmax

def get_threshold(grad, lmax):
	threshold = lmax[0]
	for z in lmax:
		if(grad[z] > grad[threshold]):
			threshold = z
	return threshold

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
		
def assign_all_reads(params):
	(	top_paths,
		reads_unzipped, 
		reads_offsets, 
		barcodes_unzipped, 
		barcodes_offsets) = params
	
	MIN_KMER_SIZE = 4
	MAX_KMER_SIZE = args['barcode_end'] - args['barcode_start']
	
	#initialize vars
	reads_assigned_db, reads_assigned_pipe = IO_utils.initialize_redis_pipeline(db=0)
	kmers_to_paths = {}
	
	#print('\tGetting kmers in paths')
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

	#print('\tAssigning reads to paths')
	pool = Pool(processes = args['threads'])	
	read_count = 0
	num_unassigned = 0
	BUFFER_SIZE = 100000	
	
	for reads_chunk, barcodes_chunk in zip(
		IO_utils.get_read_chunks(
			reads_unzipped,
			reads_offsets,
			random = False,
			BUFFER_SIZE = BUFFER_SIZE),
		IO_utils.get_read_chunks(
			barcodes_unzipped,
			barcodes_offsets,
			random = False,
			BUFFER_SIZE = BUFFER_SIZE)):
		
		read_count += len(reads_chunk)
		assignments = pool.map(assign_read, 
			zip(itertools.repeat(kmers_to_paths),
			itertools.repeat(MIN_KMER_SIZE),
			itertools.repeat(MAX_KMER_SIZE),
			reads_chunk,
			barcodes_chunk))
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
		kmers_to_paths: dict of kmer -> list of paths that contain it
		min_kmer_size
		max_kmer_size
		read: list of fastq entry lines
	returns
	
	"""
	(kmers_to_paths,
		min_kmer_size,
		max_kmer_size,
		(reads_data, reads_offset),
		(barcodes_data, barcodes_offset)) = params
	
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
		
		try:
			cell_offsets = IO_utils.get_from_db(reads_assigned_pipe, [cell])[0]
		except IndexError:
			pass
			
			
		assert len(cell_offsets) % 2 == 0, \
			'Cell offsets must contain an even number of entries'
		reads_iter = IO_utils.read_fastq(
			reads_unzipped, 
			[cell_offsets[i] for i in range(len(cell_offsets)) if i % 2 == 0])
		barcodes_iter = IO_utils.read_fastq(
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
	
	return vars(parser.parse_args())

if __name__ == '__main__':
	cmdline_args = get_args()	
	output_files, elapsed_time = run_all(cmdline_args)
	print('Done. Time elapsed: %f seconds' % elapsed_time)



