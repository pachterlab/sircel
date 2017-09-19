"""
Naive pipeline
	Count kmers (full-length barcode seqs)
	Threshold them (use vasilis approach, including pitfalls)
	Assign reads to thresholded kmers by Lev or Ham dist
"""
import numpy as np
import sys
from sircel.Split_reads import *
from sircel import IO_utils
from sircel import Plot_utils

from scipy.signal import savgol_filter as savgol
from multiprocessing import Pool
from Levenshtein import distance, hamming
from itertools import repeat

args = {}
def run_naive_pipeline(barcodes, reads, output_dir):
	global args
	output_files = {}
	
	args['barcodes'] = barcodes
	args['reads'] = reads
	args['output_dir'] = output_dir
	
	if not os.path.exists(args['output_dir']):
		os.makedirs(args['output_dir'])
	if not os.path.exists(args['output_dir'] + '/plots'):
		os.makedirs(args['output_dir'] + '/plots')
	
	args['barcode_start']	= 0
	args['barcode_end'] 		= 12
	args['umi_start']			= 12
	args['umi_end'] 			= 20
	args['threads'] 			= 4
	
	print('Unzipping files (temporary)')
	reads_unzipped = \
		IO_utils.unzip(args['reads'].split(','))
	barcodes_unzipped = \
		IO_utils.unzip(args['barcodes'].split(','))
	args['reads'] = reads_unzipped
	args['barcodes'] = barcodes_unzipped
	
	print('Counting kmers')
	kmer_counts = count_kmers(
		args['barcodes'],
		args['barcode_start'],
		args['barcode_end'])
	
	print('Thresholding kmers')
	thresholded_bcs, plt = threshold_bcs(kmer_counts, args['output_dir'])
	print('%i kmers above threshold' % len(thresholded_bcs))
	#write thresholded kmers to file
	output_files['threshold_paths'] = '%s/threshold_paths.txt' % args['output_dir']
	with open(output_files['threshold_paths'], 'w') as f:
		printer = '\n'.join([i for i in thresholded_bcs])
		writer.write(printer)

	print('Assigning reads')
	reads_assigned = assign_all_reads((
		thresholded_bcs,
		reads_unzipped, 
		barcodes_unzipped))
	
	print('Writing split fq files')
	
	output_files['split'] = write_split_fastqs(
		(reads_assigned,
		args['output_dir'],
		reads_unzipped,
		barcodes_unzipped))
	output_files['plt'] = plt
	output_files['run_outputs'] = '%s/run_outputs.json' % args['output_dir']
	
	import json
	with open(output_files['run_outputs'], 'w') as writer:
		writer.write(json.dumps(output_files, indent=3))
	
	print('Done')

def count_kmers(barcodes_unzipped, barcode_start, barcode_end):
	BUFFER_SIZE = 10000
	
	pool = Pool(processes = args['threads'])
	read_count = 0
	kmer_counts = {}
	
	barcodes_f = open(barcodes_unzipped, 'rb')
	for (chunk_num, reads_chunk) in enumerate(
		IO_utils.get_read_chunks(
			barcodes_f,
			random = False,
			BUFFER_SIZE = BUFFER_SIZE)):
					
		read_count += len(reads_chunk)
	
		for (read, _) in reads_chunk:
			kmer = read[1][barcode_start : barcode_end]
			if kmer not in kmer_counts:
				kmer_counts[kmer] = 0
			kmer_counts[kmer] += 1
	barcodes_f.close()
	pool.close()
	return kmer_counts
	
def threshold_bcs(kmer_counts, output_dir):
	LOCAL_WINDOW_LEN = 50
	
	#convert to tuple, sort
	kmer_counts_lst = sorted(
		list(kmer_counts.items()), 
		key = lambda tup: tup[1],
		reverse = True)	
	y = [tup[1] for tup in kmer_counts_lst]
	x = list(range(len(y)))
	
	print('\tComputing gradient')
	grad = [-1 * i for i in \
				local_lin_fit(np.log10(y),
				window_len=LOCAL_WINDOW_LEN)]  
	print('\tComputing second gradient') 
	second_grad = local_lin_fit(grad, window_len = LOCAL_WINDOW_LEN)
	lmax = get_lmax(second_grad, LOCAL_WINDOW_LEN)
	threshold = get_threshold(grad, lmax, None, kmer_counts_lst)
	top_paths = kmer_counts_lst[0:threshold]
	
	print('\tThreshold is %i' % threshold)
	
	plt_name = Plot_utils.plot_path_threshold(
		(output_dir, 
			y,
			grad,
			second_grad,
			lmax,
			threshold,
			LOCAL_WINDOW_LEN))
	
	return set(tup[0] for tup in top_paths), plt_name

def assign_all_reads(params):
	(	consensus_bcs,
		reads_unzipped, 
		barcodes_unzipped) = params
	
	BUFFER_SIZE = 100000
	pool = Pool(processes = args['threads'])
	reads_assigned = {}
		#key / value map of: [cell name] :-> list of line offsets
	for bc in consensus_bcs:
		reads_assigned[bc] = []
	reads_assigned['unassigned'] = []
	
	print('\tAssigning reads')
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
		
		assignments = pool.map(assign_read_levenshtein, 
			zip(
			repeat(consensus_bcs),
			reads_chunk,
			barcodes_chunk))
		
		for (assignment, offset1, offset2) in assignments:
			if(assignment == 'unassigned'):
				num_unassigned += 1
			reads_assigned[assignment].append((offset1, offset2))
		print('\tProcessed %i reads' % read_count)
	pool.close()
	print('\t%i reads could not be assigned' % num_unassigned)
	return reads_assigned
	
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

def write_split_fastqs(params):
	import gzip
	(	reads_assigned,
		output_dir,
		reads_unzipped,
		barcodes_unzipped) = params
	
	split_dir = '%s/reads_split' % output_dir
	if not os.path.exists(split_dir):
		os.makedirs(split_dir)
	output_files = {'batch' : '%s/batch.txt' % (split_dir)}
	batch_file = open(output_files['batch'], 'w')
		
	for cell, cell_offsets in reads_assigned.items():
		
		cell_name = 'cell_%s' % cell
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
			offsets = [tup[0] for tup in cell_offsets])
		barcodes_iter = IO_utils.read_fastq_random(
			barcodes_f,
			offsets = [tup[1] for tup in cell_offsets])
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
	batch_file.close()
	return output_files

	
	
	
if __name__ == '__main__':
	run_naive_pipeline(sys.argv[1], sys.argv[2], sys.argv[3])
	
	
	
	
	
	
	
	