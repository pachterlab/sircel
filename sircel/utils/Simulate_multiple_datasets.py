"""

Prepare all simulated datasets

"""

import numpy as np
np.random.seed(0)
import sys
import os
import gzip
import json
import itertools
import subprocess

from sircel.utils import IO_utils
from sircel.Sircel_master import get_args, run_all
from sircel.utils.Naive_pipeline import run_naive_pipeline

def evaluate_simulations(summary_file):
	output_dir = sys.argv[1]
	
	#read summary file
	summary_data = []
	with open (summary_file, 'r') as inf:
		header = inf.readline().strip().split('\t')
		for line in inf:
			row = line.strip().split('\t')
			summary_data.append(row)

	#for each file: get true and pred barcodes
	entry_to_col = {key:col for (col, key) in enumerate(header)}
	get_col = lambda key: entry_to_col[key]
	
	new_header = header + [
		'num_true_pos_bc', 'num_false_pos_bc', 'num_fals_neg_bc', 'bc_tpr', 'bc_fpr']
	
	summary_processed_file = '%s/summary_processed.txt' % output_dir
	writer = open(summary_processed_file, 'w')
	writer.write('\t'.join(new_header) + '\n')
	
	for simulation_entry in summary_data:
		simulation_dir = simulation_entry[get_col('Output_dir')]
		
		(num_tp,		#number of true pos barcodes
			num_fp,	#number of false pos barcodes
			num_fn,	#number of false neg barcodes
			bc_tpr,	#true positive rate (list) of barcode assignments
			bc_fpr) = eval_single_file(simulation_dir)
		#save all this data
		bc_tpr_str = ','.join([str(i) for i in bc_tpr])
		bc_fpr_str = ','.join([str(i) for i in bc_fpr])
		simulation_entry += [
			num_tp, num_fp, num_fn, bc_tpr_str, bc_fpr_str]
		printer = ('\t'.join([str(i) for i in simulation_entry]))
		writer.write(printer + '\n')
	writer.close()
	print('done')
	
	return summary_processed_file

def eval_single_file(simulation_output_dir):
	sim_dat_dir = simulation_output_dir[0:simulation_output_dir.rindex('/')]
	true_barcodes = get_barcodes_set(
		'%s/true_barcodes.txt' % (sim_dat_dir))
	pred_barcodes = get_barcodes_set(
		'%s/threshold_paths.txt' % simulation_output_dir)
	true_positives = true_barcodes & pred_barcodes
	
	num_tp = len(true_positives)				#number of false positive bc
	num_fp = len(pred_barcodes) - num_tp	#number of true positive bc
	num_fn = len(true_barcodes) - num_tp	#numberof false negative bc
	
	#print(num_tp, num_fp, num_fn)
	
	#for each true positive barcode, get fraction of correct reads
	barcodes_tpr = []
	barcodes_fpr = []
	for bc in pred_barcodes:
		(bc_tpr, bc_fpr) = get_fraction_correct_reads(
			bc, simulation_output_dir)
		barcodes_tpr.append(bc_tpr)
		barcodes_fpr.append(bc_fpr)
	
	return (num_tp,
		num_fp,
		num_fn,
		barcodes_tpr,
		barcodes_fpr)

def get_barcodes_set(true_bc_file):
	barcodes = set()
	try:
		inf = open(true_bc_file, 'r')
	except FileNotFoundError:
		return barcodes
	for line in inf:
		bc = line.strip().split('\t')[0]
		barcodes.add(bc)
	inf.close()
	return barcodes

def get_fraction_correct_reads(pred_bc, simulation_output_dir):
	fq_fname = '%s/reads_split/cell_%s_barcodes.fastq.gz' % \
		(simulation_output_dir, pred_bc)
	fq_file = gzip.open(fq_fname, 'rb')
	if not os.path.exists(fq_file):
		return (0,0)
	fq_iter = IO_utils.read_fastq_sequential(fq_file)
	
	tpr = 0.
	fpr = 0.
	for (lines, _) in fq_iter:
		read_name = lines[0]
		assigned_bc = read_name.split(':')[-1]
		true_bc = read_name.split(':')[-2].split('_')[0]	
		
		if(assigned_bc == true_bc):
			tpr += 1.
		else:
			fpr += 1.
	total_reads = tpr + fpr
	tpr /= total_reads
	fpr /= total_reads
	
	return(tpr, fpr)

def run_simulations():
	ALPHABET = ['A', 'C', 'G', 'T']
	BARCODE_LENGTH = 12
	UMI_LENGTH = 8
	NUM_READS = 100000
	NUM_BARCODES = 500
	NUM_REPS = 3
	
	abundances = ['normal', 'uniform', 'exponential']
	error_types = ['any', 'mismatch', 'insertion', 'deletion']
	poiss_errors = list(range(0, 4))
	
	output_dir = sys.argv[1]
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
		
	summary_file = '%s/summary.txt' % output_dir
	summary = open(summary_file, 'w')
	summary.write(
		'Simulation\tReplicate\tAbundances\tError_type\tPoiss_error\tOutput_dir\tPipeline\n')
	count = 1
	for rep in range(0, NUM_REPS):
		true_barcodes = get_barcodes(
			NUM_BARCODES,
			BARCODE_LENGTH,
			ALPHABET)
		
		for (vals) in itertools.product(abundances, error_types, poiss_errors):
			print('Working on simulation %i' % count)
			(abundance, error_type, poiss_error) = vals
			simulation_dir = '%s/simulation_%i' % (output_dir, count)
			if not os.path.exists(simulation_dir):
				os.makedirs(simulation_dir)
			
			entries = [
				count,
				rep,
				abundance,
				error_type,
				poiss_error,
				None,
				None]
			
			entries[5] = simulation_dir + '/sircel_kmers'
			entries[6] = 'sircel_kmers'
			summary.write('\t'.join([str(i) for i in entries]) + '\n')
			
			entries[5] = simulation_dir + '/sircel_lev'
			entries[6] = 'sircel_lev'
			summary.write('\t'.join([str(i) for i in entries]) + '\n')
			
			entries[5] = simulation_dir + '/naive'
			entries[6] = 'naive'
			summary.write('\t'.join([str(i) for i in entries]) + '\n')
			
			print('\t%s abundances\n\t%s error type\n\t%i errors per read' % vals)	
			
			print('\tSimulating ground truth')
			barcode_abundances = get_barcodes_abundance(
				NUM_BARCODES, abundance_distr = abundance)
						
			write_barcodes(true_barcodes, barcode_abundances, simulation_dir)
			print('\tSimulating random reads')
			write_reads((true_barcodes,
				barcode_abundances,
				error_type,
				poiss_error,
				NUM_READS,
				BARCODE_LENGTH,
				ALPHABET,
				UMI_LENGTH,
				simulation_dir))
			print('\tRunning sircel using kmers to assign reads')
			run_sircel_kmers(simulation_dir)
			print('\tRunning sircel using Levenshtein distance to assign reads')
			run_sircel_lev(simulation_dir)
			print('\tRunning naive pipeline to split reads')
			run_naive(simulation_dir)
			print('\n\n')
			
			count += 1
	summary.close()	
	
	return summary_file
	
def run_sircel_kmers(simulation_dir):
	bc_file = '%s/barcodes.fastq.gz' % simulation_dir
	output_dir = '%s/sircel_kmers' % simulation_dir
	
	args = get_args([ 
		'--dropseq',
		'--reads', bc_file,
		'--barcodes', bc_file,
		'--output_dir', output_dir,
		'--threads', '32'])
	_ = run_all(args)

def run_sircel_lev(simulation_dir):
	bc_file = '%s/barcodes.fastq.gz' % simulation_dir
	output_dir = '%s/sircel_lev' % simulation_dir
	
	args = get_args([ 
		'--dropseq',
		'--reads', bc_file,
		'--barcodes', bc_file,
		'--output_dir', output_dir,
		'--threads', '32',
		'--split_levenshtein', 'True'])
	_ = run_all(args)

def run_naive(simulation_dir):
	bc_file = '%s/barcodes.fastq.gz' % simulation_dir
	output_dir = '%s/naive' % simulation_dir
	
	_ = run_naive_pipeline(bc_file, bc_file, output_dir)

def get_barcodes(NUM_BARCODES, BARCODE_LENGTH, ALPHABET):
	barcodes = []
	for i in range(NUM_BARCODES):
		barcode = np.random.choice(ALPHABET, BARCODE_LENGTH)
		barcodes.append(''.join(barcode))
	return barcodes

def get_barcodes_abundance(NUM_BARCODES, abundance_distr = None):
	if(abundance_distr == None or abundance_distr == 'flat'):
		abundances = [1 for i in range(NUM_BARCODES)]
	elif(abundance_distr == 'uniform'):
		abundances = np.random.uniform(0, 1, NUM_BARCODES)
	elif(abundance_distr == 'normal'):
		abundances = np.random.normal(10, 1, NUM_BARCODES)
		offset = np.min(abundances) + 1 #ensure all values are positive
		abundances += offset
	elif(abundance_distr == 'exponential'):
		abundances = np.random.exponential(0.2, NUM_BARCODES)
		
	abundances = np.fabs(abundances)
	abundances /= np.sum(abundances)#normalize
	
	return abundances

def add_single_error(true_barcode, BARCODE_LENGTH, ALPHABET, error_type = None):
	#returns true_barcode sequence with one SNP, INSERTION, or DELETION
	mutated_barcode = list(true_barcode)
	mutation_pos = np.random.choice(BARCODE_LENGTH)
	if(error_type == 'any'):
		mutations = ['mismatch', 'insertion', 'deletion']
	else:
		mutations = [error_type]
	mutation_choice = np.random.choice(mutations)
	
	if(mutation_choice == 'mismatch'):#snp
		mutated_barcode[mutation_pos] = np.random.choice(ALPHABET)
	elif(mutation_choice == 'deletion'):#deletion
		mutated_barcode = true_barcode[0:mutation_pos] + \
										true_barcode[mutation_pos+1:] + \
										np.random.choice(ALPHABET)#insert a random nucleotide
	else:#insertion
		mutated_barcode = true_barcode[0:mutation_pos] + \
										np.random.choice(ALPHABET) + \
										true_barcode[mutation_pos:]
		mutated_barcode = mutated_barcode[0:-1]
	mutated_barcode = ''.join(mutated_barcode)
	#print(('%s\t%i\n%s\n%s\n') % \
	#	(mutation_choice, mutation_pos, true_barcode, mutated_barcode))
	return(mutated_barcode)

def add_multiple_errors(true_barcode, BARCODE_LENGTH, ALPHABET, RATE=2, error_type='any'):
	"""
	Adds a poisson number of errors 
	"""
	seq = true_barcode
	num_errors = np.random.poisson(lam=RATE)
	for i in range(num_errors):
		seq = add_single_error(seq, BARCODE_LENGTH, ALPHABET, error_type=error_type)
	return seq, num_errors

def write_reads(params):
	(barcodes,
		barcode_abundances,
		err_type,
		poiss_error,
		num_reads,
		BARCODE_LENGTH,
		ALPHABET,
		UMI_LENGTH,
		output_dir) = params
	
	barcodes_file = '%s/barcodes.fastq.gz' % output_dir
	with gzip.open(barcodes_file, 'wb') as writer:
		for i in range(num_reads):
			true_barcode = np.random.choice(barcodes, p = barcode_abundances)
		
			mutated_barcode, num_errors = add_multiple_errors(
				true_barcode, BARCODE_LENGTH, ALPHABET, RATE = poiss_error, error_type=err_type)
			umi = ''.join(np.random.choice(ALPHABET, UMI_LENGTH))
			seq = mutated_barcode + umi
			qual = ''.join(['I' for i in seq])
		
			read_name = '@ReadNum:%i_NumErr:%i_TrueBarcode:%s' % (i, num_errors, true_barcode)
			fq_line = '%s\n%s\n%s\n%s\n' % (read_name, seq, '+', qual)
		
			writer.write(fq_line.encode('utf-8'))

def write_barcodes(barcodes, abundances, output_dir):
	true_barcodes_file = '%s/true_barcodes.txt' % output_dir
	with open(true_barcodes_file, 'w') as writer:
		for (bc, abund) in zip(barcodes, abundances):
			writer.write('%s\t%f\n' % (bc, abund))


if __name__ == "__main__":
	print("Running simulations")
	#summary_file = run_simulations()
	summary_file = '%s/summary.txt' % sys.argv[1]
	print("Evaluating simulations")
	summary_processed_file = evaluate_simulations(summary_file)
	




	
	
