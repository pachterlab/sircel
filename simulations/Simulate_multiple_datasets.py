"""

Prepare all simulated datasets

"""

import numpy as np
import sys
import os
import gzip
import json
import itertools
import subprocess

def run_simulations():
	ALPHABET = ['A', 'C', 'G', 'T']
	BARCODE_LENGTH = 12
	UMI_LENGTH = 8
	NUM_READS = 100000
	NUM_BARCODES = 500
	NUM_REPS = 3
	
	abundances = ['normal', 'uniform', 'exponential']
	error_types = ['any', 'mismatch', 'insertion', 'deletion']
	poiss_errors = list(range(0, 5))
	
	output_dir = sys.argv[1]
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	summary = open('%s/summary.txt' % output_dir, 'w')
	summary.write(
		'Simulation\tReplicate\tAbundances\tError_type\tPoiss_error\tOutput_dir\n')
	count = 1
	for rep in range(0, NUM_REPS):
		np.random.seed(rep)
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
							
			printer = [str(i) for i in \
				[count, rep, abundance, error_type, poiss_error, simulation_dir]]
			summary.write('\t'.join(printer) + '\n')
			
			print('\t%s abundances\n\t%s error type\n\t%i errors per read' % vals)	
			
			print('\tSimulating ground truth')
			barcode_abundances = get_barcodes_abundance(
				NUM_BARCODES, abundance_distr = abundance)
						
			write_barcodes(true_barcodes, barcode_abundances, simulation_dir)
			print('\tSimulating random reads')
			write_reads((true_barcodes,
				barcode_abundances,
				error_type,
				poiss_errors,
				NUM_READS,
				BARCODE_LENGTH,
				ALPHABET,
				UMI_LENGTH,
				simulation_dir))
			print('\tSplitting reads')
			run_splitter(simulation_dir)
			count += 1
	summary.close()	
	
def run_splitter(simulation_dir):
	split_reads_master = sys.argv[2]
	
	bc_file = '%s/barcodes.fastq.gz' % simulation_dir
	output_dir = '%s/barcodes_split' % simulation_dir
	
	cmd = ['python3', split_reads_master, 
		'--dropseq',
		'--reads', bc_file,
		'--barcodes', bc_file,
		'--output_dir', simulation_dir,
		'--threads', '32',
		'--dropseq']
	
	split_reads = subprocess.Popen(cmd)
	_ = split_reads.communicate()[0]	

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
	if(num_errors == 0):
		return seq, num_errors
	
	for i in range(num_errors):
		seq = add_single_error(seq, BARCODE_LENGTH, ALPHABET, error_type=error_type)
	return seq, num_errors

def write_reads(params):
	(barcodes,
		barcode_abundances,
		err_type,
		poiss_errors,
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
				true_barcode, BARCODE_LENGTH, ALPHABET, RATE = poiss_errors, error_type=err_type)
			umi = ''.join(np.random.choice(ALPHABET, UMI_LENGTH))
			seq = mutated_barcode + umi
			qual = ''.join(['I' for i in seq])
		
			read_name = 'ReadNum:%i_NumErr:%i_TrueBarcode:%s' % (i, num_errors, true_barcode)
			fq_line = '%s\n%s\n%s\n%s\n' % (read_name, seq, '+', qual)
		
			writer.write(fq_line.encode('utf-8'))

def write_barcodes(barcodes, abundances, output_dir):
	true_barcodes_file = '%s/true_barcodes.txt' % output_dir
	with open(true_barcodes_file, 'w') as writer:
		for (bc, abund) in zip(barcodes, abundances):
			writer.write('%s\t%f\n' % (bc, abund))


if __name__ == "__main__":
	run_simulations()
	




	
	