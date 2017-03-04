"""

Prepare all simulated datasets

"""

import numpy as np
import sys
import os
import gzip
import json


def run_simulation(params):
	(output_dir,
	num_reads,
	num_barcodes,
	poiss_err,
	err_type,
	barcode_abundance_dist,
	seed,
	expt_name) = params
	
	ALPHABET = ['A', 'C', 'G', 'T']
	BARCODE_LENGTH = 12
	UMI_LENGTH = 8
	np.random.seed(seed)

	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	true_barcodes = get_barcodes(
		num_barcodes,
		BARCODE_LENGTH,
		ALPHABET)
	barcode_abundances = get_barcodes_abundance(
		num_barcodes, abundance_distr = barcode_abundance_dist)
	write_reads((true_barcodes,
		barcode_abundances,
		err_type,
		num_reads,
		BARCODE_LENGTH,
		ALPHABET,
		UMI_LENGTH,
		output_dir))
	write_barcodes(true_barcodes, barcode_abundances, output_dir)
	
	simulation_params = {
		'num_reads' : num_reads,
		'num_barcodes' : num_barcodes,
		'poiss_err' : poiss_err,
		'err_type' : err_type,
		'barcode_abundance_dist' : barcode_abundance_dist,
		'seed' : seed,
		'output_dir' : output_dir,
		'experiment' : expt_name}
	
	return simulation_params

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
		abundances = np.random.normal(0, 1, NUM_BARCODES)
		abundances += min(abundances) + 1 #ensure all values are positive
	elif(abundance_distr == 'exponential'):
		abundances = np.random.exponential(0, 1, NUM_BARCODES)
	
	abundances /= np.sum(abundances)#normalize
	assert(all(abundances) > 0), 'Abundance must be greater than zero'
		#this should never happen
	
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
	#print(('%s\t%i\n%s\n%s\n') % (mutation_choice, mutation_pos, true_barcode, mutated_barcode))
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
		num_reads,
		BARCODE_LENGTH,
		ALPHABET,
		UMI_LENGTH,
		output_dir) = params
	
	
	barcodes_file = '%s/barcodes.fastq.gz' % output_dir
	with gzip.open(barcodes_file, 'wb') as writer:
		for i in range(num_reads):
			true_barcode = np.random.choice(barcodes, p=barcode_abundances)
		
			mutated_barcode, num_errors = add_multiple_errors(
				true_barcode, BARCODE_LENGTH, ALPHABET, error_type=err_type)
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



base_dir = sys.argv[1]
all_simulations = []
NUM_REPLICATES = 5



#1. flat abundances, variable poisson error, all error type

print('Simulating varied errors per read')
expt_name = 'vary_poiss_error'
num_reads = 100000
num_barcodes = 1000
err_type = 'any'
abundane_distr = 'flat'

for poiss_err in range(0, 5):
	for replicate in range(NUM_REPLICATES):
		output_dir = '%s/%s_%i_rep_%i' % (base_dir, expt_name, poiss_err, replicate + 1)
		simulation_params = run_simulation(
			(output_dir,
			num_reads,
			num_barcodes,
			poiss_err,
			err_type,
			abundane_distr,
			replicate,
			expt_name))
		all_simulations.append(simulation_params)




#2. flat abundances, fixed poisson = 1, vary error types

print('Simulating varied error types')
expt_name = 'vary_error_type'
num_reads = 100000
num_barcodes = 1000
poiss_err = 1
abundane_distr = 'flat'

for err_type in ['mismatch', 'insertion', 'deleteion', 'all']:
	for replicate in range(NUM_REPLICATES):
		output_dir = '%s/%s_%s_rep_%i' % (base_dir, expt_name, err_type, replicate + 1)
		simulation_params = run_simulation(
			(output_dir,
			num_reads,
			num_barcodes,
			poiss_err,
			err_type,
			abundane_distr,
			replicate,
			expt_name))
		all_simulations.append(simulation_params)

	
	
#3. vary abundance distr, fixed poisson = 1, fixed error types = all	

print('Simulating varied barcodes abundances')
expt_name = 'vary_abundance_distr'
num_reads = 100000
num_barcodes = 1000
poiss_err = 1
err_type = 'any'

for abundance_distr in ['flat', 'normal', 'exponential']:
	for replicate in range(NUM_REPLICATES):
		output_dir = '%s/%s_%s_rep_%i' % (base_dir, expt_name, abundance_distr, replicate + 1)
		simulation_params = run_simulation(
			(output_dir,
			num_reads,
			num_barcodes,
			poiss_err,
			err_type,
			abundane_distr,
			replicate,
			expt_name))
		all_simulations.append(simulation_params)



#4. vary number of barcodes, fixed poisson = 1, fixed error types = all, fixed abundance distr = flat

print('Simulating varied number of barcodes')
expt_name = 'vary_num_barcodes'
num_reads = 100000
poiss_err = 1
err_type = 'any'
abundance_distr = 'flat'

for num_barcodes in [10, 100, 1000, 10000]:
	for replicate in range(NUM_REPLICATES):
		output_dir = '%s/%s_%i_rep_%i' % (base_dir, expt_name, num_barcodes, replicate + 1)
		simulation_params = run_simulation(
			(output_dir,
			num_reads,
			num_barcodes,
			poiss_err,
			err_type,
			abundane_distr,
			replicate,
			expt_name))
		all_simulations.append(simulation_params)

	
#5. vary number of reads, fixed poisson = 1, fixed error types = all, fixed abundance distr = flat, fixed barcodes = 100

print('Simulating varied number of reads')
expt_name = 'vary_num_reads'
num_barcodes = 1000
poiss_err = 1
err_type = 'any'
abundance_distr = 'flat'

for num_reads in [1000, 10000, 100000, 1000000]:
	for replicate in range(NUM_REPLICATES):
		output_dir = '%s/%s_%i_rep_%i' % (base_dir, expt_name, num_reads, replicate + 1)
		simulation_params = run_simulation(
			(output_dir,
			num_reads,
			num_barcodes,
			poiss_err,
			err_type,
			abundance_distr,
			replicate,
			expt_name))
		all_simulations.append(simulation_params)


#write output summary file
print(all_simulations)

with open('%s/summary_files.txt' % base_dir, 'w') as w:
	for (i, simulation) in enumerate(all_simulations):
		keys = sorted(list(simulation.keys()))
		if(i == 0):
			w.write('\t'.join(keys) + '\n')
		w.write('\t'.join(str(simulation[k]) for k in keys) + '\n')

print('Done')









	
	