"""

simulate barcode fastq

"""
import numpy as np
np.random.seed(0)
import sys
import os
import gzip
import argparse


def run_all(args):
	ALPHABET = ['A', 'C', 'G', 'T']
	
	output_dir = args['output_dir']
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	
	print('Simulating true barcodes')
	true_barcodes = get_barcodes(
		args['num_barcodes'],
		args['barcode_length'],
		args['error_type'],
		ALPHABET)
	#print('Simulating barcode abundance')
	print('Writing reads fastq.gz')
	write_reads(
		true_barcodes,
		args['num_reads'],
		args['barcode_length'],
		ALPHABET,
		args['umi_length'],
		output_dir)
	print('Saving true barcodes txt')
	write_barcodes(
		true_barcodes,
		output_dir)
	print('Done')
	

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
	
	abundances /= np.sum(abundances)
	return abundances

def add_single_error(true_barcode, BARCODE_LENGTH, ALPHABET, mutations = None):
	#returns true_barcode sequence with one SNP, INSERTION, or DELETION
	mutated_barcode = list(true_barcode)
	mutation_pos = np.random.choice(BARCODE_LENGTH)
	if(mutations == None):
		mutations = ['mismatch', 'insertion', 'deletion']
	#mutations = ['mismatch', 'deletion', 'insertion', 'none']
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

def add_multiple_errors(true_barcode, BARCODE_LENGTH, ALPHABET, RATE=2):
	"""
	Adds a poisson number of errors with constant poisson 
	"""
	seq = true_barcode
	num_errors = np.random.poisson(lam=RATE)
	for i in range(num_errors):
		seq = add_single_error(seq, BARCODE_LENGTH, ALPHABET)
	return seq, num_errors

def write_reads(barcodes, num_reads, BARCODE_LENGTH, ALPHABET, UMI_LENGTH, output_dir):
	barcodes_file = '%s/barcodes.fastq.gz' % output_dir
	writer = gzip.open(barcodes_file, 'wb')
	for i in range(num_reads):
		true_barcode = np.random.choice(barcodes)
		
		mutated_barcode, num_errors = add_multiple_errors(true_barcode, BARCODE_LENGTH, ALPHABET)
		umi = ''.join(np.random.choice(ALPHABET, UMI_LENGTH))
		seq = mutated_barcode + umi
		qual = ''.join(['I' for i in seq])
		
		read_name = 'ReadNum:%i_NumErr:%i_TrueBarcode:%s' % (i, num_errors, true_barcode)
		fq_line = '%s\n%s\n%s\n%s\n' % (read_name, seq, '+', qual)
		
		writer.write(fq_line.encode('utf-8'))
	writer.close()

def write_barcodes(barcodes, output_dir):
	true_barcodes_file = '%s/true_barcodes.txt' % output_dir
	with open(true_barcodes_file, 'w') as writer:
		writer.write('\n'.join(barcodes))




def get_args():
	parser = argparse.ArgumentParser(
		description = 'This script simulates error-prone barcodes like Dropseq')
		
	parser.add_argument('--barcode_length',
		help = 'Length of the cell barcode',
		type=int,
		default=12)
	parser.add_argument('--umi_length',
		help = 'Length of the molecular barcode',
		type=int,
		default=8)
	parser.add_argument('--num_barcodes',
		help = 'Number of cell barcodes to simulate',
		type=int,
		default=500)
	parser.add_argument('--num_reads',
		help = 'Number of reads to simulate',
		type=int,
		default=1000000)
	parser.add_argument('--poisson_errors',
		help = 'Median number of errors per read',
		type=int,
		default=2)
	parser.add_argument('--error_type',
		help = 'Type of error',
		type=str,
		default='any')
	parser.add_argument('--output_dir',
		help = 'Directory to write output files to',
		type=str,
		required=True)
	return vars(parser.parse_args())
	

if __name__ == '__main__':
	args = get_args()
	run_all(args)





