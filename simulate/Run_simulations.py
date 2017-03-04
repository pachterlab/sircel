


import sys
import os
import subprocess

summary_file = sys.argv[1]
split_reads_master = sys.argv[2]

data = []
with open (summary_file, 'r') as inf:
	header = inf.readline().strip().split('\t')
	for line in inf:
		row = line.strip().split('\t')
		data.append(row)

entry_to_col = {key:col for (col, key) in enumerate(header)}
get_col = lambda key: entry_to_col[key]
for entry in data:
	simulation_dir = entry[get_col('output_dir')]
	bc_file = '%s/barcodes.fastq.gz' % simulation_dir
	output_dir = '%s/barcodes_split' % simulation_dir
	
	cmd = ['python3', split_reads_master, 
		'--dropseq',
		'--reads', bc_file,
		'--barcodes', bc_file,
		'--output_dir', output_dir,
		'--threads', '16',
		'--dropseq']
	
	
	split_reads = subprocess.Popen(cmd)
	_ = split_reads.communicate()[0]
	

print('Done')