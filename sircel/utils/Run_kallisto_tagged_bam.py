"""
Annotated bam to split_files
	Barcodes called by Macosko et al are in "XM" tag of bam
		Tags:	XC:Z:CGACAAATTGTT	XF:Z:INTERGENIC	PG:Z:STAR	RG:Z:A	NH:i:2	NM:i:3	XM:Z:GCAGGGCA	UQ:i:37	AS:i:51
		Cell: XM:Z:GCAGGGCA
	1. Read bam and get the set of cells and 
	2. 

"""
import subprocess
from sircel import Sircel_master
import sys
import os
from collections import Counter
import io

def run_all(args):
	output_files = {}
	batch = '%s/reads_split/batch.txt'
	output_files['split'] = {}
	output_files['split']['batch'] = batch
	
	print('Getting cell barcodes')
	cell_bcs = get_consensus_cells(args['bam'])
	for bc in cell_bcs:
		fq, umi = write_cell_to_fq(bc, bam, outf)
		batch.write('%s\t%s\t%s\n' % (bc, umi, fq))
		output_files['split'][bc] = \
			{'barcodes' : None, 'umis' : umi, 'reads' : fq}
	#run kallisto
	
	#print(args['kallisto_idx'])
	print('Running kallisto')
	kallisto_dir = '%s/kallisto_outputs' % args['output_dir']
	if not os.path.exists(kallisto_dir):
		os.makedirs(kallisto_dir)
	output_files['kallisto'] = Sircel_master.run_kallisto(
		args,
		kallisto,
		kallisto_dir,
		output_files)
	print('Getting transcript compatibility counts')
	output_files['tcc'] = Sircel_master.write_transcript_compatability_counts(
		args,
		output_files,
		kallisto_dir)
	


def write_cell_to_fq(bc, bam, outf):
	out_fq = '%s/reads_split/cell_%s.fastq.gz' % (output_dir, bc)
	out_umi = '%s/reads_split/cell_%s.umi.txt' % (output_dir, bc)
	fq_writer = gzip.open(out_fq, 'wb')
	umi_writer = open(out_umi, 'w')
	
	sam_iter = read_bam(bam)
	for entries in sam_iter:
		(name, seq, phred, cell_tag, umi_tag) = entries
		if cell_tag == bc:
			fq_writer.write(b'%s\n%s\n+\n%s\n' % (name, seq, phred))
			umi_writer.write('%s\n' % umi_tag)
	fq_writer.close()
	umi_writer.close()
	return out_fq, out_umi
	

def get_consensus_cells(bam):
	cell_bcs = {}
	
	sam_iter = read_bam(bam)
	for entries in sam_iter:
		(name, seq, phred, cell_tag, umi_tag) = entries
		cell_bcs[cell_tag] = True
	print('Identified %i barcodes' % len(cell_bcs))
	return cell_bcs
		
def get_tags(row):
	tags_lst = row[11:]
	tags = {}
	for tag in tags_lst:
		entries = tag.split(':')
		tags[entries[0]] = '_'.join(entries[2:])
	return tags
		
def read_bam(bam):
	samtools_cmd = ['samtools', 'view', bam]
	samtools_view = subprocess.Popen(
		samtools_cmd,
		stdout = subprocess.PIPE)
	for sam_line in io.TextIOWrapper(samtools_view.stdout, encoding='utf-8'):		
		if sam_line[0] != '@':	#sam header. shouldn't show up normally anyway
			row = sam_line.rstrip().split()
			name = row[0]
			seq = row[8]
			phred = row[9]
			cell_tag = get_tags(row)['XC']
			umi_tag = get_tags(row)['XM']			
			yield (name, seq, phred, cell_tag, umi_tag)

if __name__ == "__main__":
	args = {'threads' : '32'}
	args['bam'] = sys.argv[1]
	args['output_dir'] = sys.argv[2]
	args['kallisto_idx'] = sys.argv[3]
	
	run_all(args)
	
	
	
	
	
	
	
	
	
	