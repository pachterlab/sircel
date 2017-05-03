"""
1. read summary file
2. for each simulation:
	Count the number of positive / negative barcode calls
		read true_barcodes as dict
		read top_paths as dict
			get intersection of keys, count them
		
"""
import sys
import IO_utils
import os

def run_all():
	#read summary file
	summary_file = sys.argv[1]
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
	
	writer = open('summary_processed.txt', 'w')
	writer.write('\t'.join(new_header) + '\n')
	
	
	
	for simulation_entry in summary_data:
		simulation_dir = simulation_entry[get_col('Output_dir')]
		
		(num_tp,		#number of true pos barcodes
			num_fp,	#number of false pos barcodes
			num_fn,	#number of false neg barcodes
			bc_tpr,	#true positive rate (list) of barcode assignments
			bc_fpr) = run_single_file(simulation_dir)	#false pos ...
		#save all this data
		bc_tpr_str = ','.join([str(i) for i in bc_tpr])
		bc_fpr_str = ','.join([str(i) for i in bc_fpr])
		simulation_entry += [
			num_tp, num_fp, num_fn, bc_tpr_str, bc_fpr_str]
		printer = ('\t'.join([str(i) for i in simulation_entry]))
		writer.write(printer + '\n')
	writer.close()
	print('done')

def run_single_file(simulation_dir):
	
	true_barcodes = get_barcodes(
		'%s/true_barcodes.txt' % simulation_dir)
	pred_barcodes = get_barcodes(
		'%s/threshold_paths.txt' % simulation_dir)
	
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
			bc, simulation_dir)
		barcodes_tpr.append(bc_tpr)
		barcodes_fpr.append(bc_fpr)
	
	return (num_tp,
		num_fp,
		num_fn,
		barcodes_tpr,
		barcodes_fpr)

def get_barcodes(true_bc_file):
	barcodes = set()
	try:
		inf= open(true_bc_file, 'r')
	except FileNotFoundError:
		return barcodes
	for line in inf:
		bc = line.strip().split('\t')[0]
		barcodes.add(bc)
	return barcodes

def get_fraction_correct_reads(pred_bc, simulation_output_dir):
	fq_file = '%s/reads_split/cell_%s_barcodes.fastq.gz' % \
		(simulation_output_dir, pred_bc)
	
	if not os.path.exists(fq_file):
		return (0,0)
	
	
	fq_iter = IO_utils.read_fastq_sequential(fq_file, gzip=True)
	
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



if __name__ == '__main__':
	run_all()






	
	
	
	