"""
Akshay Tambe
Pachter and Doudna groups
"""

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from collections import Counter

def plot_kmer_subsamp_pearson(output_dir, counts_corr_coefs, num_reads):
	fig, ax = plt.subplots(
		nrows = 1, 
		ncols = 1,
		figsize = (4,4))
	
	ax.plot(num_reads, counts_corr_coefs)
	ax.set_xlabel('Number of reads indexed')
	ax.set_ylabel('Pearson R')
	ax.set_title(\
		'Correlation between relative kmer counts' + \
		'\nas number of reads are incrementally increased')
	ax.grid()
	ax.set_ylim([0,1])
	plt.tight_layout()
	
	fname = '%s/plots/indexed_reads_correlation.png' % output_dir
	fig.savefig(fname, dpi = 300)
	return fname

def plot_path_threshold(params):
	(output_dir, 
		path_weights,
		grad,
		second_grad,
		lmax,
		threshold,
		LOCAL_WINDOW_LEN) = params
	fig, ax = plt.subplots(
		nrows = 1, 
		ncols = 3,
		figsize = (12,4))
	x = range(len(path_weights))
	ax[0].step(x, path_weights, label='Path weights')
	ax[0].set_yscale('log')
	
	x = [LOCAL_WINDOW_LEN / 2 + i for i in range(len(grad))]
	ax[1].plot(x, grad, color='b', label='First derivative (smoothed)')
	x = [LOCAL_WINDOW_LEN + i for i in range(len(second_grad))]
	ax[2].plot(x, second_grad, color = 'b', label = 'Second derivative (smoothed)')
	
	for a in ax:
		for m in lmax:
			a.axvline(m, color = 'gray', alpha = 0.25)
		a.axvline(threshold, color='k', label='Threshold')
		a.legend(loc = 1, fontsize = 6)
		a.grid()
		a.set_xscale('log')
		a.set_xlim([0, len(path_weights)])
		a.set_xlabel('Path (sorted)')
	
	plt.tight_layout()
	fname = '%s/plots/paths_plotted.png' % output_dir
	fig.savefig(fname, dpi = 300)
	return fname

def plot_nuc_content(output_dir, reads_nuc_content, barcodes_nuc_content):	
	fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (8,8))
	colors = 'kbry'
	
	for i, nuc in enumerate('ACGT'):
		y = barcodes_nuc_content[nuc]
		x = list(range(len(y)))
		ax[0].plot(x, y, color = colors[i], label = nuc)
		
		y = reads_nuc_content[nuc]
		x = list(range(len(y)))
		ax[1].plot(x, y, color = colors[i], label = nuc)
	
	ax[0].set_xlabel( \
		'Position on barcode read\n[after merging for 10x genomcs]')
	ax[1].set_xlabel('Position on RNAseq read')
	
	for a in ax:
		a.set_ylabel('Frequency')
		a.legend(fontsize = 8)
		a.grid()
		a.set_title('Nucleotide frequency by position')
	plt.tight_layout()
	
	fname = '%s/plots/nucleotide_freq_vs_pos.png' % output_dir
	fig.savefig(fname, dpi = 300)
	return fname

def plot_cell_distance_hmap(output_dir, distances, bc_len):
	"""
	distances is a
	list of tuples: (cell_name (str), lev_dist (counter), ham_dist (counter) )
	
	scatter plot for each cell: 
		average lev dist (x), 
		average ham dist(y)
		ebars
	"""
	
	fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (4,4))
	
	num_cells = len(distances)
	ham_medians = []
	lev_medians = []
	
	for i, (cell_name, lev_dist, ham_dist) in enumerate(distances):
		ham_mean = np.mean(ham_dist)
		ham_std = np.std(ham_dist)
		lev_mean = np.mean(lev_dist)
		lev_std = np.std(lev_dist)
		
		ax.scatter(ham_mean, lev_mean, alpha = 0.1)
		#ax.errorbar(ham_mean, lev_mean, xerr=ham_std, yerr = lev_std, alpha = 0.1, color = 'b')
	
	ax.plot([0,15], [0,15], color = 'grey', ls='-')
	ax.set_xlim([0, bc_len])
	ax.set_ylim([0, bc_len])
	ax.grid()
	
	ax.set_xlabel('Average Hamming distance to consensus')
	ax.set_ylabel('Average Levenshtein distance to consensus')
	
	plt.tight_layout()
	fname = '%s/plots/per_cell_edit_distance.png' % output_dir
	fig.savefig(fname, dpi = 300)
	return fname
		

		
		
		
		
		
		
		
		
		
		