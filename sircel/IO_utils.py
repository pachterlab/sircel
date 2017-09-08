"""
Akshay Tambe
Pachter and Doudna groups
"""
import subprocess
import sys
import tempfile
import numpy as np
import redis
import itertools
import shutil
import gzip as gz

np.random.seed(0)

def get_nuc_content(unzipped, num_reads):
	nuc_content = {}
	increment = 1.0 / num_reads
	seq_len = 0
	
	with open(unzipped) as in_file:
		for lines in grouper(in_file, 4):
			seq = lines[1].strip().upper()
			
			#some datasets have variable read lengths
			if len(seq) > seq_len:
				seq_len = len(seq)
				for nt, freqs in nuc_content.items():
					while len(freqs) < seq_len:
						freqs.append(0)
					nuc_content[nt] = freqs
			
			for i, nt in enumerate(seq):
				if(nt not in nuc_content):
					nuc_content[nt] = [0] * len(seq)
				nuc_content[nt][i] += increment
	return nuc_content

def get_kmers(sequence, k):
	"""
	Args:
		sequence (string)
		k (int)
	Returns list
		all subsequences of sequence of length ks
	"""
	kmers = []
	for i in range(len(sequence) - k + 1):
		kmers.append(sequence[i:i+k])
	return kmers

def get_cyclic_kmers(read, k, start, end, indel=True):
	"""
	Args
		read (list)
			a fastq entry as a list of lines
		k (int)
			size of kmer
		start (int)
			start site of barcode within read
		end (int)
			end site of barcode within read
	Returns list
			list of tuples (string kmer, string qual)
			the input read is circularized
	"""	
	seq = read[1]
	seq = seq[:start] + '$' + seq[start:]
	qual = read[3]
	qual = qual[:start] + '$' + qual[start:]
	
	alphabet=['A', 'C', 'G', 'T']
	if(end >= len(seq)):#this doesn't typically happen for drop seq
		end = len(seq)
		seq += np.random.choice(alphabet)#add something random for insertions 
	
	#cyclicized, correct length
	cyclic =				lambda s: s[start:end] + s[start:start + k - 1]
	#cyclicized, one nucleotide truncated
	cyclic_truncate =	lambda s: s[start:end - 1] + s[start:start + k - 1]
	#cyclicized, one nucleotide extended
	cyclic_extend =	lambda s: s[start:end + 1] + s[start:start + k - 1]
	
	kmers = 	get_kmers(cyclic(seq), k)
	quals = 	get_kmers(cyclic(qual), k) 
	if(indel == True):
		kmers += get_kmers(cyclic_truncate(seq), k) + \
			get_kmers(cyclic_extend(seq), k)
		quals += get_kmers(cyclic_truncate(qual), k) + \
			get_kmers(cyclic_extend(qual), k)
	
	ret = [tup for tup in zip(kmers, quals)]
	return ret

def unzip(gzipped_lst):
	"""
	Args
		gzipped (string)
	Returns
		temp file path (string), offsets (list), nuc_content (counter)
	"""
	out_file = tempfile.NamedTemporaryFile(delete=False)
	offsets = []
	offset = 0
	
	for gzipped in gzipped_lst:
		if not gzipped.endswith('.gz'):
			raise TypeError('File does not appear to be gzipped: %s' % gzipped)
		with gz.open(gzipped) as in_file:
			for lines in grouper(in_file, 4):
				
				lines = b''.join(lines)
				offset += len(lines)
				offsets.append(offset)
				out_file.write(lines)
				
	return out_file.name, offsets

def get_read_chunks(barcodes, offsets, random = False, BUFFER_SIZE = 10000):
	"""
	Args: 
		barcodes (unzipped fq file)
		reads (unzipped fq file)
		lines (list)
	yields
		data_buffer (list of tuples)
	
	simultaneously reads two fastq files (line-by-line) and stores data in a list
	yields this list (chunks of the datasets)
	"""
	barcodes_iter = read_fastq(barcodes, offsets, random = random)
	
	break_outer = False
	while(not break_outer):
		data_buffer = []#list of tuples(barcodes_data, barcodes_offset, reads_offset)
		while(True):
			try:
				barcodes_data, barcodes_offset = next(barcodes_iter)
			except StopIteration:
				break_outer = True
				break
			
			if(len(data_buffer) % BUFFER_SIZE == 0 and len(data_buffer) > 0):
				#break_inner = True
				break
			data_buffer.append((
				barcodes_data, 
				barcodes_offset))			
		yield data_buffer

def read_fastq(fq_file, offsets, random = False):
	"""
	Args:
		fq_file (path to a fastq format file), unzipped
		offsets (list)
	Yields
		tuple (lines, offset)
			lines: list of 4 lines from fastq file, ordered by offets
			offset: character offset from beginning of the file for this fq read
	"""
	fq_reader = open(fq_file, 'rb')
	if(random):
		np.random.shuffle(offsets)
	while True:
		try:
			offset = offsets.pop()
		except IndexError:
			break
		fq_reader.seek(offset, 0)
		lines = [fq_reader.readline().strip().decode('utf-8') for i in range(4)]
		yield (lines, offset)
	fq_reader.close()

def merge_barcodefiles_10x(cells_gz, umis_gz):
	""""
	10x genomics to dropseq conversion

	3 fastq.gz files from each file
		I1: read index, UMI
		R1: cell barcode(?)
		R2: RNAseq read

	I1- contains cell barcodes
	@ST-K00126:307:HFM3NBBXX:1:1101:3772:1244 1:N:0:NTCGCCCT
	NTCGCCCT
	+
	#AAAFJ-J

	R1- contains umis
	@ST-K00126:307:HFM3NBBXX:1:1101:3772:1244 1:N:0:NTCGCCCT
	NCATTTGAGTAACCCTGATGTCATAA
	+
	#AAFFJJJJJJJJJJJJJJJFJJJJJ

	R2- contains RNAseq data
	@ST-K00126:307:HFM3NBBXX:1:1101:3772:1244 2:N:0:NTCGCCCT
	NAAGCCAGTTGTGAATCATGCACATCAGCTCCTTCTGAAATGTGTTTATGGCCTAG
	+
	#<<A<FJJJFJFJJJJJJJJJFJFJJJJJJJJJJJJJJJFJJFFJJJJJAFJJFJF
	"""
	cells, cells_offsets = unzip(cells_gz)
	cells_iter = read_fastq(cells, cells_offsets)
	umis, umis_offsets = unzip(umis_gz)
	umis_iter = read_fastq(umis, umis_offsets)
	out_file = tempfile.NamedTemporaryFile(delete=False)
	
	out_file.seek(0)
	offset = 0
	offsets = []
	writer = open(out_file.name, 'wb')
	
	while(True):
		try:
			cell, _ = next(cells_iter)
			umi, _ = next(umis_iter)
		except StopIteration:
			break
		get_prefix = lambda r: r[0].split(' ')[0]
		assert(get_prefix(umis) == get_prefix(cells)), \
			'Reads are not in order\n%s\n%s' % \
			('\n'.join(umis), '\n'.join(cells))				
		combined_seq =  cells[1] + umis[1]
		combined_phred = cells[3] + umis[3]
		
		output = [
			cells[0],
			combined_seq,
			cells[2],
			combined_phred]
		output_str = ('\n'.join(output) + '\n').encode('utf-8')
		writer.write(output_str)
		offset += len(output_str)
		offsets.append(offset)
	writer.close()
	return out_file.name, offsets

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    args = [iter(iterable)] * n
    return zip(*args)

def save_paths_text(output_dir, paths, prefix=''):
	paths_file = '%s/%s_paths.txt' % (output_dir, prefix)
	with open(paths_file, 'w') as writer:
		for tup in sorted(
				paths,
				key = lambda tup: tup[1],
				reverse = True):
			writer.write('%s\t%i\t%i\t%s\n' %
				(tup[0], tup[1], tup[2], ','.join(tup[3])))
	return paths_file

def initialize_redis_pipeline(db=0):
	redis_db = redis.StrictRedis(host="localhost", port=6379, db=db)#redis
	
	redis_db.flushall()
	redis_pipe = redis_db.pipeline()
	return redis_db, redis_pipe

def get_from_db(kmer_idx_pipe, keys):
	for key in keys:
		kmer_idx_pipe.get(key)
	pipe_out = kmer_idx_pipe.execute()	
	entries = []
	for entry in pipe_out:
		#entry is a comma separated bytestring of ints. return just the list
		if(entry != None):
			offsets = [int(i) for i in entry.decode('utf-8').split(',')[0:-1]]
			entries.append(offsets)
	return entries

def get_digital_expression(tsv_file, cells_file):
	num_cells = get_num_cells(cells)
	nonzero_ec, ec_to_index = get_nonzero_ec(tsv_file)
	
	dge = np.zeros((num_cells, len(nonzero_ec)))
	cells_iter = read_tsv_by_cell(tsv_file)
	for (i, (cell, total_counts, ec_counts)) in enumerate(cells_iter):
		for (eq_class, count) in ec_counts:
			dge[i, ec_to_idx[eq_class]] = count / total_counts
	#print('Shape of DGE: %i %i' % (np.shape(dge)))
	return dge

def get_nonzero_ec(tsv_file):
	nonzero_ec = set()
	with open(tsv_file, 'rb') as inf:
		for line in inf:
			[eq_class, cell, count] = \
				line.strip().decode('utf-8').split('\t')
			nonzero_ec.add(eq_class)
	
	ec_to_index = {}
	for (i, ec) in enumerate(nonzero_ec):
		ec_to_index[int(ec)] = i
	num_ec = i + 1
	print('Total nonzero ECs: %i ' % num_ec)
	
	return nonzero_ec, ec_to_index

def get_num_cells(cells_file):
	cells = set()
	with open(cells_file, 'rb') as inf:
		for (i, line) in enumerate(inf):
			pass
	print('Total number of cells: %i' % i)
	return i

def read_tsv_by_cell(tsv_file):
	tsv_iter = open(tsv_file, 'rb')
	ec_counts = []
	prev_cell = ''
	total_counts = 0
	while(True):
		while(True):
			try:
				tsv_line = next(tsv_iter)
			except StopIteration:
				tsv_iter.close()
				raise StopIteration
				break
			[eq_class, cell, count] = [ \
				int(i) for i in tsv_line.decode('utf-8').strip().split('\t')]
			tsv_data = (eq_class, count)
			if(len(ec_counts) == 0):
				ec_counts.append(tsv_data)
				total_counts = count
				prev_cell = cell
			elif(prev_cell == cell):
				ec_counts.append(tsv_data)
				total_counts += count
			else:
				yield cell, total_counts, ec_counts
				ec_counts = [tsv_data]
				prev_cell = cell
				total_counts = count

class Logger(object):
	"""
	Copypasta from stackoverflow forums. All credit to user Eric Leschinsky
	http://stackoverflow.com/questions/14906764/
	
	Writes all stdout to file as well as printing to stdout
	Initialize with sys.stdout = Logger()
	"""
	def __init__(self, fname):
		self.terminal = sys.stdout
		self.log = open(fname, 'w')
		self.pos = 0

	def write(self, message):
		self.pos = self.log.tell()
		self.terminal.write(message)
		self.log.write(message)

	def flush(self):
		self.log.seek(self.pos)


