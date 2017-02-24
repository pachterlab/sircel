Readme.txt


Author

Akshay Tambe
Pachter and Doudna groups
UC Berkeley


Summary

This program separates reads in a fastq file based on barcode sequences that occur at known positions of reads. This is an essential first step in analyzing single-cell genomics data from experiments such as Drop-seq. It should be noted that barcode sequences from Drop-seq experiments often contain deletion errors that arise during barcode synthesis [CITE] and we have designed our barcode splitting approach with this in mind. The abundances and sequences of each barcode are detected in an unbiased manner; however the position within the read where a barcode occurs must be known beforehand. 

In order to identify barcodes in a deletion-robust manner, we begin by counting k-mers within circularized barcode regions from the reads. Circularization of the barcode portion ensures that an mismatch-containing barcode will share at least one k-mer with the true barcode sequence. This is not the case when using 'linear' barcode sequences for some larger values of k.

To address insertion and deletion errors, barcodes are also circularized after either appending or deleting one sequenced base. In the case of dropseq barcoded beadas, a deletion will result in the last sequenced position of the barcode actually arising from the first base of the UMI. Conversely an insertion will result in the first base (by position) of the UMI actually arising from the last base of the barcode. As each read might potentially contain an unknown mismatch / insertion / deltion error, each read is converted to three circularized variants. Kmers are then counted from these circularized barcodes, and a De Bruijn graph is built, where a node represents by a kmer, and edges represent two kmers that were adjacent in at least one read. Edges are weighted by how many reads connect those kmers.

In this graph, a barcode will appear as a circular path of fixed length through this graph. We define the weight of such a path to be that of the lowest-weight edge within this path. We can identify such circular paths in this graph, and such paths represent possible consensus barcode sequences. It should be noted here that reads with either 0 or 1 error can contribute to these weighted paths; reads with 2 or more errors typically do not. 

This path weight can be interpreted as the number of reads that support a given barcode at its lowest-confidence kmer. As such we identify the true consensus barcodes, by selecting the paths with highest weight. This process is typically straightforward, as we have observed (in multiple datasets) that a histogram of circular path weights contains two distinct peaks. We attribute the peak with higher capacity to true barcode sequences.

We assign each read in our dataset to a consensus barcode by how many kmers the read shares with the consensus barcode. When computing kmers for this assignment, the value of k can differ from the value used to produce the barcode De Bruijn graph. This allows us to assign reads which might include 2 or more errors.






Requirements

	python3
	numpy
	scipy
	scikit-learn (for TCCs)
	kallisto (0.43.0 or higher)
		Optional. Not needed if only splitting reads by barcode (not quantifying single-cell expression levels)
	Tested with linux (Ubuntu 14.04.5) and mac OSX (10.12.1)
		Note: the zcat command on OSX is broken (see https://github.com/dalibo/pgbadger/issues/70)




Commandline arguments

	--barcodes			Fastq.gz file from read 1 (beads / UMI)
	--reads				Fastq.gz file from read 2 (RNA-seq / 3' sequence tags)
	--barcode_start	First base of the barcode within reads 1 file. For Dropseq, this is 0
	--barcode_end		Last base of the barcode within reads 1 file. For Dropseq, this is 12
	--umi_start			First base of the UMI within reads 1 file. For Dropseq, this is 12
	--umi_end			Last base of the UMI within reads 1 file. For Dropseq, this is 20
	--kmer_size			Kmer size used to build barcode De Bruijn graph.
							Higher value is more specific, lower value might work better for low-quality data
							Suggested value: 10
	
	--phred_min			Minimum phred score to use when building De Bruijn graph. Note that all reads are assigned afterwards. 
	--kmer_min			Minimum kmer size used when assigning reads to consensus barcodes
	--kmer_max			Maximum kmer size used when assigning reads to consensus barcodes
	--threads			Number of threads for multithreading
	--kallisto_idx		[Optional]. Kallisto transcriptiome index for quantifying abundances
	--output_dir		Output files are written here. Requires read / write permission

Note about barcode and UMI coordinate indexing
Barcode and UMI positions within a read follow the same indexing convension as python strings. For example, the string BARCODEUMI would have coordinates:
	barcode_start		0
	barcode_end			7
	umi_start			7
	umi_end				10
	
This follows from the following python snippet:
>>> seq = 'BARCODEUMI'
>>> seq[0:7]
'BARCODE'
>>> seq[7:10]
'UMI'



Use example

This tutorial will run this software on an example dataset from Macosko et al (SRR1873277). According to the authors' software we expect ~570 single cells from this dataset. To speed things up we only use the first 1 million reads in the file (rather than ~400 million). Use the command below to split reads note that if you are not using kallisto to quantify expression levels, omit the --kallisto_idx option.

python3 Dropseq_complete_pipeline.py 
	--barcodes macosko/1m_reads/SRR1873277_1m_r1.fastq.gz 
	--reads macosko/1m_reads/SRR1873277_1m_r2.fastq.gz 
	--barcode_start 0 
	--barcode_end 12 
	--umi_start 12
	--umi_end 20 
	--kmer_size 10
	--phred_min 30
	--kmer_min 8
	--kmer_max 12
	--threads 32
	--kallisto_idx kallisto_idx/hgmm_kallisto.idx
	--output_dir macosko/1m_reads/TEST



Outputs
	
[STUFF GOES HERE]






Run time (1 million reads from SRR1873277)
Using python 3.4.3 on Ubuntu 14.04.5 with 32 cores

	Step							Time (s)
	Find paths					948.24
		Count kmers				49.53		Runtime increases linerly with reads
		Build graph				6.48		Runtime does not depend on the number of reads
		Find paths				889.30	Runtime does not depend on the number of reads
		Misc						2.93
		
	Split reads					108.85			
		Assign reads			50.07		Runtime increases lineraly with reads
		Write output files	54.91		Runtime increases linearly with reads
		Misc						3.87
		
	kallisto						32.30		This is embarassingly fast
	Compute TCCs				16.92
	Total							1106.31					


Known issues / to to list

Proper setup.py; none of this janky stuff	
Cannot handle 10xGenomics data (because of the interleaved data format. This will be addressed soon)
Consnesus barcodes containing long tracts of homopolymer (such as a barcode with the sequence AAAAAAAAAAAA) produce junk output. Such sequences produce self-edges in the barcode De Bruijn graph, which are currently not handled well (or rather at all). Note that this problem does not affect other consensus barcodes; output is only garbage for homopolymeric barcodes. 







