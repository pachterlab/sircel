Readme.md

	Akshay Tambe
	Pachter and Doudna groups
	UC Berkeley


Summary

sircel (pronounced "circle")  separates reads in a fastq file based on barcode sequences that occur at known positions of reads. This is an essential first step in analyzing single-cell genomics data from experiments such as Drop-Seq. Barcode sequences often contain deletion and/or mismatch errors that arise during barcode synthesis and sequencing, and we have designed our barcode recovery approach with these issues in mind. In addition to identifying barcodes in an unbiased manner, sircel also quantifies their abundances.

sircel is well-suited for Drop-Seq, in which barcode deletion errors are prevalent. In order to identify barcodes in a deletion-robust manner, sircel begins by counting k-mers in circularized barcodes extracted from the reads. Circularization of the barcode portion ensures that a mismatch-containing barcode will share at least one k-mer with the true barcode sequence. This is not the case when using 'linear' barcode sequences for some larger values of k.

To address insertion and deletion errors, barcodes are also circularized after either appending or deleting one sequenced base. In the case of dropseq barcoded beadas, a deletion will result in the last sequenced position of the barcode actually arising from the first base of the UMI. Conversely an insertion will result in the first base (by position) of the UMI actually arising from the last base of the barcode. As each read might potentially contain an unknown mismatch / insertion / deltion error, each read is converted to three circularized variants. Kmers are then counted from these circularized barcodes, and a De Bruijn graph is constructed where a node represents by a kmer, and edges represent two kmers that were adjacent in at least one read. Edges are weighted by how many reads connect those kmers.

In this graph, a barcode will appear as a circular path of fixed length through this graph. We define the weight of such a path to be that of the lowest-weight edge within this path. We can identify such circular paths in this graph, and such paths represent possible consensus barcode sequences. It should be noted here that reads with either 0 or 1 error can contribute to these weighted paths; reads with 2 or more errors typically do not. 

This path weight can be interpreted as the number of reads that support a given barcode at its lowest-confidence kmer. As such we identify the true consensus barcodes by selecting the paths with highest weight. This process is typically straightforward, as we have observed that the cumulative distribution of path weights contains an inflection point corresponding to true paths. This result is seen in a number of simulated and real datasets.

We assign each read in our dataset to a consensus barcode based on the number of kmers the read shares with the consensus barcode. When computing kmers for this assignment, the value of k can differ from the value used to produce the barcode De Bruijn graph. This allows us to assign reads which might include 2 or more errors.

Requirements

	python3
	numpy
	scipy
	redis
		In-memory database structure: https://redis.io/
	redis-py
		Python binding for redis server
	scikit-learn (for TCCs)
		Optional. Not needed if only splitting reads by barcode (not quantifying single-cell expression levels)
	kallisto (0.43.0 or higher)
		Optional. Not needed if only splitting reads by barcode (not quantifying single-cell expression levels)


Commandline arguments

	--dropseq			Use barcode / umi positions from Macosko et al.
	--10xgenomics		Use barcodes / umi coordinates from 10xGenomics version2 chemistry
	--barcodes			Fastq.gz file from bead barcodes (Reads 1 for dropseq, reads 1 for 10xGenomics)
	--reads				Fastq.gz file from RNA-seq / 3' sequence tags (Reads 2 for dropseq, reads 2 for 10xGenomics)
	--umis				UMIs file. (Index I1 file, only required for 10xGenomics )
	--output_dir		Directory where output files are written
	--kmer_size			Size of k-mer to use for building the barcode graph
	--depth				Search depth at each node [Default: 4]
	--breadth			Search breadth through the graph [Default: 1000]
	--threads			Number of threads to use [Default: 1]
	
Optional arguments

	--kallisto_idx		Path to pre-computed kallisto index.
								If this is provided the script will quantify single-cell expression profiles using kallisto
	--barcode_start	First position of cell barcode within read in --barcode file.
								Not needed if --dropseq or --10xgenomics arguments are provided
	--barcode_end		Last position of cell barcode within read in --barcode file.
								Not needed if --dropseq or --10xgenomics arguments are provided
	--umi_start			First position of unique molecule ID within read in --barcode file.
								Not needed if --dropseq or --10xgenomics arguments are provided
	--umi_end			Last position of unique molecule ID within read in --barcode file.
								Not needed if --dropseq or --10xgenomics arguments are provided
	
Note about barcode and UMI position indexing
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


Use example: Dropseq data from Macosko 2015

This tutorial will run this software on an example dataset from Macosko et al (SRR1873277). According to the authors' software we expect ~570 single cells from this dataset. To speed things up we only use the first 1 million reads in the file (rather than ~400 million).Note that if you are not using kallisto to quantify expression levels, omit the --kallisto_idx option. This command will run a complete singe-cell RNA-seq analysis, going from raw, un-split reads to single-cell expression profiles.

	python3 scripts/circular_kmers_pipeline/Split_reads_master.py
		--threads 32
		--dropseq
		--output_dir example
		--reads macosko/1m_reads/SRR1873277_1m_r2.fastq.gz
		--barcodes macosko/1m_reads/SRR1873277_1m_r1.fastq.gz
		--threads 32
		--kallisto_idx kallisto_idx/hgmm_kallisto.idx

The following output files will be produced in the example/ directory:

		all_paths.txt
			Tab-delimited file containing all circular paths (putative barcodes) found in the graph 
		merged_paths.txt
			Tab-delimited file containing circular paths after Hamming correction	
		paths_plotted.pdf
			Plots / histograms of circular path weights
		fits.txt
			Tab-delimited file containing Gaussian fits for plots above
		kallisto/
			Folder containing kallisto and TCC output files
			See Pachterlab scRNA-seq TCC pipeline for more:
				https://github.com/pachterlab/scRNA-Seq-TCC-prep
		reads_split/
			Folder containing fastq.gz files split by read.
			Also contains umi.txt and batch.txt files for kalisto single-cell
		run_log.txt
		run_outputs.json

The outputs from this command can be visualized with an included ipython notebook. Simply point the notebook to the appropriate run output file [codeblock 3] to re-compute the plots.


