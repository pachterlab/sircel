This set of scripts uses simulations to rigorously benchmarks sircel's ability to identify barcodes from error containing short reads. There are three scripts here that must be run sequentially. Simulation results can be inspected with a python notebook (provided in the notebooks directory of this repo)

	Simulate_multiple_datasets.py
		Prepares a large number of simulated barcodes.fastq.gz files
		User-defined parameters are: 
			The number of barcodes to simulate
			The number of reads to simulate
			The length of barcode sequences
			The length of UMI sequences
		The script then systematically prepares datasets for combinatios of each of the following parameters:
			Barcode abundance distribution (uniform, normal, exponential)
			Barcode error type (insertion, deletion, mismatch, any)
			Poisson error rate (1,2,3,4,5)
			
	Evaluate_simulations.py
	Run_simulations.py
	

	