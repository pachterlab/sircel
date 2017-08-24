from sircel.Sircel_master import get_args, run_all

def main():	
	args = get_args()
	output_files = run_all(args)

if __name__ == "__main__":
	main()