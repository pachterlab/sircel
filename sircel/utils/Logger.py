
import sys

class Logger(object):
	"""
	Copypasta from stackoverflow forums.
	All credit to users Amith Koujalgi and Eric Leschinsky
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
	
	def close(self):
		self.log.close()

def start(filename):
    """Start logger, appending print output to given filename"""
    sys.stdout = Logger(filename)

def stop():
    """Stop logger and return print functionality to normal"""
    sys.stdout.log.close()
    sys.stdout = sys.stdout.terminal