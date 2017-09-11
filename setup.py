"""

"""

import shlex
import sys
import os
import json

from setuptools import setup, find_packages

params = {}
args = shlex.split(' '.join(sys.argv))
print(args)

if('--kallisto' in args):
	index = args.index('--kallisto')
	assert(index < len(args)), \
		'--kallisto option requires a path'
	
	kallisto_path = args[index + 1]
	assert os.path.exists(kallisto_path), \
		'kallisto path is invalid.\n%s' % kallisto_path
	params['kallisto'] = kallisto_path
	
	sys.argv.remove('--kallisto')
	sys.argv.remove(kallisto_path)
	
else:
	params['kallisto'] = 'kallisto'

"""
prepare params.json
"""

with open('./sircel/params.json', 'w') as writer:
	writer.write(json.dumps(params, indent = 3))

setup(name='sircel',
      version = '0.1.1',
      packages = find_packages(),
      install_requires = [
			'numpy',
			'scipy',
			'python-Levenshtein',
			'matplotlib'],
      package_data = {'': ['params.json']},
      include_package_data = True,
		entry_points = {
			'console_scripts' : [
				'sircel = sircel.__main__:main']},
		
	   description = 'Identify and error correct barcodes for single-cell genomics',
	   url = 'https://github.com/pachterlab/Sircel',
	   author = 'Akshay Tambe',
	   author_email = 'akshay.tambe@berkeley.edu',
	   license = 'MIT',)

