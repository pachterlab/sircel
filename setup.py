"""

"""

import shlex
import sys
import os

from setuptools import setup

params = {}
args = shlex.split(' '.join(sys.argv))

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
	params['kallisto'] = None

if('--osx' in args):
	params['zcat'] = 'gzcat'	#zcat function is broken on mac
	sys.argv.remove('--osx')
else:
	params['zcat'] = 'zcat'




setup(name='sircel',
      version='0.1',
      description='Identify and error correct barcodes for single-cell genomics',
      url='https://github.com/pachterlab/Sircel',
      author='Akshay Tambe',
      author_email='akshay.tambe@berkeley.edu',
      license='MIT',
      packages=['sircel'],
		py_modules=['numpy', 'scipy', 'sklearn', 'redis'])



"""
prepare params.json
"""

import json 
params['sircel'] = './sircel/sircel_master.py'


with open('./sircel/params.json', 'w') as writer:
	writer.write(json.dumps(params, indent = 3))









	