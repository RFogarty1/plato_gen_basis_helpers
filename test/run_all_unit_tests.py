#!/usr/bin/python3

import os
import sys
import unittest

def main():
	startDir = os.path.split( os.path.abspath(os.getcwd(),) )[0]

	print("startDir = {}".format(startDir))
	#Find all unit-tests and run them
	loader = unittest.TestLoader()
	suite = loader.discover(startDir, pattern='*utest*.py')
	runner = unittest.TextTestRunner(verbosity=2)
	result = runner.run(suite)

	sys.exit( int(not result.wasSuccessful()) ) #return 0 on success, else 1

if __name__ == '__main__':
	main()
