import os
import pathlib
import json

""" Code to simplify dumping a load of dictionaries to a json file """

def dumpDictsToFilePath(dicts, filePath):
	""" Dumps dicts into a file (NOTE: It WRITES rather than appends). This is a convenience function which deals with creating the path if its not present and hanldes opening/closing of file
	
	Args:
		dicts: (iter of dicts) 
		filePath: (str, path) Location to write to

	"""
	outDir = os.path.split(filePath)[0]
	pathlib.Path(outDir).mkdir(parents=True, exist_ok=True)
	with open(filePath,'w') as fp:
		json.dump(dicts,fp)


