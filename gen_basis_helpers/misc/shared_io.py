
import os
import json
import pathlib

def dumpObjWithToDictToJson(inpObj, outFile):
	""" Write an object to outFile (json format) if inpObj has a .toDict() method
	
	Args:
		inpObj: The object we want to dump to a json file
		outFile: (str) Path to the output file
 
	"""
	outDir = os.path.split(outFile)[0]
	pathlib.Path(outDir).mkdir(parents=True, exist_ok=True)
	with open(outFile,"wt") as f:
		f.write(json.dumps(inpObj.toDict()))


def readObjWithFromDictFromJsonFile(inpCls, inpFile):
	""" Initializes instance of inpCls directly from JSON file (inpFile) if inpCls has a fromDict() class method (i.e. alternative initializer)
	
	Args:
		inpCls: (Class, NOT INSTANCE, with a fromDict() initializer)
		inpFile: (str) Path to the .json file

	Returns
		outInstance: Instance of inpCls with attributes taken from inpFile
 
	"""
	with open(inpFile,"rt") as f:
		currDict = json.loads(f.read())
		outObj = inpCls.fromDict(currDict)
	return outObj

