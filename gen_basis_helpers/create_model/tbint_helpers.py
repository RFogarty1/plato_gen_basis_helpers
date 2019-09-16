
import os
import subprocess


from ..shared import ch_dir as chDir


def runTbintFromCommStr(filePath:"Should really be a folder", optStr=None, startDir=None, elements:list=None):

	if startDir is None:
		startDir = os.getcwd()
	if optStr is None:
		optStr = ""

	basePath = os.path.splitext(filePath)[0]

	if os.path.isdir(basePath):
		workFolder = basePath
		fileName = None
	else:
		workFolder, fileName = os.path.split(basePath)


	if elements is None: #For backwards compat. with stupid way I orginally designed this
		if fileName is None:
			raise ValueError("{} is a folder but elements is not set".format(fileName))
		elementStr = fileName
	else:
		elementStr = " ".join([x for x in elements])

	runComm = "tbint " + optStr.strip() + " " + elementStr

	with chDir.ChDir(startDir, workFolder):
		subprocess.check_call(runComm,shell=True)


def removeTbintFilesFromFolder(folder):
	fileExts = [".adt",".bdt"]
	relevantFiles = [os.path.join(folder,x) for x in os.listdir(folder) if _isTbintOutFile(x)]
	[os.remove(x) for x in relevantFiles]



def _isTbintOutFile(filePath):
	fileExts = [".adt",".bdt"]
	for ext in fileExts:
		if filePath.endswith(ext):
			return True
	return False

