
import os
import pathlib
import shutil

from ..analyse_md import thermo_data as thermoDataHelp
from ..analyse_md import traj_core as trajHelp



def getMdThermoObjFromBaseDbFolderAndRecord(baseDbFolder, inpRecord):
	thermoPath = getMdFilePathFromRecord(baseDbFolder, inpRecord, "md_thermo_filename")
	return thermoDataHelp.readThermoDataFromFile(thermoPath)

def getMdTrajInMemFromBaseDbFolderAndRecord(baseDbFolder, inpRecord):
	trajPath = getMdFilePathFromRecord(baseDbFolder, inpRecord, "md_traj_filename")
	return trajHelp.readTrajObjFromFileToTrajectoryInMemory(trajPath)

def getMdFilePathFromRecord(baseDbFolder, inpRecord, key):
	""" Gets a full path to an MD-file (e.g. traj-dump, thermo-dump etc) 
	
	Args:
		baseDbFolder: (str) Path to the base database_input folder
		inpRecord: (dict) MD-database record, should have partly been generated with getOutDictForMDFromStdOutObj_simple
		key: (str) The key in the inpRecord which contains a file name (not the full path). Standard examples are "md_thermo_filename", 
		           "md_traj_filename", "wfn_filename", "restart_filename"

	Returns
		outPath: (str) Full path to the file
 
	"""
	mdFileFolder = os.path.join( baseDbFolder, inpRecord["db_ext_path"], inpRecord["md_path_ext"] )
	return os.path.join(mdFileFolder, inpRecord[key])

def getOutDictForMDFromStdOutObj_simple(startDir, stdOutObj, writeFiles=True):
	""" Get output dictionary (to dump to database) and write relevant files from parsed dict
	
	Args:
		startDbDir: (str) Path to a base folder for writing database related things. Files will be written to subpaths starting here
		stdOutObj: () Standard output object containing data for a single MD run. Restricted to the "parsedOutFile" type here, meaning data will be contained 
		writeFiles: (Bool, optional) If True then trajectory and restart files will be written. Information on paths will be contained in outDict (
 
	Returns
		outDict: (dict) contains information on file paths containg MD run results (e.g. temperatures, trajectory)
 
	"""
	dumper = MDFilesFromOutObjsFileDumperStandard(copyWfnFile=writeFiles, copyRestartFile=writeFiles)
	outDict = dumper.dumpFiles(stdOutObj, startDir)
	return outDict

class MDFilesFromOutObjsFileDumperStandard():
	""" Class to handle dumping of files associated with MD run from a normal std output obj

	"""

	def __init__(self, copyWfnFile=True, copyRestartFile=True):
		self.copyWfnFile = copyWfnFile
		self.copyRestartFile = copyRestartFile

	def dumpFiles(self, stdOutObj, startDir):
		outDict = self._getFilePathExtensionDict(stdOutObj)
		if self.copyRestartFile:
			self._copyRestartFiles(startDir, outDict)
		if self.copyWfnFile:
			self._dumpMdFiles(stdOutObj, startDir, outDict)
		return outDict

	def _dumpMdFiles(self, stdOutObj, startDir, pathExtDict):
		#Get thermo data and trajectory objects
		assert len(stdOutObj.data[0])==1
		thermoObj = stdOutObj.data[0][0].parsedFile.thermo_data
		trajObj = stdOutObj.data[0][0].parsedFile.trajectory
		outFolder = os.path.join(startDir, pathExtDict["md_path_ext"])
		outThermoPath = os.path.join(outFolder, "out_thermo.thermo")
		outTrajPath = os.path.join(outFolder, "out_traj.traj")

		pathlib.Path(outFolder).mkdir(parents=True,exist_ok=True)

		trajHelp.dumpTrajObjToFile(trajObj, outTrajPath)
		thermoDataHelp.dumpStandardThermoDataToFile(thermoObj, outThermoPath)

	#startDir will be where we dump the main *.json file
	def _copyRestartFiles(self, startDir, pathExtDict):
		runDir = pathExtDict["final_run_dir"]
		outDir = os.path.join(startDir, pathExtDict["md_path_ext"])
		pathlib.Path(outDir).mkdir(parents=True,exist_ok=True)

		if pathExtDict.get("restart_filename", None) is not None:
			currFileName = pathExtDict["restart_filename"]
			shutil.copy2( os.path.join(runDir, currFileName), os.path.join(outDir, currFileName) )

		if pathExtDict.get("wfn_filename", None) is not None:
			currFileName = pathExtDict["wfn_filename"]
			shutil.copy2( os.path.join(runDir, currFileName), os.path.join(outDir, currFileName) )

	def _getFilePathExtensionDict(self, stdOutObj):
		outDict = {}

		#Get path extensions for trajectory/thermo data
		assert len(stdOutObj.label)==1
		keys = [getattr(stdOutObj.label[0],x) for x in ["eleKey","structKey","methodKey"]]
		outDict["md_path_ext"] = os.path.join(*keys)
		outDict["md_thermo_filename"] = "out_thermo.thermo"
		outDict["md_traj_filename"] = "out_traj.traj"

		#names of restart files
		assert len(stdOutObj.data[0])==1
		runDir = stdOutObj.data[0][0].parsedFile.finalRunFolder
		outDict["final_run_dir"] = runDir

		fileNames = os.listdir(runDir)
#		if self.copyWfnFiles:
		if any([x.endswith("RESTART.wfn") for x in fileNames]):
			wfnFilenames = [x for x in fileNames if x.endswith("RESTART.wfn")]
			assert len(wfnFilenames)==1
			outDict["wfn_filename"] = wfnFilenames[0]

		if any([x.endswith("-1.restart") for x in fileNames]):
			restartFilenames = [x for x in fileNames if x.endswith("-1.restart")]
			assert len(restartFilenames)==1
			outDict["restart_filename"] = restartFilenames[0]

		return outDict



