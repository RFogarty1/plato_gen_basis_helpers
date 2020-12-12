
import collections
import copy
import os
import pathlib
import types

from . import file_io as fileIoHelp
from . import lammps_parsers as lammpsParsers
from ..shared import method_objs as methObjHelp

class LammpsCalcObjStandard(methObjHelp.CalcMethod):

	def __init__(self, baseFolderPath, baseFileName, dataFileOrderedDict, scriptFileOrderedDict, typeIdxToEle=None):
		""" Initializer
		
		Args:
			baseFolderPath: (str) Path to the folder to run in [doesnt need to exist currently]
			baseFileName: (str) The name of the base file, i.e. without extensions
			dataFileOrderedDict: (OrderedDict) Represents the data file, which contains things like geometry. Keys are the headers of that file while vals are the "body" text corresponding to each header
			scriptFileOrderedDict: (OrderedDict) Represents the script file, which contains all commands for lammps to run an MD simulation
			typeIdxToEle: (dict) Keys are typeIndices while values are elements these correspond to. Needed to parse geometries which include element identities
	 
		Raises:
			KeyError: writeFile will raise this if the read_data command(key) is not found in the scriptFileOrderedDict
		"""
		self.baseFolderPath = baseFolderPath
		self.baseFileName = os.path.splitext(baseFileName)[0]
		self.dataFileOrderedDict = dataFileOrderedDict
		self.scriptFileOrderedDict = scriptFileOrderedDict
		self.typeIdxToEle = typeIdxToEle

	def writeFile(self):
		pathlib.Path(self.baseFolderPath).mkdir(parents=True, exist_ok=True)
		outDataDict = copy.deepcopy(self.dataFileOrderedDict)
		outDataDict["read_data"] = os.path.split(self.dataFilePath)[0]
		fileIoHelp.writeScriptFileFromTokens(self.scriptFilePath, self.scriptFileOrderedDict)
		fileIoHelp.writeDataFileFromTokens(self.dataFilePath, self.dataFileOrderedDict)

	@property
	def nCores(self):
		raise NotImplementedError("")

	@property
	def outFilePath(self):
		raise NotImplementedError("")

	@property
	def runComm(self):
		scriptFileName = os.path.split(self.scriptFilePath)[-1]
		outComm = "lmp -in {}".format(scriptFileName)
		return outComm

	@property
	def parsedFile(self):
		outDict = dict()
		logPath = os.path.join(self.baseFolderPath,"log.lammps")
		dumpPath = os.path.join(self.baseFolderPath,"dump.lammpstrj")

		parsedLogFile = lammpsParsers.parseLammpsLogFile(logPath)
		trajObj = lammpsParsers.getTrajectoryFromLammpsDumpFile(dumpPath, timeStep=parsedLogFile["timestep"], typeIdxToEle=self.typeIdxToEle)

		#Annoyingly slow, but we need to convert all the parsed units from angstrom to bohr
		def _convertTrajStepFromAngToBohr(trajStep):
			trajStep.unitCell.convAngToBohr()
		trajObj.applyFunctToEachTrajStep( _convertTrajStepFromAngToBohr )

		outDict["md_thermo_data"] = parsedLogFile["thermo_data"]
		outDict["md_traj"] = trajObj
		return types.SimpleNamespace(**outDict)

	@property
	def scriptFilePath(self):
		return os.path.join(self.baseFolderPath, self.baseFileName) + ".in"

	@property
	def dataFilePath(self):
		return os.path.join(self.baseFolderPath,self.baseFileName) + ".data"


class ScriptFileOptionsStandard():
	
	def __init__(self, initOpts=None, setupBox=None, setupAtoms=None,
	             forceFieldOpts=None, settingsOpts=None, fixSection=None,
	             runSection=None, outputSection=None):
		""" Initializer. Order of keywords simply determines the order optDicts are written to a file (initOpts is first, outputSection is last)
		
		Args:
			*kwargs: (OrderedDict) Each keyword argument takes an OrderedDict as input. Its possible to obtain ALL functionality by just using any ONE of the kwargs (e.g. putting everything in initOpts). The range of keywords only really exist beacuase it makes it simpler for me to separate sections of the input file (they are semi-arbritrary separations regardless) and make sure they are written in the correct order
				 
		"""
		self.initOpts = collections.OrderedDict() if initOpts is None else collections.OrderedDict(initOpts)
		self.setupBox = collections.OrderedDict() if setupBox is None else collections.OrderedDict(setupBox)
		self.setupAtoms = collections.OrderedDict() if setupAtoms is None else collections.OrderedDict(setupAtoms)
		self.forceFieldOpts = collections.OrderedDict() if forceFieldOpts is None else collections.OrderedDict(forceFieldOpts)
		self.settingsSection = collections.OrderedDict() if settingsOpts is None else collections.OrderedDict(settingsOpts)
		self.fixSection = collections.OrderedDict() if fixSection is None else collections.OrderedDict(fixSection)
		self.runSection = collections.OrderedDict() if runSection is None else collections.OrderedDict(runSection)
		self.outputSection = collections.OrderedDict() if outputSection is None else collections.OrderedDict(outputSection)

	def getOutputDict(self):
		outDict = collections.OrderedDict()
		kwargs = ["initOpts", "setupBox", "setupAtoms", "forceFieldOpts", "settingsSection", "fixSection", "outputSection", "runSection"]
		for kwarg in kwargs:
			outDict.update( getattr(self,kwarg) )
		return outDict

