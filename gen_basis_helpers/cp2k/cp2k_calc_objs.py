
import re
import os
import pathlib
import types

import plato_pylib.parseOther.parse_cp2k_files as parseCP2K

from ..shared import method_objs as methodObjs
from . import cp2k_file_helpers as pyCP2KHelpers
from . import parse_md_files as parseMdHelp
from . import parse_neb_files as parseNebHelp

#NOTE: Loads of descriptors are added below (at the bottom of the file)
class CP2KCalcObj(methodObjs.CalcMethod):

	def __init__(self, pycp2kObj, basePath=None, saveRestartFile=True, md=False, postWriteHooks=None, runType=None):
		""" Initializer
		
		Args:
			postWriteHooks: (iter of f(instance)) These are called by writeFile after the main file is written. Original use is for copying restart files to a target folder
				 
		"""
		self.cp2kObj = pycp2kObj
		if basePath is None:
			self.basePath = os.path.abspath( os.path.join( os.getcwd(), "cp2k_file" ) )
		else:
			self.basePath = os.path.abspath( os.path.splitext(basePath)[0] )
		self.saveRestartFile = saveRestartFile
		self.md = md
		self.postWriteHooks = list() if postWriteHooks is None else list(postWriteHooks)
		self.runType = runType
	
	def writeFile(self):
		self.cp2kObj.project_name = os.path.split(self.basePath)[1]
		self.cp2kObj.working_directory = os.path.split(self.basePath)[0]
		pathlib.Path(self.cp2kObj.working_directory).mkdir(parents=True,exist_ok=True)
		self.cp2kObj.write_input_file()
		for f in self.postWriteHooks:
			f(self)	
	
	@property
	def outFilePath(self):
		return self.basePath + ".cpout"

	@property
	def outGeomPath(self):
		return self.basePath + "-pos-1.xyz"
	
	@property
	def nCores(self):
		return 1 #I only implement the serial version. Error should occur if trying to set this
	
	@property
	def runComm(self):
		baseName = os.path.split(self.basePath)[1]
		inpFolder= os.path.abspath(os.path.split(self.basePath)[0] )
		inpFName = baseName + ".inp"
		outFName = os.path.split(self.outFilePath)[1]
		commFmt = "cd {};cp2k.sopt {}>{}"
		if self.saveRestartFile is not True:
			commFmt += ";rm {}-RESTART.kp*".format(baseName)
		return commFmt.format(inpFolder, inpFName, outFName)
	
	@property
	def parsedFile(self):
		#Figure out runType
		if self.md:
			runType = "md"
		else:
			runType = self.runType

		#Parse accordingly
		if runType is None:
			parsedDict = parseCP2K.parseCpout(self.outFilePath)
			outObj = types.SimpleNamespace(**parsedDict)
			try:
				outCartCoords = self._getFinalCartCoordsFromOpt()
				outObj.unitCell.cartCoords = outCartCoords
			except FileNotFoundError:
				pass
			outObj.unitCell.convAngToBohr()
		elif runType.lower()=="md":
			outObj = self._parseMdStandard()
		elif runType.lower() == "band":
			outObj = self._parsedNudgedBandStandard()
		else:
			raise ValueError("{} is an invalid runType".format(runType))

		return outObj

	def _getFinalCartCoordsFromOpt(self):
		parsedXyzDict = parseCP2K.parseXyzFromGeomOpt(self.outGeomPath)
		return parsedXyzDict["all_geoms"][-1].cartCoords

	def _parseMdStandard(self):
		#This case should ONLY trigger if multiple runs wernt required
		workFolder = os.path.split(self.outFilePath)[0]
		if os.path.exists(self.outFilePath) and not(os.path.exists( os.path.join(workFolder,"run_1") )):
			cpoutPaths = [self.outFilePath]
			xyzPaths = [self.outGeomPath]
			tKindPaths = [os.path.join(workFolder,x) for x in os.listdir(workFolder) if x.endswith(".temp")]
			if len(tKindPaths) != 1:
				tKindPaths = [None]
			finalRunFolder = workFolder #So we can find restart files etc
		else:
			cpoutPaths, xyzPaths, tKindPaths = list(), list(), list()
			runDirs = [x for x in os.listdir(workFolder) if os.path.isdir( os.path.join(workFolder,x) )]
			runDirs = sorted( [x for x in runDirs if self._checkStringMatchesRunFormat(x)] )
			#NEXT: We want to check for one cpout and one xyz per path
			for runDir in runDirs:
				currDir = os.path.join(workFolder,runDir)
				currXyz = [os.path.join(currDir,x) for x in os.listdir(currDir) if x.endswith(".xyz")]
				currCpout = [os.path.join(currDir,x) for x in os.listdir(currDir) if x.endswith(".cpout")]
				currTKind = [os.path.join(currDir,x) for x in os.listdir(currDir) if x.endswith(".temp")]
				assert len(currXyz)==1
				assert len(currCpout)==1
				cpoutPaths += currCpout
				xyzPaths += currXyz
				if len(currTKind)==1:
					tKindPaths += currTKind
				else:
					tKindPaths += [None]
				finalRunFolder = os.path.join(workFolder,runDirs[-1])

		assert len(cpoutPaths) == len(xyzPaths)

		parsedDict = parseMdHelp.parseMdInfoFromMultipleCpoutAndXyzPaths(cpoutPaths,xyzPaths, tempKindPaths=tKindPaths)
		return types.SimpleNamespace(finalRunFolder=finalRunFolder,**parsedDict)

	def _parsedNudgedBandStandard(self):
		parsedDict = parseNebHelp.parseNudgedBandCalcStandard(self.outFilePath, convAngToBohr=True)
		return types.SimpleNamespace(**parsedDict)

	def _checkStringMatchesRunFormat(self, inpStr):
		pattern = "run_[0-9]+"
		matchObj = re.match(pattern, inpStr)
		if matchObj is None:
			return False
		if inpStr != inpStr[slice(*matchObj.span())]:
			return False
		return True

#Optional descriptors that can be added
def addInpPathDescriptorToCP2KCalcObjCLASS(inpCls):
	attrName = "inpPath"
	setattr(inpCls, attrName, InpPathDescriptorCP2K(attrName))

def addQuickStepEpsDefDescriptorToCP2KCalcObjCLASS(inpCls, attrName="qs_eps_def"):
	setattr(inpCls, attrName, QuickStepEpsDefaultCP2K(attrName))

def addRelGridDescriptorToCP2KCalcObjCLASS(inpCls, attrName="relGridCutoff"):
	units = "eV"
	setattr(inpCls, attrName, RelGridDescriptorCP2K(attrName,units))

def addAbsGridDescriptorToCP2KCalcObjCLASS(inpCls, attrName="absGridCutoff"):
	units="eV"
	setattr(inpCls, attrName, AbsGridDescriptorCP2K(attrName,units))

def addMaxScfDescriptorToCP2KCalcObjCLASS(inpCls, attrName="maxScf"):
	setattr(inpCls, attrName, MaxScfDescriptorCP2K(attrName) )

def addAddedMOsDescriptorToCP2KCalcObjCLASS(inpCls, attrName="addedMOs"):
	setattr(inpCls, attrName, AddedMOsDescriptorCP2K(attrName))


#Defining an inpPath property for CP2KCalcObj
#Note - it is VERY HARD to include a docsting on these (id need to return a special class with a repr well defined instead of a str from get, the docstirng will be accessed on the return class
class InpPathDescriptorCP2K(object): 

	def __init__(self, attrName):
		self.name = attrName

	def __get__(self, instance, owner):
		return instance.basePath + ".inp"
  
	def __set__(self, instance, value):
		raise NotImplementedError("Cannot set attribute {}".format(self.name))

#Defining quickstep_eps_default parameter for CP2KCalcObj
class QuickStepEpsDefaultCP2K(object):

	def __init__(self,attrName):
		self.name = attrName

	def __get__(self, instance, owner):
		return instance.cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.QS.Eps_default

	def __set__(self, instance, value):
		assert( len(instance.cp2kObj.CP2K_INPUT.FORCE_EVAL_list) == 1 )
		instance.cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.QS.Eps_default = value

#Defining relative grid spacing in the CP2K object
class RelGridDescriptorCP2K():

	def __init__(self, attrName, units):
		self.name = attrName
		self.units = units.lower()

	def __get__(self, instance, owner):
		strVal = instance.cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.MGRID.Rel_cutoff
		self._checkUnitsCorrect(strVal)
		floatVal = float(strVal.strip().split()[-1])
		return floatVal

	#We can probably ALWAYS insist the units are in eV, since that makes it easier to our acceptor. If user wants difference
	#units for the interface we can convert them internally (will implement if ever actually needed)
	def _checkUnitsCorrect(self,strVal):
		assert self.units=="ev", "For now only ev units are supported for RelGridDescriptor"
		assert "ev" in strVal.lower(), "Units eV need to be present in relative cutoff str:{}".format(strVal) 
			
	#Note this can only set in eV for now. Easiest modificatoin is probably to always convert units if others requested
	def __set__(self, instance, value):
		assert len(instance.cp2kObj.CP2K_INPUT.FORCE_EVAL_list) == 1
		modDict = {"gridCutRel":value}
		pyCP2KHelpers.modCp2kObjBasedOnDict(instance.cp2kObj, modDict)

class AbsGridDescriptorCP2K():

	def __init__(self, attrName, units):
		self.name = attrName
		self.units = units.lower()

	def __get__(self, instance, owner):
		strVal = instance.cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.MGRID.Cutoff
		self._checkUnitsCorrect(strVal)
		floatVal = float(strVal.strip().split()[-1])
		return floatVal

	def __set__(self, instance, value):
		assert len(instance.cp2kObj.CP2K_INPUT.FORCE_EVAL_list) == 1 
		modDict = {"gridCutAbs":value}
		pyCP2KHelpers.modCp2kObjBasedOnDict(instance.cp2kObj, modDict)

	def _checkUnitsCorrect(self,strVal):
		assert self.units=="ev", "For now only ev units are supported for RelGridDescriptor"
		assert "ev" in strVal.lower(), "Units eV need to be present in relative cutoff str:{}".format(strVal) 


class MaxScfDescriptorCP2K():

	def __init__(self, attrName):
		self.name = attrName

	def __get__(self, instance, owner):
		strVal = instance.cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.SCF.Max_scf
		outVal = int( strVal.strip().split()[-1] )
		return outVal

	def __set__(self, instance, value):
		assert len(instance.cp2kObj.CP2K_INPUT.FORCE_EVAL_list) == 1
		modDict = {"maxscf":value}
		pyCP2KHelpers.modCp2kObjBasedOnDict(instance.cp2kObj, modDict)



class AddedMOsDescriptorCP2K():

	def __init__(self, attrName):
		self.name = attrName

	def __get__(self, instance, owner):
		return int( instance.cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.SCF.Added_mos )

	def __set__(self, instance, value):
		assert len(instance.cp2kObj.CP2K_INPUT.FORCE_EVAL_list) == 1
		modDict = {"addedMOs".lower():value}
		pyCP2KHelpers.modCp2kObjBasedOnDict(instance.cp2kObj, modDict)



#Add all the defined descriptors. At first the plan was only to add these when they were needed,
#but testing just got too confusing and annoying that way
addInpPathDescriptorToCP2KCalcObjCLASS(CP2KCalcObj)
addQuickStepEpsDefDescriptorToCP2KCalcObjCLASS(CP2KCalcObj)
addRelGridDescriptorToCP2KCalcObjCLASS(CP2KCalcObj)
addAbsGridDescriptorToCP2KCalcObjCLASS(CP2KCalcObj)
addMaxScfDescriptorToCP2KCalcObjCLASS(CP2KCalcObj)
addAddedMOsDescriptorToCP2KCalcObjCLASS(CP2KCalcObj)






