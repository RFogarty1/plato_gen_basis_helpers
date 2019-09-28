
import os

import numpy as np

from ..shared import calc_methods as calcMethods
from . import dos_helpers as dosHelp
from . import dos_base_objs as baseObjs
import plato_pylib.utils.job_running_functs as jobRun



class DosLabel():
	"""Label holding strings used to differentiate different DosRunner/Analyser objects

	Attributes:
		eleKey (str): key representing the element (or possible elements) involved.
		structKey (str): key representing the structure involved
		methodKey (str): key representing the calculation method used

	"""
	def __init__(self, eleKey=None, structKey=None, methodKey=None):

		setattr(self, "eleKey", eleKey)
		setattr(self, "structKey", structKey)
		setattr(self, "methodKey", methodKey)

		self.reqArgs = ["eleKey", "structKey", "methodKey"]
		for x in self.reqArgs:
			if getattr(self,x) is None:
				raise ValueError("{} is a required argument for DosLabel".format(x))



	@property
	def components(self):
		""" List of components contributing to identity of this label. At time of writing this is eleKey, structKey and methodKey. Consistent ordering not garanteed.
		
		"""
		return [getattr(self,x) for x in self.reqArgs]


	def __eq__(self, other):
		for attr in self.reqArgs:
			try:
				if getattr(self,attr) != getattr(other,attr):
					return False
			except AttributeError: #If other is an inappropriate type (e.g an integer) this should trigger
				return False
		return True


	def __hash__(self):
		return hash( tuple([getattr(self,x) for x in self.reqArgs]) )


class DosRunnerComposite(baseObjs.DosRunnerBase):

	def __init__(self, inpObjs):
		self.objs = inpObjs
		self._ensureNoDuplicateLabels()


	@property
	def dosComms(self):
		outList = list()
		for x in self.objs:
			outList.extend(x.dosComms)
		return outList

	@property
	def singlePointEnergyComms(self):
		outList = list()
		for x in self.objs:
			outList.extend(x.singlePointEnergyComms)
		return outList

	@property
	def label(self):
		outList = list()
		for x in self.objs:
			outList.extend(x.label)
		return outList

	def runDosGeneratingCalcs(self, nCores=1):
		jobRun.executeRunCommsParralel(self.dosComms,nCores)

	def runSinglePointEnergyCalcs(self, nCores=1):
		jobRun.executeRunCommsParralel(self.singlePointEnergyComms, nCores)

	def _ensureNoDuplicateLabels(self):
		if len(self.label) != len(set(self.label)):
			raise ValueError("Duplicate labels detected")

	def createAnalyser(self):
		inpObjs = list()
		for x in self.objs:
			inpObjs.append( x.createAnalyser() )
		return DosAnalyserComposite(inpObjs)

class DosRunnerPlato(baseObjs.DosRunnerBase):

	def __init__(self, platoCalcObj, smearWidth, stepSize, label):
		""" Initialiser for DosRunnerPlato
		
		Args:
			platoCalcObj (calc_methods.CalcObj): object which holds details for running the SPE calculation
			smearWidth (float): Smearing width for the DoS in eV
			stepSize (float): size of steps to take in eV (e.g. 0.1 means output a DoS value for every 0.1 eV)
			label (DosLabel): object which holds labels to help identify/classify/separate different jobs
				
		"""
		self._calcObj = platoCalcObj 
		self._smearWidth = smearWidth
		self._stepSize = stepSize
		self._label = label

	@property
	def dosComms(self):
		return [dosHelp.getDosPlotData(self._baseFilePath +'.occ', self._smearWidth, self._stepSize)]

	@property
	def singlePointEnergyComms(self):
		return [self._calcObj.getRunComm()]

	@property
	def _baseFilePath(self):
		filePath = self._calcObj.filePath
		basePath = os.path.splitext(filePath)[0]
		return basePath

	@property
	def label(self):
		return self.label

	def runSinglePointEnergyCalcs(self, nCores=1):
		self._calcObj.writeFile()
		jobRun.executeRunCommsParralel( self.singlePointEnergyComms, nCores )

	def runDosGeneratingCalcs(self, nCores=1):
		jobRun.executeRunCommsParralel( self.dosComms, nCores )

	def _getDosDosDataDictAfterRunningCalcs(self):
		return dosHelp.getDosPlotData(self._baseFilePath + ".occ", self._smearWidth, self._stepSize, runDos=False)


	def createAnalyser(self):
		dosDataDict = self._getDosDosDataDictAfterRunningCalcs()
		dosData = dosDataDict["dosdata"]
		eFermi = dosDataDict["efermi"]
		return DosAnalyserStandard(dosData, eFermi, self._label)



class DosAnalyserComposite(baseObjs.DosAnalyserBase):

	def __init__(self, inpObjs):
		self.objs = inpObjs



class DosAnalyserStandard(baseObjs.DosAnalyserBase):

	def __init__(self, dosData, eFermi, label):
		""" Initialiser for DosAnalyser for plato
		
		Args:
			dosData (np array, nx2): first column is energies (assumed eV), second column is density of states
			eFermi (float): Fermi energy in same units as energies (which are first col of dosData)
			label (DosLabel object):
	
		"""
		self.data = np.array(dosData)
		self.refData = None
		self.eFermi = eFermi
		self.label = label
		self.dataPlotter = None #TODO: Will make a standard dataplotter for DoS plots

	def getObjectsWithComponents(self, components, caseSensitive=True):
		if self._analyserConsistentWithInpComponentLabels(components,caseSensitive=caseSensitive):
			return [self]
		else:
			return list()

	def attachRefData(self, dosData, components, errorIfNoMatches=True, caseSensitiveComponents=True):
		if self._analyserConsistentWithInpComponentLabels(components, caseSensitiveComponents):
			self.refData = np.array(dosData)
		else:
			if errorIfNoMatches:
				raise ValueError("Reference data with components {} not matched".format(components))


	def _analyserConsistentWithInpComponentLabels(self, components, caseSensitive=True):
		inpComponents = list(components)
		allComps = list(self.label.components)
		if not caseSensitive:
			inpComponents = [x.lower() for x in inpComponents]
			allComps = [x.lower() for x in allComps]

		for x in inpComponents:
			if x not in allComps:
				return False

		return True


	def plotData(self, **kwargs):
		return None

#TODO: Probably accept label as an input argument
class DosOptions():
	def __init__(self, methodStr, structKey, eleKey, kpts, startFolder, uCellStruct, smearWidth, stepSize, modelDataFolder, integGrid = None):
	   
		self.eleKey = eleKey
		self.structKey= structKey 
		self.methodStr = methodStr
		self.kpts = kpts
		self.startFolder = startFolder
		self.smearWidth = smearWidth
		self.stepSize = stepSize
		self.integGrid = integGrid
		self.struct = uCellStruct
		
		#Create the method
		self.platoMethodObj = calcMethods.createPlatoMethodObj(self.methodStr)
		self.platoMethodObj.kpts = kpts
		self.platoMethodObj.dataSet = modelDataFolder
		if integGrid is not None:
			self.platoMethodObj.integGrid = integGrid

			
	def createDosObj(self):
		inpPath = os.path.abspath( os.path.join(self.startFolder,self.methodStr,"spe_for_dos.in") )
		strDict = self.platoMethodObj.getStrDictWithStruct(self.struct)
		runCommFunct = self.platoMethodObj.runCommFunction
		calcObj =  calcMethods.getPlatoCalcObjFromInpPathAndStrDictAndRunCommFunction(inpPath, strDict, runCommFunct)
		label = baseObjs.DosLabel(eleKey=self.eleKey, structKey=self.structKey, methodKey=self.methodStr)
		return DosRunnerPlato(calcObj, self.smearWidth, self.stepSize, label)




