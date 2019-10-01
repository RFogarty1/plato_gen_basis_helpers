
import os

import numpy as np

from ..shared import calc_methods as calcMethods
from . import dos_helpers as dosHelp
from . import dos_base_objs as baseObjs
from . import dos_data_plotter as dosPlotter
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

	def __init__(self, inpObjs, runDosGenerating=True, runEnergy=True):
		self.objs = inpObjs
		self._ensureNoDuplicateLabels()
		self._runDos = runDosGenerating
		self._runEnergy = runEnergy

	@property
	def runDosGenerating(self):
		return self._runDos

	@runDosGenerating.setter
	def runDosGenerating(self,value):
		assert( isinstance(value,bool) )
		self._runDos = value

	@property
	def runEnergy(self):
		return self._runEnergy

	@runEnergy.setter
	def runEnergy(self, value):
		assert( isinstance(value,bool) )
		self._runEnergy = value

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


	#NOTE: Currently writeFiles() isnt garanteed in the interface definition
	def writeFiles(self):
		for x in self.objs:
			x.writeFiles()

	def runDosGeneratingCalcs(self, nCores=1):
		if self.runDosGenerating:
			jobRun.executeRunCommsParralel(self.dosComms,nCores)

	def runSinglePointEnergyCalcs(self, nCores=1):
		if self.runEnergy:
			self.writeFiles()
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

	def __init__(self, platoCalcObj, smearWidth, stepSize, label, runDosGenerating=True, runEnergy=True):
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
		self._runDos = runDosGenerating
		self._runEnergy = runEnergy


	@property
	def runDosGenerating(self):
		return self._runDos

	@runDosGenerating.setter
	def runDosGenerating(self,value):
		assert( isinstance(value,bool) )
		self._runDos = value

	@property
	def runEnergy(self):
		return self._runEnergy

	@runEnergy.setter
	def runEnergy(self, value):
		assert( isinstance(value,bool) )
		self._runEnergy = value

	@property
	def dosComms(self):
		return [dosHelp.getDosRunComm_plato(self._baseFilePath +'.occ', self._smearWidth, self._stepSize)]

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
		return [self._label]

	def writeFiles(self):
		self._calcObj.writeFile()

	def runSinglePointEnergyCalcs(self, nCores=1):
		if self.runEnergy:
			self.writeFiles()
			jobRun.executeRunCommsParralel( self.singlePointEnergyComms, nCores )

	def runDosGeneratingCalcs(self, nCores=1):
		if self.runDosGenerating:
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

	def getObjectsWithComponents(self, components, caseSensitive=True):
		outObjs = list()
		for x in self.objs:
			outObjs.extend( x.getObjectsWithComponents(components,caseSensitive=caseSensitive) )
		return outObjs


	def attachRefData(self, dosData, components, errorIfNoMatches=True, caseSensitiveComponents=True):

		if errorIfNoMatches:
			allMatched = self.getObjectsWithComponents(components,caseSensitive=caseSensitiveComponents)
			if (len(allMatched) == 0):
				raise ValueError("Couldnt find analyser object with components {}".format(components))

		for x in self.objs:
			x.attachRefData(dosData,components,errorIfNoMatches=False,caseSensitiveComponents=caseSensitiveComponents)


	def plotData(self, **kwargs):
		outHandles = list()
		for x in self.objs:
			currHandles = x.plotData(**kwargs)
			outHandles.extend(currHandles)
		return outHandles


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
		self.dataPlotter = dosPlotter.DataPlotterDos.fromDefaultPlusKwargs()

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
		thisFunctKwargs = ["inclRefData", "extraData", "shiftDataEFermiToZero", "extraDataBeforeRef"]
		plotData = [self.data]
		inclRefData = kwargs.get("inclRefData",True)

		#Add ref/extra data to plotData
		extraData = kwargs.get("extraData",list())
		if (self.refData is not None) and (inclRefData):
			refData = [self.refData]
		else:
			refData = list()		

		if kwargs.get("extraDataBeforeRef",True):
			plotData.extend(extraData)
			plotData.extend(refData)
		else:
			plotData.extend(refData)
			plotData.extend(extraData)

		#Apply shift to normal data if requested
		if kwargs.get("shiftDataEFermiToZero", False):
			shiftedData = np.array(self.data)
			shiftedData[:,0] = shiftedData[:,0] - self.eFermi
			plotData[0] = shiftedData


		#Set the title str (if usr hasnt)
		titleStr = "{}-{}-{}".format(self.label.eleKey, self.label.structKey, self.label.methodKey)
		if "titleStr" in kwargs:
			pass
		else:
			kwargs["titleStr"] = titleStr

		#Want to try to catch any wrong keywords and throw an error here, instead of in call to createPlot
		allKwargs = list(thisFunctKwargs)
		allKwargs.extend( self.dataPlotter.registeredKwargs)
		for key in kwargs:
			if key not in allKwargs:
				raise KeyError("{} is an invalid keyword arg. valid args are {}".format(key, allKwargs))

		#Need to only pass kwargs that are relevant to dataPlotter
		for key in thisFunctKwargs:
			if key in kwargs:
				kwargs.pop(key)

		return [self.dataPlotter.createPlot(plotData, **kwargs)]

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
		label = DosLabel(eleKey=self.eleKey, structKey=self.structKey, methodKey=self.methodStr)
		return DosRunnerPlato(calcObj, self.smearWidth, self.stepSize, label)




