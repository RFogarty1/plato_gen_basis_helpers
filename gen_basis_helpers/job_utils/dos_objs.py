
import os

import numpy as np

from ..shared import calc_methods as calcMethods
from ..shared import label_objs as labelsBase
from ..shared import misc_utils as misc
from . import dos_helpers as dosHelp
from . import dos_base_objs as baseObjs
from . import dos_data_plotter as dosPlotter
import plato_pylib.utils.job_running_functs as jobRun




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
			label (labelsBase.StandardLabel): object which holds labels to help identify/classify/separate different jobs
				
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

	@misc.getObjectsWithComponentsInstanceWrapper(isComposite=True)
	def __init__(self, inpObjs):
		self.objs = inpObjs

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

	@misc.getObjectsWithComponentsInstanceWrapper(isComposite=False)
	def __init__(self, dosData, eFermi, label):
		""" Initialiser for DosAnalyser for plato
		
		Args:
			dosData (np array, nx2): first column is energies (assumed eV), second column is density of states
			eFermi (float): Fermi energy in same units as energies (which are first col of dosData)
			label (labelsBase.StandardLabel object):
	
		"""
		self.data = np.array(dosData)
		self.refData = None
		self.eFermi = eFermi
		self.label = [label]
		self.dataPlotter = dosPlotter.DataPlotterDos.fromDefaultPlusKwargs()


	def attachRefData(self, dosData, components, errorIfNoMatches=True, caseSensitiveComponents=True):
		if len(self.getObjectsWithComponents(components, caseSensitiveComponents)) == 1:
			self.refData = np.array(dosData)
		else:
			if errorIfNoMatches:
				raise ValueError("Reference data with components {} not matched".format(components))



	def plotData(self, **kwargs):
		thisFunctKwargs = ["inclRefData", "extraData", "shiftDataEFermiToZero", "refDataFirst"]
		ourData = np.array(self.data)
		inclRefData = kwargs.get("inclRefData",True)
		refDataFirst = kwargs.get("refDataFirst",False)


		#Apply shift to normal data if requested
		if kwargs.get("shiftDataEFermiToZero", False):
			shiftedData = np.array(ourData)
			shiftedData[:,0] = shiftedData[:,0] - self.eFermi
			ourData = shiftedData

		plotData = [ourData]

		#self.data
		if (self.refData is not None) and (inclRefData):
			if refDataFirst:
				plotData = [self.refData, ourData]
			else:
				plotData = [ourData, self.refData]

		#Add ref/extra data to plotData
		extraData = kwargs.get("extraData",list())
		plotData.extend(extraData)



		#Set the title str (if usr hasnt)
		titleStr = "{}-{}-{}".format(self.label[0].eleKey, self.label[0].structKey, self.label[0].methodKey)
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
		label = labelsBase.StandardLabel(eleKey=self.eleKey, structKey=self.structKey, methodKey=self.methodStr)
		return DosRunnerPlato(calcObj, self.smearWidth, self.stepSize, label)




