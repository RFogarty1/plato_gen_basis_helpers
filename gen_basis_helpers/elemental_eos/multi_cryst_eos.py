
''' Purpose of these objects are to hold data on the eos fits to multiple crystal structures '''

from collections import OrderedDict
import numpy as np 

from ..shared.label_objs import StandardLabel
import gen_basis_helpers.shared.plot_functions as pltFuncts
from . import data_plotter_energy_breakdowns as dPlotter

class EosDataHolderOneElementAndMethod():
	""" Holds data for equation of state fits. Can include fits to multiple structures together(but single element/method)

	"""

	@property
	def e0(self):
		""" Dict, keys=struct labels, values = energy of equilibrium structure """
		raise NotImplementedError

	@property
	def deltaE0(self):
		""" Dict; keys are structure labels (for one type of crystal); values are relative energy of the equilibrium structure
		"""
		raise NotImplementedError("")

	def getPlotData(self):
		""" Dict; keys are struct labels, values are nx2 arrays (or similar) containing volume vs energy data"""
		raise NotImplementedError("")


	@property
	def label(self):
		""" iter of Label object (see BaseLabel class)"""
		raise NotImplementedError("")


	@property
	def numbCrystStructs(self):
		""" Int: Number of crystal structures """
		raise NotImplementedError("")

	@property
	def getTableData(self):
		""" Dict: keys=struct-keys, values = list, each entry containing the value in tableHeadings for that key """
		raise NotImplementedError("")

	@property
	def tableHeadings(self):
		""" Str List: Headings for the table """
		raise NotImplementedError("")



class MultiCrystEosResult( EosDataHolderOneElementAndMethod ):

	def __init__(self,singleCrysts:"SingleCrystEosResult objs"):
		self._singCrysts = list(singleCrysts)

		self._ensureStructLabelsAllUnique()
		self._ensureMethodLabelsAllTheSame()
		self._ensureElementLabelsAllTheSame()

	def _ensureStructLabelsAllUnique(self):
		allStructLabels = [x.structKey for x in self.label]
		if len(allStructLabels) != len(set(allStructLabels)):
			raise ValueError("Duplicate structLabels present in list {}".format(allStructLabels))

	def _ensureMethodLabelsAllTheSame(self):
		methodLabels = [x.methodKey for x in self.label]
		numbUnique = len(set(methodLabels))
		if numbUnique != 1:
			raise ValueError("{} unique method labels found in input objects".format(numbUnique)) 

	def _ensureElementLabelsAllTheSame(self):
		allElementLabels = [x.eleKey for x in self.label]
		numbUnique = len(set(allElementLabels))
		if numbUnique != 1:
			raise ValueError("{} unique element labels found in input objects".format(numbUnique))

	@property
	def label(self):
		outLabels = list()
		for x in self._singCrysts:
			outLabels.extend(x.label)
		return outLabels

	def _getCombinedDictFromAllInputObjs(self, prop):
		outDict = OrderedDict()
		for x in self._singCrysts:
			currDict = getattr(x,prop)
			outDict.update(currDict)
		return outDict

	def getTableData(self,numbDp=3):
		outDict = OrderedDict()
		for x in self._singCrysts:
			currData = x.getTableData(numbDp=numbDp)
			outDict.update(currData)

		return outDict

	def getPlotData(self):
		outDict = OrderedDict()
		for x in self._singCrysts:
			outDict.update(x.getPlotData())
		return outDict

	@property
	def e0(self):
		return self._getCombinedDictFromAllInputObjs("e0")

	@property
	def deltaE0(self):
		outDict = self.e0
		minEnergy = min( outDict.values() )
		for key in outDict:
			outDict[key] -= minEnergy
		return outDict

	@property
	def tableHeadings(self):
		return ["label","v0","b0","e0"]


#Note SingleCrystEosData is a subset of EosDataHolderOneElementAndMethod
class SingleCrystEosResult( EosDataHolderOneElementAndMethod ):
	def __init__(self,v0=None, e0=None, b0=None, fitData=None, actData=None,
				 structLabel=None, elementLabel=None, methodLabel=None):
		self._v0 = v0
		self._b0 = b0
		self._e0 = e0
		self.fitData = fitData
		self.actData = actData
		self._label = StandardLabel(eleKey=elementLabel, structKey=structLabel, methodKey=methodLabel)

	def getPlotData(self):
		return {self._label.structKey: np.array(self.actData)}

	def getTableData(self,numbDp=3):
		numbFormat = "{:." + str(numbDp) + "f}"
		outData = [self._label.methodKey] + [numbFormat.format(x) for x in [self._v0, self._b0, self._e0]]
		return {self._label.structKey: outData}

	@property
	def elementLabel(self):
		return self._label.eleKey

	@elementLabel.setter
	def elementLabel(self,val):
		self._label.eleKey = val

	@property
	def e0(self):
		return {self._label.structKey:self._e0}

	@property
	def v0(self):
		return {self._label.structKey:self._v0}

	@property
	def b0(self):
		return {self._label.structKey:self._b0}

	@property
	def label(self):
		return [self._label]

	@property
	def tableHeadings(self):
		return ["Method","v0","b0","e0"]





class GroupedMultiCrystEosForOneElement():
	def __init__(self, multiCrystResults, dataPlotter=None):
		self._multiCrystResults = multiCrystResults
		self._ensureAllElementKeysTheSame()
		if dataPlotter is None:
			self._createDefaultDataPlotter()
		else:
			self.dataPlotter = dataPlotter


	def _createDefaultDataPlotter(self):
		kwargsDict = {"ylabel": "$\Delta$E per atom (eV)",
		              "lineStyles": ["-","--"],
		              "showTitle": True,
		              "legend": True,
		              "lineColors": ['b','g','orange'],
		              "lineMarkers": ['x'] }

		self.dataPlotter = dPlotter.EosEnergyDataPlotter.fromDefaultPlusKwargs(**kwargsDict)
		self.dataPlotter.methodProps.remove("dataLabels")
		self.dataPlotter.dataSeriesProps.append("dataLabels")

	def _ensureAllElementKeysTheSame(self):
		pass
	
	def createTables(self):
		allTables = list()
		for currStruct in self._uniqueStructLabels:
			currTable = [[currStruct]]
			headings = self._multiCrystResults[0].tableHeadings
			currTable.append(headings)
			for currObj in self._multiCrystResults:
				currTabContrib = self._getTableDataSingleMultiCrystUsingDeltaE(currObj)[currStruct]
				currTable.append(currTabContrib)
			allTables.append(currTable)
		return allTables
	
	
	def _getTableDataSingleMultiCrystUsingDeltaE(self, multiCryst):
		outTableDict = multiCryst.getTableData(numbDp=3)
		minE = min([float(x[-1]) for x in outTableDict.values()])
		for key in outTableDict.keys():
			outTableDict[key][-1] = "{:.3f}".format( float(outTableDict[key][-1]) - minE )
		return outTableDict
	
	def createPlots(self, refStr=None, structOrder=None, deltaE0=True, methodStr=None):
		if refStr is None:
			pass
		else:
			refMethod = refStr

		if methodStr is None:
			methods = set([x.label[0].methodKey for x in self._multiCrystResults])
		else:
			methods = [methodStr]
	
		allPlots = list()
		for currObj in self._multiCrystResults:
			methodStr = currObj.label[0].methodKey
			if refStr is None:
				refMethod = methodStr
			if methodStr in methods:
				currFig = self._createPlotOneMethod(methodStr,refStr=refMethod, structOrder=structOrder,deltaE0=deltaE0)
				allPlots.append(currFig)
		return allPlots
	
	def _createPlotOneMethod(self, methodLabel,refStr=None, structOrder=None, deltaE0=True):
		plotDataDict = self._getPlotDataDictOneMethod(methodLabel,deltaE0=deltaE0)
		refDataDict = self._getPlotDataDictOneMethod(refStr,deltaE0=deltaE0) if refStr is not None else copy.deepcopy(plotDataDict)
		plotData, refData = list(), list()
		if structOrder is None:
			structOrder = plotDataDict.keys()
			
		for struct in structOrder:
			plotData.append(plotDataDict[struct])
			refData.append(refDataDict[struct])
		
		title = self.elementLabel.replace("_"," ") + " " + methodLabel.replace("_"," ")
		structLabels = [x.replace("_"," ") for x in structOrder]


		currFig = self.dataPlotter.createPlot([refData,plotData], titleStr=title, dataLabels=structLabels)
		return currFig
	
	def _getPlotDataDictOneMethod(self,methodLabel,deltaE0=True):
		#Step 1 = get the object for this method
		relObjs = list()
		for x in self._multiCrystResults:
			if x.label[0].methodKey==methodLabel:
				relObjs.append(x)
				
		assert(len(relObjs)==1),"Only 1 set of results can be present for a method."
	
		#Step 2 = get the plot data dict
		dataDict = relObjs[0].getPlotData()
		if deltaE0:
			e0ValDict = relObjs[0].e0
			minE0 = min(e0ValDict.values())
			for key in dataDict.keys():
				currDataMin = min( dataDict[key][:,1] )
				minE0 = min([minE0,currDataMin])

			for key in dataDict.keys():
				dataDict[key][:,1] -= minE0
		return dataDict
	
	
	@property
	def elementLabel(self):
		return self._multiCrystResults[0].label[0].eleKey
	
	
	@property
	def _uniqueStructLabels(self):
		allStructLabels = list()
		for x in self._multiCrystResults:
			currLabels = [y.structKey for y in x.label]
			allStructLabels.extend(currLabels)
		return list(set(allStructLabels))
 
