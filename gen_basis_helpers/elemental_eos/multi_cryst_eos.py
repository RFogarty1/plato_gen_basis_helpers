
''' Purpose of these objects are to hold data on the eos fits to multiple crystal structures '''

import numpy as np 

import gen_basis_helpers.shared.plot_functions as pltFuncts


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


	#TODO: Maybe remove the label property
	@property
	def label(self):
		""" Str: label representing the object """
		raise NotImplementedError("")


	@property
	def structLabels(self):
		""" List of strs: Labels for all the structures """
		raise NotImplementedError("")

	@property
	def elementLabel(self):
		""" Str: The element that the pure-crystals are; simply needs to be the same for all structs so it should
		be possible to use non-elemental structures too (just put a compound label here)
		"""
		raise NotImplementedError("")

	@property
	def methodLabel(self):
		""" Str: Label representing the method """
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

	def __init__(self,singleCrysts:"SingleCrystEosResult objs", label):
		self._singCrysts = list(singleCrysts)
		self._label = label

		self._ensureStructLabelsAllUnique()
		self._ensureMethodLabelsAllTheSame()
		self._ensureElementLabelsAllTheSame()

	def _ensureStructLabelsAllUnique(self):
		allStructLabels = self.structLabels
		if len(allStructLabels) != len(set(allStructLabels)):
			raise ValueError("Duplicate structLabels present in list {}".format(allStructLabels))

	def _ensureMethodLabelsAllTheSame(self):
		methodLabels = [x.methodLabel for x in self._singCrysts]
		numbUnique = len(set(methodLabels))
		if numbUnique != 1:
			raise ValueError("{} unique method labels found in input objects".format(numbUnique)) 

	def _ensureElementLabelsAllTheSame(self):
		allElementLabels = [x.elementLabel for x in self._singCrysts]
		numbUnique = len(set(allElementLabels))
		if numbUnique != 1:
			raise ValueError("{} unique element labels found in input objects".format(numbUnique))

	@property
	def label(self):
		return self._label

	@property
	def elementLabel(self):
		return self._singCrysts[0].elementLabel

	@elementLabel.setter
	def elementLabel(self,value):
		for x in self._singCrysts:
			x.elementLabel = value

	@property
	def methodLabel(self):
		return self._singCrysts[0].methodLabel

	@property
	def structLabels(self):
		outList = list()
		for x in self._singCrysts:
			outList.extend(x.structLabels)
		return outList

	def _getCombinedDictFromAllInputObjs(self, prop):
		outDict = dict()
		for x in self._singCrysts:
			currDict = getattr(x,prop)
			outDict.update(currDict)
		return outDict

	def getTableData(self,numbDp=3):
		outDict = dict()
		for x in self._singCrysts:
			currData = x.getTableData(numbDp=numbDp)
			outDict.update(currData)

		return outDict

	def getPlotData(self):
		outDict = dict()
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


class SingleCrystEosResult( EosDataHolderOneElementAndMethod ):
	def __init__(self,v0=None, e0=None, b0=None, fitData=None, actData=None,
				 structLabel=None, elementLabel=None, methodLabel=None):
		self._v0 = v0
		self._b0 = b0
		self._e0 = e0
		self.fitData = fitData
		self.actData = actData
		self._structLabel = structLabel
		self._elementLabel = elementLabel
		self._methodLabel = methodLabel

	def getPlotData(self):
		return {self._structLabel: np.array(self.actData)}

	def getTableData(self,numbDp=3):
		numbFormat = "{:." + str(numbDp) + "f}"
		outData = [self.methodLabel] + [numbFormat.format(x) for x in [self._v0, self._b0, self._e0]]
		return {self._structLabel: outData}

	@property
	def e0(self):
		return {self._structLabel:self._e0}

	@property
	def v0(self):
		return {self._structLabel:self._v0}

	@property
	def b0(self):
		return {self._structLabel:self._b0}

	@property
	def structLabels(self):
		return [self._structLabel]

	@property
	def methodLabel(self):
		return self._methodLabel

	@property
	def elementLabel(self):
		return self._elementLabel

	@elementLabel.setter
	def elementLabel(self,value):
		self._elementLabel = value

	@property
	def tableHeadings(self):
		return ["Method","v0","b0","e0"]

#	@property
#	def deltaE0(self):
#		return {self.structLabel:0.0}









class GroupedMultiCrystEosForOneElement():
    def __init__(self, multiCrystResults):
        self._multiCrystResults = multiCrystResults
        self._ensureAllElementKeysTheSame()
    
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
    
    def createPlots(self, refStr=None, structOrder=None, deltaE0=True):
        if refStr is None:
            pass
        else:
            refMethod = refStr
        
        allPlots = list()
        for currObj in self._multiCrystResults:
            methodStr = currObj.methodLabel
            if refStr is None:
                refMethod = methodStr
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
        
        xLabel = "Volume per atom / $a_{0}^{3}$"
        yLabel = "$\Delta$E per atom / eV"
        title = self.elementLabel.replace("_"," ") + " " + methodLabel.replace("_"," ")
        structLabels = [x.replace("_"," ") for x in structOrder]
        modelLabels = [x.replace("_"," ") for x in [refStr,methodLabel]]
        currFig = pltFuncts.createEnergyVsVolCurves_multiStructEachMethod( [refData, plotData], title=title,structLabels=structLabels,
                                                                           modelLabels=modelLabels,xlabel=xLabel,ylabel=yLabel )
        return currFig
    
    def _getPlotDataDictOneMethod(self,methodLabel,deltaE0=True):
        #Step 1 = get the object for this method
        relObjs = list()
        for x in self._multiCrystResults:
            if x.methodLabel==methodLabel:
                relObjs.append(x)
                
        assert(len(relObjs)==1),"Only 1 set of results can be present for a method."
    
        #Step 2 = get the plot data dict
        dataDict = relObjs[0].getPlotData()
        if deltaE0:
            e0ValDict = relObjs[0].e0
            minE0 = min(e0ValDict.values())
            print("minE0 = {}".format(minE0))
            for key in dataDict.keys():
                dataDict[key][:,1] -= minE0
        return dataDict
    
    
    @property
    def elementLabel(self):
        return self._multiCrystResults[0].elementLabel
    
    
    @property
    def _uniqueStructLabels(self):
        allStructLabels = list()
        for x in self._multiCrystResults:
            allStructLabels.extend(x.structLabels)
        return list(set(allStructLabels))
 
