
import copy
import itertools as it
import math

from ..shared import creator_resetable_kwargs as resetableKwargFactory
from ..shared import data_plot_base as dPlotBase
from ..shared import table_maker_base as tableMakerHelp

class ThermoDashboardSimple():
	""" Class for creating summary data from a ThermoDataStandard instance

	"""
	def __init__(self, props, startIdx=0, movingAvg=True, templatePlotter=None, centralWindowMovingAvg=False,
	             centralWindowWidthEachSide=1):
		""" Initializer
		
		Args:
			props: (iter of str) The properties to include in the dashboard outputs (e.g. ["temp"])
			startIdx: (int) The step index (not value) at which to start plotting/calculating averages
			movingAvg: (Bool) If True calculate and plot a moving average (no effect on the tables)
			templatePlotter: (dPlotBase.DataPlotterStandard object) Used as a template to plot
			centralWindowMovingAvg: (Bool) Whether
			centralWindowWidthEachSide: (int) The number of steps each side to use for the central window moving average. Only relevant if centralWindowMovingAvg=True
				 
		"""
		self.props = props
		self.startIdx = startIdx #Apply even to non-avg vals
		self.movingAvg = movingAvg
		self.templatePlotter = self._getDefaultPlotter() if templatePlotter is None else templatePlotter
		self.centralWindowMovingAvg = centralWindowMovingAvg
		self.centralWindowWidthEachSide = centralWindowWidthEachSide

	def _getDefaultPlotter(self):
		return dPlotBase.DataPlotterStandard(lineStyles=['none','-','-'], lineMarkers=['x','',''], xlabel="Time", lineColors=["r","k","b"])

	def getPlotFactories(self, thermoObj):
		outFactories = list()
		for prop in self.props:
			currFactory = self._getPlotFactoryOneProp(thermoObj, prop)
			outFactories.append(currFactory)
		return outFactories

	def _getPlotFactoryOneProp(self, thermoObj, prop):
		#1) Get the data needed
		outData = list()
		time = thermoObj.dataDict["time"][self.startIdx:]
		yData = thermoObj.dataDict[prop][self.startIdx:]
		rawData = [ [x,y] for x,y in it.zip_longest(time,yData)	]
		outData.append(rawData)
		if self.movingAvg:
			movingAvgYData = getMovingAverageFromThermoData(thermoObj, prop, startIdx=self.startIdx)
			avgData = [ [x,y] for x,y in it.zip_longest(time,movingAvgYData) ]
			outData.append(avgData)

		if self.centralWindowMovingAvg:
			currArgs = [thermoObj, prop, self.centralWindowWidthEachSide]
			centMovingAvg = getSimpleCentralWindowAverageFromThermoData(*currArgs, startIdx=self.startIdx)
			centMovingData = [ [x,y] for x,y in it.zip_longest(time, centMovingAvg) if x is not None]
			outData.append(centMovingData)

		#2) Create the output factory
		outFactory = copy.deepcopy(self.templatePlotter)
		outFactory.data = outData
		outFactory.ylabel = prop
		return outFactory

	def getThermoSummaryTable(self, thermoObj, fmt="html"):
		#1) Get headers
		headers = self._getThermoSummaryTableHeaders()
	
		#2) Get data
		getStatsObj = GetStatsForThermoProps(props=self.props, startIdx=self.startIdx)
		statsDict = getStatsObj.create(data=thermoObj)
		outData = list()
		for prop in self.props:
			currRow = [prop] + self._getTableStrListFromStatsObjOneProp(statsDict[prop])
			outData.append(currRow)
	
		#3) Create table
		tableMaker = tableMakerHelp.TableMakerStandard(fmt=fmt, headers=headers, mapInputToRowsFunct=lambda x:x)
		return tableMaker.createTable(outData)

	def _getTableStrListFromStatsObjOneProp(self, statsObj):
		attrs = ["mean", "standardDev", "minVal", "maxVal", "range", "nVals"]
		outList = list()
		for attr in attrs:
			currVal = getattr(statsObj,attr)
			currStr = "{:.0f}".format(currVal)
			outList.append(currStr)
		return outList

	def _getThermoSummaryTableHeaders(self):
		return ["Property", "Mean", "std-dev", "min", "max", "range", "nVals"]


def getMovingAverageFromThermoData(thermoData, prop, startIdx=0):
	""" 
	
	Args:
		thermoData: ThermoDataStandard object
		prop: (str) the property the moving average is needed for
		startIdx: (int) The idx in thermoData to start calculating the average value at
			 
	Returns
		movingAvg: (float) Values of the moving average for each idx after startIdx (inclusive)
 
	"""
	currSum = 0
	outVals = list()
	for idx, val in enumerate(thermoData.dataDict[prop][startIdx:], start=1):
		currSum += val
		currVal = currSum/idx
		outVals.append(currVal)
	return outVals

def getSimpleCentralWindowAverageFromThermoData(thermoData, prop, widthEachSide, fillVal=None, startIdx=0):
	""" 
	
	Args:
		thermoData: ThermoDataStandard object
		prop: (str) the property the moving average is needed for
		widthEachSide: (int) Number of data points to take each side when calculating the mean
		fillVal: The value we output when we cant calculate a central moving average (e.g. we cant get a moving average for the first or last data points)
		startIdx: (int) The idx in thermoData to start calculating the average value at
			 
	Returns
		movingAvg:  (float) Values of the moving average for each idx after startIdx (inclusive). Note: Will contain "fillVal" values
 
	"""
	inpIter = thermoData.dataDict[prop][startIdx:]
	args = [inpIter, widthEachSide]
	kwargs = {"fillValWhenNotCalc":fillVal}
	return _getSimpleCentralWindowAverageFromIter(*args, **kwargs)

def _getSimpleCentralWindowAverageFromIter(inpIter, widthEachSide, fillValWhenNotCalc=None):
	""" Gets a central-window moving average for the input data (each mean value is calculated from n values either side) 
	
	Args:
		inpIter: (iter of Numbers) We take the moving averages of these numbers
		widthEachSide: (int) Number of data points to take each side when calculating the mean
		fillValWhenNotCalc: The value we output when we cant calculate a central moving average (e.g. we cant get a moving average for the first or last data points)

	Returns
		outVals: (iter of float) Central moving average. Each data point is an average of 2*widthEachSide + 1 data points
 
	"""
	stack = list()
	outVals = list()
	lenIter = len(inpIter)
	reqStackSize = 2*widthEachSide + 1

	#Initialise our stack
	stack = [x for x in inpIter[:widthEachSide]]

	#Calculate moving averages
	for idx,val in enumerate(inpIter):

		#Deal with the stack
		if len(stack) == reqStackSize:
			stack.pop(0)
		
		if idx < lenIter-widthEachSide:
			stack.append(inpIter[idx+widthEachSide])

		#Calculate moving average
		if len(stack) < reqStackSize:
			outVals.append(fillValWhenNotCalc)
		else:
			outVals.append( sum(stack)/len(stack) )
		
	return outVals



class GetStatsForThermoProps(resetableKwargFactory.CreatorWithResetableKwargsTemplate):
	registeredKwargs = set(resetableKwargFactory.CreatorWithResetableKwargsTemplate.registeredKwargs)

	registeredKwargs.add("startTime") #Inclusive
	registeredKwargs.add("endTime") #Inclusive i guess
	registeredKwargs.add("timeTol") #Absolute tolerance when checking if two times are equal 
	registeredKwargs.add("startStep") #Inclusive; overrides startTime
	registeredKwargs.add("endStep") #Inclusive; overrides endTime
	registeredKwargs.add("startIdx") #Inclusive; overrides startStep
	registeredKwargs.add("endIdx") #Inclusive; overrides endStep
	registeredKwargs.add("props")
	registeredKwargs.add("data") #ThermoObj obviously

	def _createFromSelf(self):
		eqTol = 1e-5
		props = self.props if self.props is not None else ["temp","pressure"]
		startIdx, endIdx = self._getIndices()

		outDict = dict()
		for prop in props:
			currRawData = self.data.dataDict[prop][startIdx:endIdx]
			outDict[prop] = self._getStatsObjForOneProp( currRawData )
		return outDict	

	def _getIndices(self):
		steps = self.data.dataDict["step"]
		time = self.data.dataDict["time"]
		assert len(steps)==len(time)
		timeTolerance = 1e-6 if self.timeTol is None else self.timeTol
		startIdx,endIdx = 0, len(steps)

		#Set based on time THEN on steps; means steps will overwrite time if both set
		#(which is the behaviour i want)
		if self.startTime is not None:
			deltaT = [x-self.startTime for x in time]
			startIdx =  len(time) - len([x for x in deltaT if x>-1*timeTolerance])
		if self.endTime is not None:
			deltaT = [x-self.endTime for x in time]
			endIdx = len(time) - len([x for x in deltaT if x>timeTolerance]) 

		#steps override time
		if self.startStep is not None:
			deltaStep = [x-self.startStep for x in steps]
			startIdx = len(steps) - len([x for x in deltaStep if x>0])
		if self.endStep is not None:
			deltaStep = [x-self.endStep for x in steps]
			endIdx = len(steps) - len([x for x in deltaStep if x>0])

		#Directly passed startIdx,endIdx overrides all others
		if self.startIdx is not None:
			startIdx = self.startIdx
		if self.endIdx is not None:
			endIdx = self.endIdx

		return startIdx,endIdx

	def _getStatsObjForOneProp(self, inpData):
		nVals = len(inpData)
		mean = sum(inpData)/nVals
		minVal, maxVal = min(inpData), max(inpData)
		stdDev = math.sqrt( sum([(x-mean)**2 for x in inpData])/nVals )
		kwargDict = {"mean":mean, "standardDev":stdDev, "minVal":minVal, "maxVal":maxVal, "nVals":nVals}
		return StandardStatsObject(**kwargDict)


class StandardStatsObject():

	def __init__(self, **kwargs):
		""" Initializer
		
		kwargs (All default to None):
			mean: (float)
			standardDev: (float) 
			minVal: (float) Lowest value obtained
			maxVal: (float) Highest value obtained
			nVals: (int) Number of values used to calculate properties
				 
		"""
		self._eqTol = 1e-5
		self.posVals = ["mean", "standardDev", "minVal", "maxVal", "nVals"]
		for attr in self.posVals:
			setattr(self, attr, kwargs.get(attr,None))


	@property
	def range(self):
		if (self.minVal is None) or (self.maxVal is None):
			outVal = None
		else:
			outVal = abs(self.maxVal-self.minVal)
		return outVal

	def __eq__(self, other):
		eqTol = min(self._eqTol, other._eqTol)

		for attr in self.posVals:
			valA, valB = getattr(self,attr), getattr(other,attr)
			if (valA is None) and (valB is None):
				pass
			elif (valA is None) or (valB is None):
				return False
			elif abs(valA-valB) > eqTol:
				return False

		return True

