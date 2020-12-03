
import math

from ..shared import creator_resetable_kwargs as resetableKwargFactory

class GetStatsForThermoProps(resetableKwargFactory.CreatorWithResetableKwargsTemplate):
	registeredKwargs = set(resetableKwargFactory.CreatorWithResetableKwargsTemplate.registeredKwargs)

	registeredKwargs.add("startTime") #Inclusive
	registeredKwargs.add("endTime") #Inclusive i guess
	registeredKwargs.add("timeTol") #Absolute tolerance when checking if two times are equal 
	registeredKwargs.add("startStep") #Inclusive; overrides startTime
	registeredKwargs.add("endStep") #Inclusive; overrides endTime
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

		if self.startStep is not None:
			deltaStep = [x-self.startStep for x in steps]
			startIdx = len(steps) - len([x for x in deltaStep if x>0])
		if self.endStep is not None:
			deltaStep = [x-self.endStep for x in steps]
			endIdx = len(steps) - len([x for x in deltaStep if x>0])

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

