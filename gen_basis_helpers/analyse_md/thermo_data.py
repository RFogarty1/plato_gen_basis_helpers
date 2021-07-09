
import itertools as it
import json
import numpy as np


from . import shared_misc as miscHelp


class ThermoDataInterface():
	""" Goal  of this object is to store thermodynamic data (e.g. pressure/temperature) for varying steps in the simulation

	"""

	@property
	def props(self):
		""" The properties available. Order not gauranteed to be consistent between calls """
		raise NotImplementedError("")

	def getPropsArray(self, props):
		""" Get an array with desired properties
		
		Args:
			props: (iter of str) Each str should correspond to an entry in props
				 
		Returns
			outArray: nxlen(props) np array. Each column corresponds to data of one attribute in props
	 
		"""
		raise NotImplementedError("")


#This implementation assumes we can just read all the data in one go (i.e. that it fits into memory). However, the interface should be
#general enough for when this isnt the case
class ThermoDataStandard(ThermoDataInterface):

	def __init__(self, propDataDict):
		""" Initializer
		
		Args:
			propDataDict: (dict) Keys are strs like "temp", "pressure" etc. Values are iters containing NUMERIC values. All should be the same length; though this wont be checked
	 
		"""
		self._eqTol = 1e-5
		self.dataDict = dict(propDataDict)

	@classmethod
	def fromStdKwargs(cls, step=None, time=None, temp=None, eTotal=None, eKinetic=None,
	                  ePot=None, pressure=None):
		outDict = {"step":step, "time":time, "temp":temp, "eTotal":eTotal, "eKinetic":eKinetic,
		           "ePot":ePot, "pressure":pressure}

		dictKeys = list(outDict.keys())
		for key in dictKeys:
			if outDict[key] is None:
				outDict.pop(key)

		return cls(outDict)

	#Just to satisfy an interface
	@classmethod
	def fromDict(cls, inpDict):
		return cls(inpDict)

	def toDict(self):
		return dict(self.dataDict)

	@property
	def props(self):
		return self.dataDict.keys()

	def getPropsArray(self, props):
		inpList = [self.dataDict[prop] for prop in props]
		outVals = np.array( (inpList) ).transpose()
		return outVals

	@property
	def dataListLengthsAllEqual(self):
		listLengths = [len( self.dataDict[k] ) for k in self.props]
		if all( [x==listLengths[0] for x in listLengths] ):
			return True
		else:
			return False

	def __eq__(self, other):
		eqTol = min(self._eqTol, other._eqTol)
		#compare properties, make sure order doesnt matter
		ourProps = self.props
		otherProps = other.props
		if len(ourProps)!=len(otherProps):
			return False

		for prop in ourProps:
			if prop not in otherProps:
				return False

		#Get the numerical values are equivalent
		arrayA = self.getPropsArray(ourProps)
		arrayB = other.getPropsArray(ourProps)
		if np.allclose(arrayA,arrayB, atol=eqTol) is False:
			return False

		return True

def dumpStandardThermoDataToFile(thermoDataObj, outFile):
	""" Dump MD thermodynamic info (e.g. temperatures/pressures) to a file. File is simply json
	
	Args:
		thermoDataObj: (ThermoDataStandard)
		outFile: (str) Path to the output file
			 
	"""
	outDict = thermoDataObj.toDict()
	with open(outFile,"w") as f:
		json.dump(outDict, f)

def readThermoDataFromFile(inpFile):
	""" Reads a thermoData file (format defined by dumpStandardThermoDataToFile) into a ThermoDataStandard object
	
	Args:
		inpFile: (str) Path to the file containing thermo data
			 
	Returns
		thermoDataObj: (ThermoDataStandard)
	
	"""
	with open(inpFile,"r") as f:
		inpDict = json.load(f)
	return ThermoDataStandard.fromDict(inpDict)


def getMergedStandardThermoData(dataList, overlapStrat="simple", trimStrat="simple"):
	""" Returns a ThermoDataStandard object made by merging those in dataList
	
	Args:
		dataList: (iter of ThermoDataStandard objs) They should all have the same properties
		overlapStrat: Keyword for how to deal with overlapping steps; details handled by miscHelp.getSlicesForMergingTrajectories
		trimStrat: (str or None) How to handle case where trajectory steps overlap (e.g steps=[0,5,10], [5,10,15])

	WARNING:
		a) This function doesnt involve COPYING anything, since that would be too inefficient for many typical cases. Thus modifying the output from this function will also modify the trajSteps in trajList.
		b) The input trajectories are "trimmed" in place. Which may be confusing if you plan on using the untrimmed trajectories
		c) "Trimming" is applied before the overlap strat is. This may affect what you get depending on the options used
 
	Returns
		 outData: Single ThermoDataStandard object. Note we DONT COPY ANY DATA. Thus, care should be taken if objects in dataList are kept around
 
	"""

	#Step 1 = order by step number.
	startSteps = [ min(x.dataDict["step"]) for x in dataList ]
	endSteps = [ max(x.dataDict["step"]) for x in dataList ]

	orderedIdxVsStartSteps = sorted( [x for x in enumerate(startSteps)], key=lambda x:x[1] )

	orderedTrajs = list()
	stepIndices = list()

	for idx,unused in orderedIdxVsStartSteps:
		orderedTrajs.append( dataList[idx] )
		stepIndices.append( (startSteps[idx],endSteps[idx]) )

	nSteps = [len(x.dataDict["step"]) for x in orderedTrajs]

	#Step 1.5 - check all ThermoData is consistent
	if not all(  [set(x.props)==set(dataList[0].props) for x in orderedTrajs] ):
		raise ValueError("Not all objects in datalist have the same .prop values")
	for obj in dataList:
		assert obj.dataListLengthsAllEqual

	#Now trim if required
	trajSteps = [ dataList[idx].dataDict["step"] for idx,unused in orderedIdxVsStartSteps ]
	sliceIndices = miscHelp.getSliceIndicesForTrimmingTrajectories(trajSteps, trimStrat=trimStrat)
	for thermoObj, sliceIdx in it.zip_longest(dataList, sliceIndices):
		_trimThermoDataBasedOnSlice(thermoObj, slice(*sliceIdx) )

	#Re figure out start/end steps. TODO: Remove duplication
	startSteps = [ min(x.dataDict["step"]) for x in dataList ]
	endSteps = [ max(x.dataDict["step"]) for x in dataList ]

	orderedIdxVsStartSteps = sorted( [x for x in enumerate(startSteps)], key=lambda x:x[1] )

	orderedTrajs = list()
	stepIndices = list()

	for idx,unused in orderedIdxVsStartSteps:
		orderedTrajs.append( dataList[idx] )
		stepIndices.append( (startSteps[idx],endSteps[idx]) )

	nSteps = [len(x.dataDict["step"]) for x in orderedTrajs]

	#Step 2 = figure out based on step number
	stepSlices = miscHelp.getSlicesForMergingTrajectories(stepIndices, nSteps, overlapStrat=overlapStrat)

	#Step 3 = create a new object with merged trajectories
	outKwargDict = dict()
	for prop in orderedTrajs[0].props:
		currArray = list()
		for thermoObj,stepSlice in zip(orderedTrajs,stepSlices):
			currArray.extend( thermoObj.dataDict[prop][slice(*stepSlice)] )
		outKwargDict[prop] = currArray

	return ThermoDataStandard(outKwargDict)


def _trimThermoDataBasedOnSlice(thermoData, sliceObj):
	for prop in thermoData.dataDict.keys():
		thermoData.dataDict[prop] = thermoData.dataDict[prop][sliceObj]


def getThermoDataObjSampledEveryN(inpObj, sampleEveryN):
	""" Samples ThermoDataStandard every N steps (step is defined as a data point in the ThermoDataObj).
	
	Args:
		inpObj: (ThermoDataStandard)
		sampleEveryN: (int) How frequently to sample. e.g. if set to 10 we take every 10 steps 
			 
	Returns
		outThermo: (ThermoDataStandard)

	NOTES:
		a) Nothing gets copied by default. This isnt a problem if using the data as read-only; but be careful if planning to modify data at any point
 
	"""
	#1) figure out the indices
	assert inpObj.dataListLengthsAllEqual
	outIndices = list()
	allProps = list(inpObj.props)
	currPropVals = inpObj.dataDict[allProps[0]]  #Doesnt matter which one; since lengths are all equal

	for idx,val in enumerate(currPropVals):
		if idx%sampleEveryN == 0:
			outIndices.append(idx)

	#2) Use these indices to sample each list in dataDict
	outDataDict = dict()
	for prop in allProps:
		allVals = inpObj.dataDict[prop]
		filteredVals = [ allVals[idx] for idx in outIndices ]
		outDataDict[prop] = filteredVals

	return ThermoDataStandard(outDataDict)

