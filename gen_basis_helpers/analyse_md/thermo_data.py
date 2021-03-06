
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



def getMergedStandardThermoData(dataList, overlapStrat="simple"):
	""" Returns a ThermoDataStandard object made by merging those in dataList
	
	Args:
		dataList: (iter of ThermoDataStandard objs) They should all have the same properties
		overlapStrat: Keyword for how to deal with overlapping steps; details handled by miscHelp.getSlicesForMergingTrajectories	
 
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


