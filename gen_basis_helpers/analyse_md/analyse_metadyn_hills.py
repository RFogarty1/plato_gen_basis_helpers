
import itertools as it
import math
import numpy as np


def getMergedMetadynHillsInfoInstance(inpInstances):
	""" Takes iter of MetadynHillsInfo instances and returns a merged version
	
	Args:
		inpInstances: (iter of MetadynHillsInfo objects)
			 
	Returns
		outObj: (MetadynHillsInfo instance) Merged version of inputs.

	NOTE:
		Default implementation doesnt involve copying anything; this could cause issues if you keep using inpInstances after doing this merge (since multiple objects will be effectively storing the same lists of data)
 
	"""
	#Combine lists
	attrsToCombo = ["times", "positions", "scales", "heights"]
	outIters = list()
	for attr in attrsToCombo:
		currList = getattr(inpInstances[0], attr)
		for instance in inpInstances[1:]:
			currList.extend( getattr(instance, attr) )
		outIters.append(currList)

	#Create output object
	kwargDict = {k:v for k,v in it.zip_longest(attrsToCombo, outIters)}
	return MetadynHillsInfo(**kwargDict)


class MetadynHillsInfo():
	""" Class meant to store info on metadynamics hills and allow access for hills over various times; i.e. this is the class meant for manipulating the hills """

	def __init__(self, times=None, positions=None, scales=None, heights=None, sortTimes=False):
		""" Initializer
		
		Args:
			times: (iter of floats) Each is the time at which the hill was added
			positions: (iter of len-n iter of floats) Where n is number of collective variables. e.g [ [1,2], [3,4], [5,6] ] for 2 collective vars
			scales: (iter of len-n iter of floats) Where n is number of collective variables
			heights: (iter of len-n iter of floats) Where n is number of collective variables
			sortTimes: (Bool) If true then re-order inputs internally such that they are sorted in order of the times

		"""
		self.times = times
		self.positions = positions
		self.scales = scales
		self.heights = heights
		if sortTimes:
			self.sortTimes()

	def sortTimes(self):
		""" Sorts the values of self.times in ascending order, and updates over attrs accordingly
		"""
		indicesAndTimesSorted = sorted( [[idx,t] for idx,t in enumerate(self.times)], key=lambda x:x[1] )
		sortedIndices = [x[0] for x in indicesAndTimesSorted]
		attrsToRearrange = ["times", "positions", "scales", "heights"]
		for attr in attrsToRearrange:
			startVals = getattr(self,attr)
			newVals = [ startVals[idx] for idx in sortedIndices ]
			setattr(self, attr, newVals)

	def multiplyHeightsByFactor(self, factor):
		""" Multiply the stored heights by a factor; example use is unit conversion
		
		Args:
			factor: (float) The factor to multiply heights for
				 
		"""
		outHeights = list()
		for heightsOneHill in self.heights:
			currHeights = [x*factor for x in heightsOneHill]
			outHeights.append(currHeights)
		self.heights = outHeights

	def createGroupedHills(self, timeRange=None, timeTol=1e-4):
		""" Function to create "GroupedMultiDimGaussHills" object, which is useful for plotting the potential from a set of hills
		
		Args:
			timeRange: (len-2 iter of floats) [minTime, maxTime] Default is to use each hill regardless of time spawned
			timeTol: (float) Two times are considered the same if their values are within "timeTol"; purpose is to deal with float errors

		Returns
			GroupedMultiDimGaussHills: 
	 
		Raises:
			Errors
		"""
		multiDimHills = self.createMultiDimHills(timeRange=timeRange, timeTol=timeTol)
		groupedHills = GroupedMultiDimGaussHills(multiDimHills)
		return groupedHills

	def createMultiDimHills(self, timeRange=None, timeTol=1e-4):
		""" 
		
		Args:
			timeRange: (len-2 iter of floats) [minTime, maxTime] Default is to use each hill regardless of time spawned
				 
		Returns
			outHills: (iter of MultiDimGaussHill)
	 
		Raises:
			Errors
		"""
		timeIndices = self._getIndicesWithinTimeRange(timeRange, timeTol)
		positions = [self.positions[idx] for idx in timeIndices]
		scales = [self.scales[idx] for idx in timeIndices]
		heights = [self.heights[idx] for idx in timeIndices]

		args = [heights, scales, positions]
		outVals = _getIterOfMultiDimHills(*args)
		return outVals

	def getTimesWithinRange(self, timeRange=None, timeTol=1e-4):
		""" Gets times stored on the object within timeRange, in the same order as stored on the object.  Useful for doing things like relating functions from self.create* methods to a specific time a function for each time
		
		Args:
			timeRange: (len-2 iter of floats) [minTime, maxTime] Default is to use each hill regardless of time spawned
				 
		Returns
			outTimes: (iter of floats) Times within timeRange; in same order as self.times.e
	 
		Raises:
			Errors
		"""
		timeIndices = self._getIndicesWithinTimeRange(timeRange, timeTol)
		return [self.times[idx] for idx in timeIndices]

	def _getIndicesWithinTimeRange(self, minAndMaxTime, timeTol):
		outIndices = list()

		if minAndMaxTime is None:
			return [x for x in range(len(self.times))]

		minTime, maxTime = minAndMaxTime
		if maxTime < minTime:
			raise ValueError("minTime={}, maxTime={}; maxTime should be greater than minTime".format(minTime,maxTime))


		for idx,t in enumerate(self.times):
			if ((t-minTime) > -1*abs(timeTol)) and ( (maxTime-t) > -1*abs(timeTol) ) :
				outIndices.append(idx)

		return outIndices


	def __eq__(self, other):
		_eqTol = 1e-5
		if len(self.times) != len(other.times):
			return False

		for timeA, timeB in it.zip_longest(self.times,other.times):
			if abs(timeA-timeB) > _eqTol:
				return False


		hillsA, hillsB = self.createMultiDimHills(), other.createMultiDimHills()

		if hillsA != hillsB:
			return False

		return True



def _getIterOfMultiDimHills(heights, scales, positions):
	""" Gets an iterable of MultiDimGaussHill objects, with one per hill spawned
	
	Args:
		heights: (iter of len-n iter of floats) Where n is number of collective variables
		scales: (iter of len-n iter of floats) Where n is number of collective variables
		positions: (iter of len-n iter of floats) Where n is number of collective variables
			 
	Returns
		outHills: (iter of MultiDimGaussHill; same length as heights/scales/positions) 
 
	Raises:
		AssertionError: If different length heights/scales/positions are passed

	"""
	assert all([len(x)==len(heights) for x in [scales,positions]])
	outObjs = list()
	for h,s,p in it.zip_longest(heights, scales, positions):
		currObj = MultiDimGaussHill.fromIters(heights=h, scales=s, positions=p)
		outObjs.append(currObj)

	return outObjs

#TODO: Implement times as an optional argument to init?
class GroupedMultiDimGaussHills():
	""" Represents the sum of individual multi-dimensional Gaussian hills """

	def __init__(self, multiDimGaus):
		""" Initializer
		
		Args:
			multiDimGaus: (iter of MultiDimGaussHill objects)
				 
		"""
		self.multiDimGaus = multiDimGaus

	def getContribsAtPositions(self, posVals):
		""" Gets an array of contributions (1 per hill) at each position in posVal 
	
		Args:
			posVals: (Iter of len-n floats; where n is number of dimensions) For example, for two dimensions this may be [ [1,2], [3,4], [5,6] ] which will evaluate for 3 positions
 
		Returns
			outContribs: (nxm array) 1 row per position, one column per hill
	 
		"""
		#Get the values
		calcdVals = list()
		for hill in self.multiDimGaus:
			currVals = hill.evalFunctAtVals(posVals)
			calcdVals.append(currVals)

		#Organise in a way i like
		outVals = list()
		for idx in range(len(posVals)):
			currVals = [calcdVals[x][idx] for x in range(len(calcdVals))]
			outVals.append(currVals)

		return outVals
		
	def evalFunctAtVals(self, posVals):
		""" Get the function evaluated at a set of position values
		
		Args:
			posVals: (Iter of len-n floats; where n is number of dimensions) For example, for two dimensions this may be [ [1,2], [3,4], [5,6] ] which will evaluate for 3 positions
				 
		"""
		individualContribs = self.getContribsAtPositions(posVals)
		outVals = list()
		for row in individualContribs:
			outVals.append( sum(row) )
		return outVals

	def __eq__(self, other):
		if len(self.multiDimGaus) != len(other.multiDimGaus):
			return False

		for mDimGauA, mDimGauB in zip(self.multiDimGaus, other.multiDimGaus):
			if mDimGauA != mDimGauB:
				return False

		return True


class MultiDimGaussHill():
	""" Combines a series of 1-dimensional GaussianHill functions and uses their product over multiple dims. e.g. if evaluating at xy we use f(x,y) = G(x)G(y) where G(x) and G(y) are the 1-dimensional Gaussian functions over x and y """

	def __init__(self, oneDimHills):
		""" Initializer
		
		Args:
			oneDimHills: (iter of OneDimGaussianHill; one per dimension)
				 
		"""
		self.oneDimHills = list(oneDimHills)

	@classmethod
	def fromIters(cls, heights=None, scales=None, positions=None):
		""" Alternative initializer
		
		Args [NOTE: All are required, despite using keywords, length of each equals n-dimensions]:
			heights: (iter of float) Paramters related to height of Gaussians
			scales: (iter of float) Parameters related to width of Gaussians
			positions: (iter of float) Where each Gaussian is centred
 
		"""
		assert all([len(x)==len(heights) for x in [heights,scales, positions]])
		outOneDimFuncts = list()
		for height, scale, position in zip(heights, scales, positions):
			currGau = OneDimGaussianHill(height=height, scale=scale, pos=position)
			outOneDimFuncts.append(currGau)
		return cls(outOneDimFuncts)

	def evalFunctAtVals(self, posVals):
		""" Get the function evaluated at a set of position values
		
		Args:
			posVals: (Iter of len-n floats; where n is number of dimensions) For example, for two dimensions this may be [ [1,2], [3,4], [5,6] ] which will evaluate for 3 positions
				 
		"""
		#Get the values for each
		valLists = list()
		for idx,gauFunct in enumerate(self.oneDimHills):
			valLists.append( np.array(gauFunct.evalFunctAtVals([x[idx] for x in posVals])) )

		#Multiple the individuals at each point
		outArray = valLists[0]
		for idx in range(1,len(valLists)):
			outArray = outArray*valLists[idx]

		outVals = outArray.tolist()

		return outVals

	def __eq__(self, other):
		if len(self.oneDimHills) != len(other.oneDimHills):
			return False

		for hillA, hillB in zip(self.oneDimHills, other.oneDimHills):
			if hillA != hillB:
				return False

		return True


class OneDimGaussianHill():
	""" Represents a one-dimensional hill function used in metadynamics. While this is a simple Gaussian function the interface reflects the terminology (heights/scales) used in metadynamics """

	def __init__(self, height=None, scale=None, pos=None):
		""" Initializer. Gaussian function is height*exp( -0.5* (\Delta pos^2)/(scale^2) )
		
		Args: [ALL ARE REQUIRED; despite the use of keyword args]
			height: (float) See formula above; larger means greater height and area
			scale: (float) See formula above; larger value means a wider hill
			pos: (float) See formula above; where the Gaussian is centred
				 
		"""
		self._eqTol = 1e-5
		self.height = height
		self.scale = scale
		self.pos = pos
		assert all([x is not None for x in [height, scale, pos]])

	def evalFunctAtVals(self, posVals):
		outList = [None for x in posVals]
		scaleFactor = self.scale**2
		for idx,x in enumerate(posVals):
			distSqr = (posVals[idx]-self.pos)**2
			outList[idx] = self.height*math.exp(-0.5*distSqr*(1/scaleFactor))
		return outList

	def __eq__(self, other):
		eqTol = min(self._eqTol, other._eqTol)
		floatAttrs = ["height", "scale", "pos"]

		for attr in floatAttrs:
			valA, valB = getattr(self,attr), getattr(other,attr)
			if abs(valA-valB)>eqTol:
				return False

		return True



