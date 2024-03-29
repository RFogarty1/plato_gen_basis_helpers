
import copy
import itertools as it
import math
import numpy as np



def getPESFromOneDimPosAndForceSimple(posVsForce):
	""" Gets an estimated potential energy surface from a list of positions and forces
	
	Args:
		posVsForce: (iter of len-2 iters) For each element x[0] is the position while x[1] is the force
			 
	Returns
		posVsEnergy: (iter of len-2 iters) For each element x[0] is the position while x[1] is the energy. Units match the energy units used in the forces

	Notes:
		Simply does linear interpolation for everything.
 
	"""
	outVals = [ [posVsForce[0][0], 0] ]#Set the zero to the first data point passed

	for idx,(pos,force) in enumerate(posVsForce[1:],start=1):
		dist = pos - posVsForce[idx-1][0]
		forceSum = force + posVsForce[idx-1][1]
		averageForce = forceSum / 2
		extraEnergy = averageForce*dist
		currEnergy = outVals[idx-1][1] + extraEnergy
		outVals.append( [pos,currEnergy] )

	return outVals

def getEstimateOfUnderlingForcesFromOneDimCVPositionData(posVals, timeVals, hillsInfo, massCV, frictionCoeff=0):
	""" Gets an estimate of the underlying forces at each position of the collective variable. Originally made for a collective variable involving an atom moving from a fixed plane
	
	Uses the formula [dE/dr = frictionCoeff*v + massCV*accel - dE_{p}/dr] where dE/dr is the force from the unmodified potential energy surface, v is the velocity of the collective variable, and E_{p} is the energy of the added potential (contained in hillsInfo) 

	Args:
		posVals: (iter of floats) Positions of the collective variable
		timeVals: (iter of floats) Same length as posVals; the times relating to each position
		hillsInfo: (MetadynHillsInfo object) Contains information on the potential added to the simulation at any one point
		massCV: (float) The mass associated with the collective variable
		frictionCoeff: (float) Friction coefficient used in Langevin dynamics; this will cause a drag force that we factor out. Not sure how to set for more complicated collective variables than movement of a single atom but....
			 
	Returns
		outData: (iter of len-2 iters) x-values are positions, y-values are estiamtes of the force at that position. Order determined from input otder
 
	NOTES:
		UNITS: The unit of mass will depend on units for position, hill energies and time used I guess. For hills in eV, positions in bohr and time in femtoseconds you should use atomicMass*protonMassInKg*joulesToEv* (metreToBohr**-2)*(secondToFs**2) )
		LENGTH(Output): We get velocities and accelerations from posVals/times
		VELOCITIES/ACCELS: We use the average values of these over a timestep
		TimeVals: I assume these are all +ve (or zero) and this is a sort of key assumption. So dont have negative times please

	"""
	#Get velocities/accelerations/times for as much data as possible
	vels, accels = _getVelocitiesAndAccelFromPositionsAndTimes(posVals, timeVals)
	trimmedPositions = posVals[2:]
	trimmedTimes = timeVals[2:]

	#Figure out positions vs potential for these values
	outForces = list()
	for idx in range(len(accels)):
		currPos, currTime = trimmedPositions[idx], trimmedTimes[idx]
		currPotObj = hillsInfo.createGroupedHills(timeRange=[-1,currTime])
		currPenaltyGradient = currPotObj.evalGradientAtVals( [[currPos]] )[0][0] 
		currForce = massCV*accels[idx] - frictionCoeff*vels[idx] - currPenaltyGradient
		outForces.append(currForce)

	return [ [pos,force] for pos,force in it.zip_longest(trimmedPositions,outForces) ]


def _getVelocitiesAndAccelFromPositionsAndTimes(inpPositions, inpTimes):
	allVel, accel = list(), list()
	for idx in range(1,len(inpPositions)):
		deltaT = inpTimes[idx] - inpTimes[idx-1]
		velocity = (inpPositions[idx] - inpPositions[idx-1]) / deltaT
		allVel.append(velocity)
		if idx > 1:
			acceleration = (allVel[-1] - allVel[-2]) / deltaT
			accel.append(acceleration)

	return allVel[1:], accel


def getTimeVsPotAtValsForEachHillAdded(inpObj, inpVals, timeRange=None, timeTol=1e-4, minTimeDiff=1e-3, oneDimOutput=True):
	""" Gets [ [x,....z], [val] ] for the potential at each time a hill was added.
	
	Args:
		inpObj: (MetadynHillsInfo obj)
		inpVals: ( iter of len-n iters) n is the number of dimensions where we want to evaluate the Hills at
		timeRange: (len-2 iter) Range of times to include potentials for; default is to just include all potentials spawned between -1 and np.inf
		timeTol: (float) Tolerance parameter for timeRange
		minTimeDiff: (float) If multiple hills were added between t and t+minTimeDiff, then they are considered as all being added at time t
		oneDimOutput: (Bool) If True AND this is a 1-dimensional case then we will return [ [xA,valA], [xB,valB] ] in data fields (outVals) instead of [ [[xA],valA], [[xB],valB] ]; which is less awkward to plot 

	Returns
		outVals: e.g. [ [timeA,dataA], [timeB,dataB] ]. Length of iterable is number of hills added. Each element represents the potential at one step (i.e. time when 1 hill is added). Each data is a len-n iter. Within each, the first elements are positions in the n-dimensional collective variable space; i.e. the length of this iter will be the number of collective variable dimensions. An example would be dataA= [ [[xA,yA], valA], [[xB,yB], valB] ]
 
	"""
	#Get time + potential objs
	timeVsPotObjs = _getTimeVsPotObjsForEachTimeHillAdded(inpObj, timeRange=timeRange, timeTol=timeTol, minTimeDiff=1e-3)
	outTimes = [ x[0] for x in timeVsPotObjs ]
	outPotObjs = [ x[1] for x in timeVsPotObjs ]

	#Calculate the actual potentials at each time
	outData = list()
	for pot in outPotObjs:
		currData = pot.evalFunctAtVals(inpVals)
		if len(inpVals[0])==1 and oneDimOutput:
			currData = [ [x[0],y] for x,y in it.zip_longest(inpVals,currData) ]
		else:
			currData = [ [x,y] for x,y in it.zip_longest(inpVals,currData) ]

		outData.append(currData)

	return [ [t,data] for t,data in it.zip_longest(outTimes,outData)]

#This is going to be much slower than it COULD be. I could alternatively get the pot obj for the full time and sum its contribs
#in certain ways to get total potential at earlier times; that would be faster but probably trickier to test
def _getTimeVsPotObjsForEachTimeHillAdded(inpObj, timeRange=None, timeTol=1e-4, minTimeDiff=1e-3):
	""" Returns iter of [time,potObj] for each time a hill was added in a metadynamics simulation
	
	Args:
		inpObj: (MetadynHillsInfo obj)
		timeRange: (len-2 iter) Range of times to include potentials for; default is to just include all potentials spawned between -1 and np.inf
		timeTol: (float) Tolerance parameter for timeRange
		minTimeDiff: (float) If multiple hills were added between t and t+minTimeDiff, then they are considered as all being added at time t
		 
	Returns
		timeVsPot: (iter of len-2 iters) Each element has [time, potObj] where potObj is a GroupedMultiDimGaussHills representing all potentials added up to that time

	NOTE:
		I've not really tested how timeTol and minTimeDiff interact. I wouldnt try to do anything too clever with them tbh and minTimeDiff>timeTol is probably advisable
		Also; setting timeTol=0 will probably totally break things; so keep it sensible
 
	"""
	useObj = copy.deepcopy(inpObj)
	useObj.sortTimes()

	#Step 0: Figure out the time range
	if timeRange is None:
		timeRange = [-1, np.inf]

	#Step 1: We need to figure out the time ranges for each
	minTime, maxTime = timeRange
	spawnTimes = [t for t in useObj.times if t>minTime and t<maxTime] #this will be in order



	#Find the INDICES corresponding to times where hills were essentially spawned
	spawnIndices, outTimeVals = [0], [spawnTimes[0]]
	for idx,t in enumerate(spawnTimes[1:],1):
		if t>outTimeVals[-1]+minTimeDiff:
			spawnIndices.append(idx)
			outTimeVals.append(spawnTimes[idx])

	#Now use minTime to get timeRanges (we need to do a load of -1 things tho)
	outTimeRanges = list()
	for idx in spawnIndices[1:]:
		outTimeRanges.append( [minTime,spawnTimes[idx-1]] )

	#We need to figure out the last one now; which just runs from minTime to maxTime
	outTimeRanges.append([minTime, maxTime])

	#Step 2: get potential objs from these times
	outObjs = list()
	for tRange in outTimeRanges:
		currObj = useObj.createGroupedHills(timeRange=tRange, timeTol=timeTol)
		outObjs.append(currObj)

	return [ [t,obj] for t,obj in it.zip_longest(outTimeVals,outObjs)]


def evalPotAddedOverTimeRangeForHillsInfoObj(inpObj, evalAtVals, timeRange=None, timeTol=1e-3):
	""" Calculates the potential added at a set of colvar coords
	
	Args:
		inpObj: (MetadynHillsInfo) Contains all info on metadynamics hills and the times they were added
		evalAtVals: (iter of len-n iters) The values at which to calculate the potential. For 1 dimension this may be [ [0], [1], [2] ], for 2-dimension it may be [ [1,2], [3,4], [5,6] ]
		timeRange: (len-2 iter) Range of times to include potentials for; default is to just include all potentials
		timeTol: (float) Tolerance parameter for timeRange

	Returns
		outVals: (iter of floats) Length is the same as len(evalAtVals); This gives the potential at each input point
 
	"""
	gauFunctGrouped = inpObj.createGroupedHills(timeRange=timeRange, timeTol=timeTol)
	outVals = gauFunctGrouped.evalFunctAtVals(evalAtVals)
	return outVals



def addBlankDimensionToMetadynHillsInfoInstance(inpInstance, scaleVals=0, posVals=0):
	""" Adds a dimension of blank valued peaks (heights=same as first dim) to a MetadynHillsInfo instance. Original purpose was to use hills from a 1-dimension-colvar simulation in a 2-dimension-colvar simulation
	
	Args:
		inpInstance: (MetadynHillsInfo object)
		scaleVals: (float) The scale of the blank hills to add. CP2K interprets scale=0 as an infinitely flat Gaussian, which presumably stops the blank dimension interfering with the other one

	Returns
		Nothing; works in place
 
	"""
	attrsToAppendTo = ["heights", "positions", "scales"]

	for idx in range( len(inpInstance.times) ):
		inpInstance.heights[idx].append( inpInstance.heights[idx][0] )
		inpInstance.positions[idx].append( posVals )
		inpInstance.scales[idx].append( scaleVals )


def getMergedMetadynHillsInfoInstance(inpInstances, copyData=True):
	""" Takes iter of MetadynHillsInfo instances and returns a merged version
	
	Args:
		inpInstances: (iter of MetadynHillsInfo objects)
		copyData: (Bool) If True then the new merged object will be totally independent of those it comes from. Safer option usually.

	Returns
		outObj: (MetadynHillsInfo instance) Merged version of inputs.
 
	"""
	#Combine lists
	attrsToCombo = ["times", "positions", "scales", "heights"]
	outIters = list()
	for attr in attrsToCombo:
		if copyData:
			currList = copy.deepcopy( getattr(inpInstances[0], attr) )
		else:
			currList = getattr(inpInstances[0], attr)		


		for instance in inpInstances[1:]:
			if copyData:
				currExt = copy.deepcopy( getattr(instance, attr) )
			else:
				currExt = getattr(instance, attr)
			currList.extend( currExt )
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


	def addHill(self, time=None, positions=None, scales=None, heights=None):
		""" Add a new hill to the object
		
		Args: (ALL REQUIRED; despite being keyword args)
			time: (float) The time this hill was added to the simulation
			positions: (iter of floats) Number of dimensions should match self.positions
			scales: (iter of floats) Number of dimensions should match self.scales
			heights: (iter of floats) Number of dimensions should match self.heights

		Returns
			What Function Returns
	 
		Raises:
			AttributeError: If any of the input args are None
			ValueError: If input hill is the wrong dimension in any of positions/scales/heights
			ValueError: If input positions/scales/heights dont all have the same dimensionality
		"""
		#Check all arguments have been set to something
		inpArgs = [time,positions,scales,heights]
		inpLabels = ["time", "positions", "scales", "heights"]
		if any([x is None for x in [time,positions,scales,heights]]):
			idx = inpArgs.index(None)
			raise AttributeError("None is an invalid value for {}".format(inpLabels[idx]))

		#Check the number of dimensions is correct
		currObjNumbDims = len(self.positions[0])
		newDims = [len(positions), len(scales), len(heights)]
		if not all([x==newDims[0] for x in newDims]):
			raise ValueError("Dimensions for positions/scales/heights should all be the same but are {},{},{}".format(*newDims))

		if newDims[0] != currObjNumbDims:
			raise ValueError("Cant add an {}-dimension hill to a {}-dimension object".format(newDims[0],currObjNumbDims))

		self.times.append(time)
		self.positions.append(positions)
		self.scales.append(scales)
		self.heights.append(heights)

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

	def createNewObjFromLimitedTimeRange(self, timeRange=None, timeTol=1e-4):
		""" Creates a new instance with peaks outside timeRange filtered out
		
		Args:
			timeRange: (len-2 iter) [minTime, maxTime]
			timeTol: (float) Two times are considered the same if their values are within "timeTol"; purpose is to deal with float errors

		Returns
			outObj: (MetadynHillsInfo) Object with filtered hills
	 
		"""
		outIndices = self._getIndicesWithinTimeRange(timeRange, timeTol)
		outKeys = ["times", "positions", "scales", "heights"]
		outDict = dict()
		for key in outKeys:
			currVals = getattr(self, key)
			outDict[key] = [ currVals[idx] for idx in outIndices ]

		return MetadynHillsInfo(**outDict)

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

	@classmethod
	def fromDict(cls, inpDict):
		return cls(**inpDict)


	def toDict(self):
		outAttrs = ["times", "positions", "scales", "heights"]
		outDict = dict()
		for attr in outAttrs:
			outDict[attr] = getattr(self, attr)
		return outDict


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

	def evalGradientAtVals(self, posVals):
		""" Gets the gradient vector at an iter of positions
		
		Args:
			posVals: (Iter of len-n floats; where n is number of dimensions) For example, for two dimensions this may be [ [1,2], [3,4], [5,6] ] which will evaluate for 3 positions
				 
		Returns
			grads: (iter of len-n floats; where n is the number of dimension) Each element is one gradient vector. Each element of the gradient vector is the derivative in one of the dimensions
	
		WARNING:
			I havent tested this for multiple dimensions
 
		"""
		#Get the individual gradient vectors
		individualGradientVectors = list()
		for hill in self.multiDimGaus:
			currGradVectors = hill.evalGradientAtVals(posVals)
			individualGradientVectors.append( currGradVectors )

		#Combine vectors
		nDims = len(currGradVectors[0])
		outVectors = list()
		for pIdx, pos in enumerate(posVals):
			currGradVector = list()
			for dimIdx in range(nDims):
				currVal = sum([x[pIdx][dimIdx] for x in individualGradientVectors])
				currGradVector.append(currVal)
			outVectors.append(currGradVector)

		return outVectors

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

	def evalGradientAtVals(self, posVals):
		""" Gets the gradient vector at an iter of positions
		
		Args:
			posVals: (Iter of len-n floats; where n is number of dimensions) For example, for two dimensions this may be [ [1,2], [3,4], [5,6] ] which will evaluate for 3 positions
				 
		Returns
			grads: (iter of len-n floats; where n is the number of dimension) Each element is one gradient vector. Each element of the gradient vector is the derivative in one of the dimensions
	 
		"""
		derivVals = list()
		for idx,oneDimHill in enumerate(self.oneDimHills):
			currPosVals = [x[idx] for x in posVals]
			derivVals.append( oneDimHill.evalNthDerivativeAtVals(currPosVals,n=1) )

		outGrads = [ [*vals] for vals in it.zip_longest(*derivVals) ]
		return outGrads

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

	#TODO: Use hermite polynomials to easily extend this to higher order derivatives
	def evalNthDerivativeAtVals(self, posVals, n=1):
		""" Get values for the n-th derivative at a set of input parameters
		
		Args:
			posVals: (iter of floats) Each represents an input position at which to evaluate the derivative
			n: (int,Optional) The n-th derivative to take. 0 means return the function, 1 is the first derivative etc.
				 
		Returns
			outVals: (iter of float)  Values of the n-th derivative at the input posVals
	 
		"""
		zerothDerivVals = self.evalFunctAtVals(posVals)
		if n==0:
			outVals = zerothDerivVals
		elif n==1:
			prefactor = -1/(self.scale**2)
			outVals = [ prefactor*(x-self.pos)*gau for x,gau in it.zip_longest(posVals,zerothDerivVals) ]

		else:
			raise NotImplementedError("")

		return outVals

	def __eq__(self, other):
		eqTol = min(self._eqTol, other._eqTol)
		floatAttrs = ["height", "scale", "pos"]

		for attr in floatAttrs:
			valA, valB = getattr(self,attr), getattr(other,attr)
			if abs(valA-valB)>eqTol:
				return False

		return True



