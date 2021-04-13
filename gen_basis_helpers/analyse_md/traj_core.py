
import copy
import itertools as it
import json
import os
import pathlib
import plato_pylib.shared.ucell_class as uCellHelp

from . import shared_misc as miscHelp

class TrajectoryBase():
	"""Object representing the total trajectory of a simulation. Main job is to yield a generator/iterable of TrajStepBase objects

	"""

	def __iter__(self):
		raise NotImplementedError("")


class TrajectoryInMemory(TrajectoryBase):

	def __init__(self, trajSteps):
		""" Initializer
		
		Args:
			trajSteps: (iter of TrajStepBase) Each object represents one step of an MD run
				 
		"""
		self.trajSteps = trajSteps

	def __iter__(self):
		return iter(self.trajSteps)

	@classmethod
	def fromDict(cls, inpDict):
		trajSteps = [TrajStepBase.fromDict(x) for x in inpDict["trajSteps"]]
		return cls(trajSteps)

	def toDict(self):
		outDict = {"trajSteps":[x.toDict() for x in self.trajSteps]}
		return outDict

	def applyFunctToEachTrajStep(self, funct):
		""" Applies a function to each trajectory step. Useful for doing things like converting length units
		
		Args:
			funct: f(TrajStepBase) Modifies the object in place
				 
		"""
		for tStep in self.trajSteps:
			funct(tStep)

	def __eq__(self, other):
		if self.trajSteps!=other.trajSteps:
			return False
		return True


class TrajStepBase():
	""" Represents a single step in a trajectory. Used to define the minimal expected attributes

	Attributes:
		unitCell: (plato_pylib UnitCell) Contains geometry at this step (central cell atoms only)
		step: (int) Contains the step number
		time: (float) Contains the time at this step number
	"""
	def __init__(self, unitCell=None, step=None, time=None):
		self._eqTol = 1e-5
		self.normalCmpAttrs = ["unitCell", "step"]
		self.numericalCmpAttrs = ["time"]
		self.numericalArrayCmpAttrs = list()
		self.unitCell = unitCell
		self.step = step
		self.time = time

	@classmethod
	def fromDict(cls, inpDict):
		useDict = {k:v for k,v in inpDict.items() if k!="unitCell"}
		if inpDict.get("unitCell",None) is not None:
			useDict["unitCell"] = uCellHelp.UnitCell.fromDict(inpDict["unitCell"])
		return cls(**useDict)

	def toDict(self):
		outDict = {"step":self.step, "time":self.time, "unitCell":self.unitCell.toDict()}
		return outDict

	def __eq__(self, other):
		eqTol = min(self._eqTol, other._eqTol)

		if sorted(self.normalCmpAttrs)!=sorted(other.normalCmpAttrs):
			return False

		if sorted(self.numericalCmpAttrs)!=sorted(other.numericalCmpAttrs):
			return False

		for attr in self.normalCmpAttrs:
			if getattr(self,attr)!=getattr(other,attr):
				return False

		for attr in self.numericalCmpAttrs:
			valA, valB = getattr(self,attr), getattr(other,attr)
			if (valA is None) and (valB is None):
				pass
			elif (valA is None) or (valB is None):
				return False
			elif abs(valA-valB)>eqTol:
				return False

		#This section is for sub-classes really
		if sorted(self.numericalArrayCmpAttrs) != sorted(other.numericalArrayCmpAttrs):
			return False

		#Numerical array comparisons	
		for attr in self.numericalArrayCmpAttrs:
			arrA, arrB = getattr(self,attr), getattr(other,attr)
			if self._checkTwoNumericalArraysEqual(arrA, arrB, eqTol) is False:
				return False

		return True

	#Needed for more general classes that inherit from this
	def _checkTwoNumericalArraysEqual(self, arrA, arrB, eqTol):

		if (arrA is None) and (arrB is None):
			pass

		elif (arrA is None) or (arrB is None):
			return False

		else:
			if len(arrA) != len(arrB):
				return False
			
			for rowA, rowB in it.zip_longest(arrA, arrB):
				if len(rowA) != len(rowB):
					return False
	
				for valA, valB in it.zip_longest(rowA, rowB):
					absDiff = abs(valA-valB)
					if absDiff > eqTol:
						return False

		return True


class TrajStepFlexible(TrajStepBase):

	def __init__(self, unitCell=None, step=None, time=None, extraAttrDict=None):
		""" Initializer. NOTE: This is a pretty general/abstract initializer which isnt really meant to be called directly at high level

			unitCell: (plato_pylib UnitCell) Contains geometry at this step (central cell atoms only)
			step: (int) Contains the step number
			time: (float) Contains the time at this step number
			extraAttrDict: (dict of dicts) See below

		extraAttrDict: This allows new attributes to be specified for each traj step. The keys are the names of the attributes you want (e.g. "velocities") The keys that should be set in the other dict are all optional but likely VERY IMPORTANT:

			value: The value to set this attribute on upon initiation. Default=None	
			cmpType: (str) Values indicate the type of equality comparison to use. Allowed values are "normal","numerical", "numericalArray". Default is "normal", meaning self.attr==other.attr gets tested directly (works for ints but not floats)

		"""
		super().__init__(unitCell=unitCell, step=step, time=time)
		self.reservedAttrs = ["unitCell", "step", "time", "_eqTol", "normalCmpAttrs", "numericalCmpAttrs", "numericalArrayCmpAttrs"]
		extraAttrDict = dict() if extraAttrDict is None else extraAttrDict
		self.extraAttrs = set()

		self.addExtraAttrDict(extraAttrDict)

	def _appendCmpTypeToRelevantList(self, attrName, cmpType):
		if cmpType=="normal":
			self.normalCmpAttrs.append(attrName)
		elif cmpType=="numerical":
			self.numericalCmpAttrs.append(attrName)
		elif cmpType=="numericalArray":
			self.numericalArrayCmpAttrs.append(attrName)
		else:
			raise ValueError("cmpType = {} for attrName={} is an invalid value".format(cmpType, attrName))

	def addExtraAttrDict(self, extraAttrDict):
		""" Adds extra attributes to the object. Used (for example) to add velocities to steps in a trajectory
		
		extraAttrDict: This allows new attributes to be specified for each traj step. The keys are the names of the attributes you want (e.g. "velocities") The keys that should be set in the other dict are all optional but likely VERY IMPORTANT:

			value: The value to set this attribute on upon initiation. Default=None	
			cmpType: (str) Values indicate the type of equality comparison to use. Allowed values are "normal","numerical", "numericalArray". Default is "normal", meaning self.attr==other.attr gets tested directly (works for ints but not floats)

		"""
		if extraAttrDict is not None:

			#Check extraAttrDict doesnt have any reservedAttr vals
			for key in extraAttrDict.keys():
				if key in self.reservedAttrs:
					raise ValueError("Cant set {} from extraAttrDict; this is a reseverd attribute name".format(key))

			#Set the values for any new attrs
			for key in extraAttrDict.keys():
				currDict = extraAttrDict[key]
				setattr(self, key, currDict.get("value",None))
				self.extraAttrs.add(key)
				if currDict.get("cmpType") is not None:
					self._appendCmpTypeToRelevantList(key, currDict.get("cmpType"))
				else:
					self._appendCmpTypeToRelevantList(key, "normal")


	def toDict(self):
		outDict = super().toDict()
		for key in self.extraAttrs:
			currVal = getattr(self, key)
			try:
				outDict[key] = currVal.toDict()
			except AttributeError:
				outDict[key] = currVal

		#Dump info on how to compare the keys
		outDict["normalCmpAttrs"] = self.normalCmpAttrs
		outDict["numerical"] = self.numericalCmpAttrs
		outDict["numericalArray"] = self.numericalArrayCmpAttrs

		return outDict

	@classmethod
	def fromDict(cls, inpDict):
		#Deal with the standard stuff
		outDict = dict()
		if inpDict.get("unitCell",None) is not None:
			outDict["unitCell"] = uCellHelp.UnitCell.fromDict(inpDict["unitCell"])

		for key in ["step","time"]:
			outDict[key] = inpDict.get(key,None)

		#Deal with extra attributes
		cmpAttrs = ["normalCmpAttrs", "numerical", "numericalArray"]

		extraAttrDict = dict()
		for key in inpDict.keys():
			if (key not in outDict) and (key not in cmpAttrs):
				extraAttrDict[key] = {"value": inpDict[key]}
				if key in inpDict["normalCmpAttrs"]:
					extraAttrDict[key]["cmpType"] = "normal"
				elif key in inpDict["numerical"]:
					extraAttrDict[key]["cmpType"] = "numerical"
				elif key in inpDict["numericalArray"]:
					extraAttrDict[key]["cmpType"] = "numericalArray"


		return cls(extraAttrDict=extraAttrDict,**outDict)


def dumpTrajObjToFile(trajObj, outFile):
	""" Dump TrajectoryBase to file. Format involves writing each step as a dict (i.e. JSON notation).
	
	Args:
		trajObj: (TrajectoryBase object) Contains all steps in the trajectory
		outFile: (str) Path to the output file
	 
	"""
	outDir = os.path.split(outFile)[0]
	pathlib.Path(outDir).mkdir(parents=True, exist_ok=True)
	with open(outFile,"wt") as f:
		for step in trajObj:
			f.write(json.dumps(step.toDict()))
			f.write("\n")

def readTrajObjFromFileToTrajectoryInMemory(inpFile):
	""" Reads a trajectory file (format defined by dumpTrajObjToFile) and returns a TrajectoryInMemory object
	
	Args:
		inpFile: (str) Path to file containing the trajectory
			 
	Returns
		 trajObj: (TrajectoryInMemory) Holds the full trajectory in memory. Simple/flexible but may not be ideal for large simulations

	"""
	outTrajObjs = list()
	with open(inpFile,"rt") as f:
		for line in f:
			currDict = json.loads(line)
#			currObj = TrajStepBase.fromDict(currDict)
			currObj = TrajStepFlexible.fromDict(currDict)
			outTrajObjs.append(currObj)

	return TrajectoryInMemory(outTrajObjs)

def readLastTrajStepFromFile(inpFile):
	""" Reads only the final trajectory step from a trajectory file (format defined by dumpTrajObjToFile). Faster than reading the whole thing and taking the last step
	
	Args:
		inpFile: (str) Path to the *.traj file
			 
	Returns
		outStep: (TrajStepBase) Contains information (including geometry) for the final step in an MD trajectory
 
	"""
	with open(inpFile,"rt") as f:
		line = getFinalNLinesFromFileObj(f)[-1]
		currDict = json.loads(line)
#		outObj = TrajStepBase.fromDict(currDict)
		outObj = TrajStepFlexible.fromDict(currDict)
	return outObj


#Taken from:
#https://stackoverflow.com/questions/136168/get-last-n-lines-of-a-file-similar-to-tail
#Original function was called tail
def getFinalNLinesFromFileObj(f, lines=1, _buffer=4098):
    """Tail a file and get X lines from the end"""
    # place holder for the lines found
    lines_found = []

    # block counter will be multiplied by buffer
    # to get the block size from the end
    block_counter = -1

    # loop until we find X lines
    while len(lines_found) < lines:
        try:
            f.seek(block_counter * _buffer, os.SEEK_END)
        except IOError:  # either file is too small, or too many lines requested
            f.seek(0)
            lines_found = f.readlines()
            break

        lines_found = f.readlines()

        # we found enough lines, get out
        # Removed this line because it was redundant the while will catch
        # it, I left it for history
        # if len(lines_found) > lines:
        #    break

        # decrement the block counter to get the
        # next X bytes
        block_counter -= 1

    return lines_found[-lines:]


#TODO: Add an option for checking the timestep spacing stays constant (finalTime/nSteps) == time(step1) - time(step0)
def getMergedTrajInMemory(trajList, overlapStrat="simple", trimStrat="simple"):
	""" Returns a TrajectoryInMemory which is the merge of trajectories in trajList
	
	Args:
		trajList: (iter of TrajectoryInMemory objects)
		mergeStrat: (str or None) How to handle the cases where some steps are present in multiple trajectories. None means just throw an error if this is the case
		trimStrat: (str or None) How to handle case where trajectory steps overlap (e.g steps=[0,5,10], [5,10,15])

	mergeStrat values:
		"simple": Allows the start step of one trajectory to be EQUAL or greater than the end step of the previous. For example for steps=[0,5], [5,10] we use step 5 from only the second object 

	trimStrat values:
		"simple": From start->end traj removes steps from early trajectories which would overlap the next one. (e.g. steps=[0,5,10], [5,10,15] becomes [0],[5,10,15])
 
	Returns
		outTraj: (TrajectoryInMemory) Contains all the ordered trajectories in trajList
 
	WARNING:
		a) This function doesnt involve COPYING anything, since that would be too inefficient for many typical cases. Thus modifying the output from this function will also modify the trajSteps in trajList.
		b) The input trajectories are "trimmed" in place. Which may be confusing if you plan on using the untrimmed trajectories
		c) "Trimming" is applied before the overlap strat is. This may affect what you get depending on the options used

	Raises:
		 ValueError: If step numbers in trajList overlap between two trajectories (e.g. if theres a "step 5" in two of the input trajectories) AND mergestrat=None
	"""
	#Step 1 = order by step number.
	startSteps = [ min(x.trajSteps,key=lambda a:a.step).step for x in trajList ]
	endSteps = [ max(x.trajSteps, key=lambda a:a.step).step for x in trajList ]

	orderedIdxVsStartSteps = sorted( [x for x in enumerate(startSteps)], key=lambda x:x[1] )


	#Step 1.5 = Trim trajectories if needed
	orderedTrajs = [trajList[idx] for idx,startStep in orderedIdxVsStartSteps]
	miscHelp.trimTrajectoriesIfRequired(orderedTrajs, trimStrat)


	#Step 1.9 - figure out new start/end steps
	startSteps = [ min(x.trajSteps,key=lambda a:a.step).step for x in trajList ]
	endSteps = [ max(x.trajSteps, key=lambda a:a.step).step for x in trajList ]

	orderedIdxVsStartSteps = sorted( [x for x in enumerate(startSteps)], key=lambda x:x[1] )

	#Step 2 - figure out start/end step indices for trimmed trajectories
	orderedTrajs = list()
	stepIndices = list()

	for idx,unused in orderedIdxVsStartSteps:
		orderedTrajs.append( trajList[idx] )
		stepIndices.append( (startSteps[idx],endSteps[idx]) )

	nSteps = [len(traj.trajSteps) for traj in orderedTrajs]

	#Step 2 = Deal with overlapping steps between trajectories
	stepSlices = miscHelp.getSlicesForMergingTrajectories(stepIndices, nSteps, overlapStrat=overlapStrat)

	#Step 3 = create a new object with merged trajectories
	outTrajSteps = list()
	for stepSlice,currTraj in zip(stepSlices,orderedTrajs):
		outTrajSteps.extend( currTraj.trajSteps[slice(*stepSlice)] )

	return TrajectoryInMemory(outTrajSteps)



def getTimeAndGeomClosestToInpTimeForInpTraj(inpTime, inpTrajInMem, equiDist=1e-2, prioritiseLate=True, convAngToBohr=False):
	""" Extracts a geometry from a trajectory closest to the input time
	
	Args:
		inpTime: (float) The time you want to get a geometry for
		inpTrajInMem: (TrajectoryInMemory) Input trajectory. MUST BE IN ORDER (small time to large time)
		equiDist: (float) Two steps are considered EQUAL if they differ from inpTime by less than this amount
		prioritiseLate: (Bool) If True then we take the latest step when two are equal.
		convAngToBohr: (Bool) If True then apply .convAngToBohr() to the output geometry
 
	Returns
		outTime: (float) The timestep 
		outGeom: (plato_pylib UnitCell obj) Geometry at outTime
 
	Notes:
		TrajectoryInMemory must be ordered by time; which SHOULD always be the case anyway but..
	"""
	#Find the index
	for idx,tStep in enumerate(inpTrajInMem):
		currDiff = abs(inpTime - tStep.time)

		if idx==0:
			minDiff, maxDiff, minIndices = currDiff, currDiff, [idx]
		else:
			
			if currDiff < minDiff+equiDist:
				#Case 1: We're definitely closer than before
				if abs(currDiff)+equiDist < abs(minDiff):
					minDiff, minIndices = currDiff, [idx]
				else:
					minIndices.append(idx)

			#If we're getting furthe away from the input time, we're never gonna get any closer
			elif (currDiff > minDiff):
				break

			#Seems very unlikely to trigger (unless two steps are equidistant)
			else:
				pass


	minIdx = minIndices[-1] if prioritiseLate else minIndices[0]

	#Use the indices to get the geometry; this only works this efficiently in trajInMemory
	# (in other cases we'd need to consume a new iterator until we reached this step)
	outGeom = copy.deepcopy( inpTrajInMem.trajSteps[minIdx].unitCell )
	outTime = inpTrajInMem.trajSteps[minIdx].time
	if convAngToBohr:
		outGeom.convAngToBohr()
	
#	print("Quenched at {} fs".format(outTime))
	return outTime, outGeom

