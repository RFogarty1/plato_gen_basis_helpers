
import itertools as it


def getSlicesForMergingTrajectories(stepIndices, nSteps, overlapStrat="simple"):
	""" Takes a list of ORDERED [start,end] indices and returns the slice indices needed to merge trajectories
	
	Args:
		stepIndices: Ordered iter of [start,end] where start and end are indices. For example you might have a trajectory with steps [0,5,10] in which case start=0, end=10
		nSteps: (iter of int) Number of values in each trajectory. For the example above (steps=[0,5,10] this would equal 3 

	overlapStrat options:
		None: Will throw an error if theres ANY overlap between adjacent stepIndices
		"simple": Allows the start step of one trajectory to be EQUAL or greater than the end step of the previous. For example for steps=[0,5], [5,10] we use step 5 from only the second object 

 
	Returns
		 outSlices: (iter of len-2 iters) These are the slices needed for each trajectory. For example, two trajectories with steps [1,2],[3,4,5] would return [0:1], [0:2]
 
	Raises:
		 Errors
	"""

	stepSlices = list()
	if overlapStrat is None:
		return _getSlicesForNoOverlapCase(stepIndices, nSteps)
	elif overlapStrat=="simple":
		return _getSlicesForSimpleOverlapCase(stepIndices, nSteps)
	else:
		raise ValueError("{} is an invalid value for overlapStrat".format(overlapStrat))
	return stepSlices

def _getSlicesForNoOverlapCase(stepIndices, nSteps):
	stepSlices = [ [0,nStep] for nStep in nSteps]

	for idx,unused in enumerate(stepIndices[1:],start=1):
		if (stepIndices[idx][0]-stepIndices[idx-1][1]) < 1:
			stepA, stepB = stepIndices[idx-1], stepIndices[idx] 
			raise ValueError("Overlapping trajectories appear; the problem is in steps {} and {}".format(stepA, stepB))
	return stepSlices

def _getSlicesForSimpleOverlapCase(stepIndices, nSteps):
	stepSlices = list()
	
	stepSlices.append( [0,nSteps[0]] )

	for idx,unused in enumerate(stepIndices[1:],start=1):
		currStart, prevEnd = stepIndices[idx][0], stepIndices[idx-1][1] 
		if (stepIndices[idx][0]-stepIndices[idx-1][1]) < 0:
			stepA, stepB = stepIndices[idx-1], stepIndices[idx] 
			raise ValueError("Overlapping trajectories appear; the problem is in steps {} and {}".format(stepA, stepB))
		elif (stepIndices[idx][0]-stepIndices[idx-1][1])==0:
			stepSlices[-1][1] -= 1
			currSlices = [0, nSteps[idx]]
			stepSlices.append(currSlices)
		else:
			currSlices = [0, nSteps[idx]]
			stepSlices.append(currSlices)

	return stepSlices



def trimTrajectoriesIfRequired(orderedTrajs, trimStrat):
	if trimStrat is None:
		return None

	#Figure out the indices we need to keep
	trajObjs = [x.trajSteps for x in orderedTrajs]
	trajSteps = list()
	for trajObj in trajObjs:
		currSteps = [x.step for x in trajObj]
		trajSteps.append(currSteps)
	sliceIndices = getSliceIndicesForTrimmingTrajectories(trajSteps, trimStrat=trimStrat)

	#Trim the trajectories	
	for sIndices, trajObj in it.zip_longest(sliceIndices, orderedTrajs):
		trajObj.trajSteps = trajObj.trajSteps[slice(sIndices[0],sIndices[1])]


def getSliceIndicesForTrimmingTrajectories(trajSteps, trimStrat="simple"):
	""" Gets the slices required for trimming trajectories such that they dont overlap
	
	Args:
		trajSteps: (iter of len-N iters) 

	trimStrat values:
		"simple": From start->end traj removes steps from early trajectories which would overlap the next one. (e.g. steps=[0,5,10], [5,10,15] becomes [0],[5,10,15])
			 
	Returns
		idxSlices: iter of len-2 iters. Each contains the slice of step indices we need for one trajectory e.g. [[0,4],[0,3]] means take the first 5 for the first traj, and the first 4 for the second traj
 
	"""
	if trimStrat is None:
		return [[0,len(x)] for x in trajSteps]

	if trimStrat=="simple":
		return _getSliceIndicesForTrimming_simpleStrat(trajSteps)
	else:
		raise ValueError("{} is an invalid value for trimStrat".format(trimStrat))


def _getSliceIndicesForTrimming_simpleStrat(trajSteps):

	outSlices = list()
	for idx, tSteps in enumerate(trajSteps):
		if idx==len(trajSteps)-1:
			outSlices.append( [0, len(tSteps)] )
		else:
			maxIdx = trajSteps[idx+1][0]
			for endIdx,stepIdx in enumerate( reversed(tSteps) ):
				if stepIdx<maxIdx:
					sliceEnd = len(tSteps)-endIdx
					break
			outSlices.append( [0, sliceEnd] )

	return outSlices

