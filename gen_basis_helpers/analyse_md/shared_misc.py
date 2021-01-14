

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

