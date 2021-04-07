
import copy

from . import traj_core as trajCoreHelp

def getTrajSampledEveryNSteps(inpTraj, sampleEveryN, inPlace=False, createView=True):
	""" Samples input trajectory every n-steps; useful for making trajectories smaller before analysis
	
	Args:
		inpTraj: (TrajectoryInMemory) MD trajectory
		sampleEveryN: (int) How frequently to sample. e.g. if set to 10 we take every 10 steps 
		inPlace: (Bool) Whether to do this in place or not; if False a new trajectory is returned, else the trajectory is updated in place
		createView: (Bool) Whether the output trajectory references the same steps as inpTraj. This is more memory efficient and cheaper to do (no copying involved). If True then altering a step from outTraj will also alter inpTraj, which may lead to bugs if the output trajectory needs modifying later
			 
	Returns
		outTraj: (TrajectoryInMemory) The trajectory sampled every N steps. Only returned if inPlace is False
 
	Raise:
		ValueError: If both createView and inPlace are set to True

	"""
	if inPlace and createView:
		raise ValueError("inPlace and createView are both set to True; at least one of these needs to be set to False for now")


	outSteps = list()
	for idx,step in enumerate(inpTraj):
		if idx%sampleEveryN == 0:
			if createView:
				outSteps.append(step)
			else:
				outSteps.append(copy.deepcopy(step))


	if inPlace:
		inpTraj.trajSteps = outSteps
	else:
		return trajCoreHelp.TrajectoryInMemory(outSteps)



def getTrajSplitIntoEqualSections(inpTraj, nStepsEach, createView=True):
	""" Get input trajectory split into equal sections
	
	Args:
		inpTraj: (TrajectoryInMemory) MD trajectory
		nStepsEach: Number of steps in each output trajectory object
		createView: (Bool) Whether the output trajectories reference the same steps as inpTraj. This is more memory efficient and cheaper to do (no copying involved).
			 
	Returns
		outTrajs: (iter of TrajectoryInMemory objs) Each has nStepsEach individual steps; if the total number of steps doesnt divide cleanly then the final steps will just be ignored (i.e. the combination of outTrajs usually wont be as long as inpTraj)
 
	"""
	outSteps = list()

	allSteps = inpTraj.trajSteps

	for idx, step in enumerate(allSteps):
		if (idx%nStepsEach == 0) and (idx!=0):
			if createView:
				outSteps.append( allSteps[idx-nStepsEach:idx]  )
			else:
				outSteps.append( copy.deepcopy(allSteps[idx-nStepsEach:idx])  )

	return [trajCoreHelp.TrajectoryInMemory(steps) for steps in outSteps]


