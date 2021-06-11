
import copy
import itertools as it
import numpy as np

import plato_pylib.shared.unit_convs as uConvHelp
import plato_pylib.shared.ucell_class as uCellHelp

from . import traj_core as trajCoreHelp
from ..shared import simple_vector_maths as vectHelp


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


def addVelocitiesToTrajInMemNVT(inpTraj, posConvFactor=1, timeConvFactor=1, velKey="velocities_from_pos"):
	""" Calculates velocities based on differences in positions between steps and adds to the trajectory object. Note that this means that the last traj step wont have a velocity associated with it
	
	Args:
		inpTraj: (TrajectoryInMemory with TrajStepFlexible)
		posConvFactor: (float, Optional) Multiply delta(positions) by this factor when calculating velocities
		timeConvFactor: (float, Optional) Multiply delta(time) by this factor when calculating velocities
		velKey: (str, Optional) The key under which these velocities are stored

	Returns
		Nothing; works in place. inpTraj is modified so that all steps (except the last one) contain a velKey attribute. For each step this is a list of velocities matching the order of atoms in unitCell.cartCoords
 
	"""
	trajSteps = inpTraj.trajSteps
	convFactor = posConvFactor / timeConvFactor

	for idx in range(1, len(trajSteps)):
		#Figure out the velocities
		currStartCoords = trajSteps[idx-1].unitCell.cartCoords
		currEndCoords = trajSteps[idx].unitCell.cartCoords
		currDeltaT = trajSteps[idx].time - trajSteps[idx-1].time
		currVels = _getCoordsAMinusCoordsB(currEndCoords, currStartCoords, multByFactor=convFactor/currDeltaT)

		#Add the current velocities to the traj steps
		trajSteps[idx-1].addExtraAttrDict( {velKey: {"value":currVels, "cmpType":"numericalArray"}} ) 


def _getCoordsAMinusCoordsB(coordsA, coordsB, multByFactor=1):
	assert len(coordsA)==len(coordsB)
	outCoords = np.zeros( [len(coordsA), 3] )

	for rIdx in range(len(outCoords)):
		for cIdx in range(3):
			outCoords[rIdx][cIdx] = (coordsA[rIdx][cIdx] - coordsB[rIdx][cIdx]) * multByFactor

	return outCoords.tolist()
	

def addAtomicTempsToTraj(inpTraj, atomTempKey="atomic_temps",velKey="velocities_from_pos", massDict=None):
	""" Calculates and attaches atomic temperatures for each step in inpTraj for which velKey is present (i.e. when velocities are already present on the trajectory)
	
	Args:
		inpTraj: (TrajectoryInMemory with TrajStepFlexible)
		velKey: (Str) The key under which atomic velocities are stored; see addVelocitiesToTrajInMemNVT for how to generate these.
		atomTempKey: (Str) The key we store atomic temperatures under
		massDict: (dict) Keys are element symbols, values are the masses to use. Default uses sensible values for elements. Values of mass are in atomic units (i.e. multiply by proton mass to get the absolute mass of an atom)

	Returns
		Nothing; works in place
 
	UNITS:
		Assumes that velocities are stored in ms^{-1} and aims to return temperatures in Kelvin

	"""
	for step in inpTraj:
		try:
			atomicTemps = _getAtomicTempArrayForTrajStep(step, velKey=velKey, massDict=massDict)
		except AttributeError:
			pass
		else:
			step.addExtraAttrDict( {atomTempKey: {"value":atomicTemps, "cmpType":"numericalArray"}} )


def _getAtomicTempArrayForTrajStep(trajStep, velKey="velocities_from_pos", massDict=None):
	#Sort out default args
	massDict = uCellHelp.getEleKeyToMassDictStandard() if massDict is None else massDict

	#Get the temperature values
	fractCoords = trajStep.unitCell.fractCoords
	velocities = getattr(trajStep, velKey)
	boltzConstant = uConvHelp.BOLTZMANN_CONSTANT_JOULE_PER_KELVIN * (1/uConvHelp.DALTON_TO_KG)
	atomicMasses = [massDict[ coords[-1] ] for coords in fractCoords] 
	velSqrd = [ vectHelp.getDotProductTwoVectors(vel,vel) for vel in velocities]
	atomicTemps = [(1/(3*boltzConstant))*mass*vSqrd for mass,vSqrd in it.zip_longest(atomicMasses,velSqrd)]

	return atomicTemps


def getTimeVsTempForTraj(inpTraj, inpIndices=None, atomTempKey="atomic_temps"):
	""" Extracts times vs temperature from an inpTraj which contains atomic temperatures in atomTempKey. These can be generated separately using addVelocitiesToTrajInMemNVT and addAtomicTempsToTraj
	
	Args:
		inpTraj: (TrajectoryInMemory with velKey in the relevant trajSteps)
		inpIndices: (iter of ints, Optional) Select indices to use to calculate temperature. Default is to use all atoms in the geometry
		atomTempKey: (str) The key in which to find the atomic temperatures
		massDict: (dict) Keys are element symbols, values are the masses to use. Default uses sensible values for elements. Values of mass are in atomic units (i.e. multiply by proton mass to get the absolute mass of an atom)
 
	Returns
		timesVsTemps: (iter of len-2 iters) [time,temp] in each element. 
 
	UNITS:
		Temperatures in Kelvin if input velocities are in metre per second.

	"""
	outTimes, outTemps = list(), list()
	for currStep in inpTraj:
		try:
			currTemps = getattr(currStep, atomTempKey)
		except AttributeError:
			pass
		else:
			currTime = currStep.time
			if inpIndices is None:
				currAverageTemp = sum(currTemps)/len(currTemps)
			else:
				tempsToAverage = [currTemps[idx] for idx in inpIndices]
				currAverageTemp = sum(tempsToAverage)/len(tempsToAverage)
			outTimes.append(currTime), outTemps.append(currAverageTemp)


	return [[t,val] for t,val in zip(outTimes,outTemps)]

