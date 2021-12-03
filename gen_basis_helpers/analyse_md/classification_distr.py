
from . import atom_combo_opts_obj_maps as atomComboOptObjMapHelp

#TODO: I'll probably want to use and iter-of-iter of optsObj interface soon enough....
#Since it could save tons of time. Though initially this "single-obj" method is more convenient + can easily be refactored to use
#more general interface as a backend
def getClassDistrCountsOverTraj(inpTraj, optsObj):
	""" Gets the number of X for each classification for each step in inpTraj
	
	Args:
		inpTraj: (TrajectoryInMemory) (Actually probably also works for TrajectoryBase)
		optsObj: (Generally CalcDistribOptionsBase obj) Really any options object leading to a binval getters that return len-1 iters for each traj step
		
	Returns
		outCounts: (NxM iter of iters) Contains counts for each type at each traj step. N is the number of traj-steps, M is the number of possible classifications
 
	"""
	sparseMatrixCalculator = atomComboOptObjMapHelp.getSparseMatrixCalculatorFromOptsObjIter([optsObj])
	binValGetter = atomComboOptObjMapHelp.getMultiDimBinValGetterFromOptsObjs([optsObj])
	outList = list()

	for trajStep in inpTraj:
		sparseMatrixCalculator.calcMatricesForGeom(trajStep.unitCell)
		currVals = binValGetter.getValsToBin(sparseMatrixCalculator)
		assert len(currVals)==1
		outList.append( currVals[0] )

	return outList



def getFilteredScatterDataOverTime(inpTraj, optsObj, foldOneDimData=True):
	""" Gets data [ [timeA,valsA], [timeA,valsB], [timeB,valsC] ] for filterered groups. Original goal was to track planar positions of self-ionised water groups (hydroxyl/hydronium) which spontaneously appeared in an MD run
	
	Args:
		inpTraj: (TrajectoryInMemory) (Actually probably also works for TrajectoryBase)
		optsObj: Provides options for what values we want. In original example we use "WaterDerivativeFilteredOptsObj_simple"
		foldOneDimData: (Bool, Optional) If True, for one-dimensional output data (e.g. just finding planar position of filtered) instead of returning an iter of [time, [value]] we return [time,value], where value is a float of interest

	Returns
		outVals: Returns [time,values] pair for each. Note: individual time values appear N times, where N is the number of relevant groups at that time (e.g. timestep=3000fs may appear 5 times if 5 hydroxyl groups are found in that one, in the case we want to track hydroxyls) 
 
	"""
	sparseMatrixCalculator = atomComboOptObjMapHelp.getSparseMatrixCalculatorFromOptsObjIter([optsObj])
	binValGetter = atomComboOptObjMapHelp.getMultiDimBinValGetterFromOptsObjs([optsObj])
	outList = list()

	for trajStep in inpTraj:
		sparseMatrixCalculator.calcMatricesForGeom(trajStep.unitCell)
		currTime = trajStep.time
		currBinnedVals = binValGetter.getValsToBin(sparseMatrixCalculator)
		currVals = [ [currTime, val] for val in currBinnedVals ]
		outList.extend(currVals)

	if foldOneDimData:
		if len(outList)>0:
			if len(outList[0][1])==1:
				outList = [ [t,val[0]] for t,val in outList ]



	return outList



