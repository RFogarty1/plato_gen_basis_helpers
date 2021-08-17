
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

