

import copy
import itertools as it

from . import calc_angular_distrib_impl as calcAngularDistribHelp
from . import calc_distrib_core as calcDistribCoreHelp
from . import calc_radial_distrib_impl as calcRadialDistribImplHelp

#Note: Only the three functions it calls are covereed by tests; not this function itself
def getAverageDistribFunctForEachInpTraj(inpOptions, inpTrajs):
	""" Calculates an averaged distribution function for each trajectory in inpTrajs. Best explained with an example. If we want to know the how a O-H bondlengths change over time (e.g. to check for dissociation events) we make an inpOption for EACH O-H bond; the output will then be an average of each OH bond (we calculate a distribution for each, and then average; this can NOT be replicated with any single distrib function).
	
	Args:
		inpOptions: (iter of CalcDistribOptionsBase objects) All from the same class (e.g. all CalcRdfOptions). The distrib. function (e.g. rdf) is calculated for each of these and the average is output for each trajectory
		inpTrajs: (iter of TrajectoryInMemory objects) 
 
	Returns
		 outVals: iter of iters. For each traj in inpTrajs we get an iter of len-2 iters (i.e an nx2 matrix) contaninig bin-centres (distances) vs the average value of the distribution function their (e.g. the average rdf value)

	Raises:
		 Errors
	"""
	allInpOpts = _getInpOptionsForEachTraj(inpOptions, inpTrajs)
	_populateInpOptsUsingTrajs(allInpOpts, inpTrajs)
	return _getAveragedDistribsFromPopulatedBinsMultiTraj(allInpOpts)

def _getInpOptionsForEachTraj(inpOptions, inpTrajs):
	outOpts = list()
	for traj in inpTrajs:
		currOpts = [copy.deepcopy(x) for x in inpOptions]
		outOpts.append(currOpts)
	return outOpts

def _populateInpOptsUsingTrajs(inpOptions, inpTrajs):
	firstOptObj = inpOptions[0][0]
	flattenedOptObjs = [x for x in it.chain(*inpOptions)]

	if isinstance(firstOptObj, calcDistribCoreHelp.CalcRdfOptions):
		assert all( [isinstance(x,calcDistribCoreHelp.CalcRdfOptions) for x in flattenedOptObjs] )
		for traj,opts in it.zip_longest(inpTrajs, inpOptions):
			calcDistribCoreHelp.populateRdfValsOnOptionObjs(traj, opts)

	elif isinstance(firstOptObj, calcRadialDistribImplHelp.CalcPlanarRdfOptions):
		assert all( [isinstance(x,calcRadialDistribImplHelp.CalcPlanarRdfOptions) for x in flattenedOptObjs] )
		for traj,opts in it.zip_longest(inpTrajs, inpOptions):
			calcRadialDistribImplHelp.populatePlanarRdfsFromOptionsObjs(traj, opts)

	elif isinstance(firstOptObj, calcAngularDistribHelp.CalcInteratomicAngularDistribOptions):
		assert all( [isinstance(x,calcAngularDistribHelp.CalcInteratomicAngularDistribOptions) for x in flattenedOptObjs] )
		for traj,opts in it.zip_longest(inpTrajs, inpOptions):
			calcAngularDistribHelp.populateInteratomicAngularDistribsFromOptionsObjs(traj, opts)

	else:
		raise NotImplementedError("")

def _getAveragedDistribsFromPopulatedBinsMultiTraj(inpOptions):
	averagedRdfData = list()
	for trajIdx, optsIter in enumerate(inpOptions):
		currDists = [x[0] for x in optsIter[0].getDistVsDistribData()]
		trajRdfs = list()
		for currOpts in optsIter:
			currRdf = [x[1] for x in currOpts.getDistVsDistribData()]
			trajRdfs.append(currRdf)

		currAverages = [sum(x)/len(x) for x in it.zip_longest(*trajRdfs)]
		averagedRdfData.append( [ [dist,rdf] for dist,rdf in it.zip_longest(currDists,currAverages) ] )

	return averagedRdfData




