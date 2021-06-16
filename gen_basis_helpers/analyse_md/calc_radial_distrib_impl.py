

""" In theory a group of convenience functions for calculating various radial distribution functions """

import itertools as it

from . import binned_res as binResHelp
from . import calc_dists as calcDistHelp
from . import calc_distrib_core as calcDistribCoreHelp

from ..shared import cart_coord_utils as cartHelp
from ..shared import plane_equations as planeEqnHelp


def populatePlanarRdfsFromOptionsObjs(inpTraj, optionsObjs):
	""" Populates bin result objs sotred on CalcPlanarRdfOptions instances
	
	Args:
		inpTraj: (TrajectoryInMemory object)
		optionsObjs: (iter of CalcRdfOptions)
			 
	Returns
		Nothing; works in place on x.binResObj in optionsObjs
 
	"""
	binResObjs = [x.binResObj for x in optionsObjs]
	indices = [x.indices for x in optionsObjs]
	planeEqns = [x.planeEqn for x in optionsObjs]
	volumes = [ x.volume for x in optionsObjs]

	assert all([x==planeEqns[0] for x in planeEqns])

	_populateBinsWithPlanarRdfVals(inpTraj, binResObjs, indices, planeEqn=planeEqns[0], volumes=volumes)


class CalcPlanarRdfOptions():
	""" Object containing options to calculate a specific radial distribution function """
	
	def __init__(self, binResObj, indices, planeEqn=None, volume=None):
		""" Initializer
		
		Args:
			binResObj: (BinnedResultsStandard object) Note that this may get modified in place
			indices: (iter of ints) Contains the iter of indices to calculate the rdf for
			planeEqn: (None or ThreeDimPlaneEquation) The plane equation to calculate rdf from
			volume: (None or float) The volume to use for calculating the rdf. None generally means use the full unit cell volume (may not be sensible for slabs/multi-phase cells)
				 
		"""
		self.binResObj = binResObj
		self.indices = indices
		self.planeEqn = planeEqn
		self.volume = volume

	def getDistVsRdfData(self, offset=0):
		""" Returns [ [x1,rdf1], [x2,rdf2],... ] from self.binResObj; assuming its been populated
		
		Args:
			offset: (float, Optional) Optionally pass a value to be summed to the rdf data. Useful for making plots which data shifted along y
				 
		Returns
			outData: (iter of len-2 iters) x values are bin centres (calculated from the edges); y values are the rdf values
	 
		"""
		binEdgePairs = binResHelp.getBinEdgePairsFromBinResObj(self.binResObj)
		distVals = [ sum(edges)/len(edges) for edges in binEdgePairs ]
		rdfVals = [x+offset for x in self.binResObj.binVals["rdf"]]
		return [ [dist,rdf] for dist,rdf in it.zip_longest(distVals,rdfVals) ]

def _populateBinsWithPlanarRdfVals(inpTraj, binResObjs, indices, planeEqn=None, volumes=None):
	""" Gets planar rdf values for multiple binResObjs/atom groups separately. This can be used to efficiently calculate, for example, the effect of different bin sizing on the rdfs
	
	Args:
		inpTraj: (TrajectoryInMemory object)
		binResObjs: (iter of BinnedResultsStandard) These contain info on bin widths etc. and are where the results are added to
		indices: (iter of iter of ints) Each element corresponds to one binResObj. Within that each element contains the indices of atoms to get an rdf for
		volumes: (iter of floats) The total cell volume to assume for each case; Default is to use the unit cell volume. Using the whole cell may not be sensible when calculating for slab geometries.
		planeEqn: (a SINGLE ThreeDimPlaneEquation object) Default is to use the surface plane equation with d=0 (i.e. the bottom of the cell)

	NOTE:
		planeEqn MUST be parralel to axb; but this likely wont be checked

	Returns
		Nothing; works in place
 
	"""
	#Sort defaults
	nBins = len(binResObjs)
	volumes = calcDistribCoreHelp._getVolumesFromTrajAndInpVolumesArg(inpTraj, nBins, volumes)
	planeEqn = cartHelp.getABPlaneEqnWithNormVectorSameDirAsC_uCellInterface(inpTraj.trajSteps[0].unitCell) if planeEqn is None else planeEqn

	#Create the binner objects
	singleBinners = list()
	for binResObj, idxList in it.zip_longest(binResObjs, indices):
		currBinner = _PlanarRdfBinnerFixedIndices(resBins=binResObj, indices=idxList, planeEqn=planeEqn)
		singleBinners.append(currBinner)
	multiBinner = _PlanarRdfMultiBinnerFixedIndices( singleBinners )

	#Figure out the counts for each case
	nSteps = 0
	for trajStep in inpTraj:
		multiBinner.updateCountsFromTrajStep(trajStep)
		nSteps += 1

	#Attach the rdf info
	surfAreaAll = inpTraj.trajSteps[0].unitCell.volume / inpTraj.trajSteps[0].unitCell.lattParams["c"]
	for resObj, idxList, vol in it.zip_longest(binResObjs, indices, volumes):
		nIndices = len(idxList)
		_addPlanarRdfInfoToBinValsForBinsWithCounts(resObj, vol, nIndices, nSteps, surfAreaAll)


#Big duplication with _addRdfToBinValsForBinsWithCounts
def _addPlanarRdfInfoToBinValsForBinsWithCounts(binRes, volTotal, nIndices, nSteps, surfArea, countKey="counts"):
	allBinEdges = binResHelp.getBinEdgePairsFromBinResObj(binRes)
	outVols, outRdf = list(), list()
	prefactor = volTotal/nIndices
	for idx, currEdges in enumerate(allBinEdges):
		currCounts = binRes.binVals["counts"][idx]
		avCounts = currCounts / nSteps
		currWidth = currEdges[1]-currEdges[0]
		currVol = currWidth * surfArea
		currRdf = avCounts*prefactor/currVol

		outVols.append(currVol), outRdf.append(currRdf)

	binRes.binVals["volume"] = outVols
	binRes.binVals["rdf"] = outRdf


class _PlanarRdfMultiBinnerFixedIndices():

	def __init__(self, singleRdfBinners):
		self.singleBinners = singleRdfBinners

	def updateCountsFromTrajStep(self, trajStep):
		allIndices = [idx for idx in it.chain(*[binner.indices for binner in self.singleBinners])]
		uniqueIndices = sorted(set(allIndices))

		self._checkAllPlaneEquationsAreTheSame()
		usePlaneEqn = self.singleBinners[0].planeEqn
		allRelevantDists = calcDistHelp.calcDistancesFromSurfPlaneForCell(trajStep.unitCell, indices=uniqueIndices, planeEqn=usePlaneEqn)

		#Map these distances to an array which has a slot for each atom index
		useDists = list()
		maxIdx = max(uniqueIndices)
		idxInUniqueIndices = 0
		for idx in range(0,maxIdx+1):
			if idx in uniqueIndices:
				useDists.append( allRelevantDists[idxInUniqueIndices] )
				idxInUniqueIndices += 1
			else:
				useDists.append( 0 )

		#Bin relevant distances
		for binner in self.singleBinners:
			binner.updateCountsFromAllDists(useDists)


	def _checkAllPlaneEquationsAreTheSame(self):
		planeEqns = [x.planeEqn for x in self.singleBinners]
		for planeEqn in planeEqns:
			assert planeEqn==planeEqns[0]


class _PlanarRdfBinnerFixedIndices():

	def __init__(self, resBins=None, indices=None, planeEqn=None):
		""" Initializer
		
		Args:
			resBins: (BinnedResultsStandard)
			indices: (iter of ints)
			planeEqn: (ThreeDimPlaneEquation) Keeping as None will probably mean the surface plane equation is used (the one on the bottom of the cell)
				 
		"""
		self.resBins = resBins
		self.indices = indices
		self.planeEqn = planeEqn

	def updateCountsFromAllDists(self, allDists):
		valsToBin = [allDists[idx] for idx in self.indices]
		binResHelp.binCountsFromOneDimDataSimple(valsToBin, self.resBins)


