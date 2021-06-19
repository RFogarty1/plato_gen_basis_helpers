
import collections
import itertools as it
import math

from . import binned_res as binResHelp
from . import calc_dists as calcDistsHelp

""" Module to provide core (generally backend) functions to help calculation of radial/angular distributions """


def populateRdfValsOnOptionObjs(inpTraj, optionsObjs):
	""" Populates bin results objs stored on CalcRdfOptions instances
	
	Args:
		inpTraj: (TrajectoryInMemory object)
		optionsObjs: (iter of CalcRdfOptions)
			 
	Returns
		Nothing; works in place on x.binResObj in optionsObjs
 
	"""
	#
	binResObjs = [x.binResObj for x in optionsObjs]
	indicesA = [x.indicesA for x in optionsObjs]
	indicesB = [x.indicesB for x in optionsObjs]
	volumes = [x.volume for x in optionsObjs]

	#
	_populateBinsWithRdfBetweenAtomGroups(inpTraj, binResObjs, indicesA, indicesB, volumes=volumes)




class CalcDistribOptionsBase():
	""" Base class for Objects representing options for calculating distribution functions (e.g. radial distribution functions) """


	def getDistVsDistribData(self, offset=0):
		""" Returns [ [x1,rdf1], [x2,rdf2],... ] from self.binResObj; assuming its been populated
		
		Args:
			offset: (float, Optional) Optionally pass a value to be summed to the rdf data. Useful for making plots which data shifted along y
				 
		Returns
			outData: (iter of len-2 iters) x values are bin centres (calculated from the edges); y values are the rdf values
	 
		"""
		binEdgePairs = binResHelp.getBinEdgePairsFromBinResObj(self.binResObj)
		distVals = [ sum(edges)/len(edges) for edges in binEdgePairs ]
		rdfVals = [x+offset for x in self.binResObj.binVals.get(self.distribKey)]
		return [ [dist,rdf] for dist,rdf in it.zip_longest(distVals,rdfVals) ]


class CalcRdfOptions(CalcDistribOptionsBase):
	""" Object containing options to calculate a specific rdf function """

	def __init__(self, binResObj, indicesA, indicesB, volume=None):
		""" Initializer
		
		Args:
			binResObj: (BinnedResultsStandard object) Note that this may get modified in place
			indicesA: (iter of ints) Contains the indices of atoms to get an rdf FROM (e.g. for g_{AB} indicesA should contain all indices of atom type A)
			indicesB: (iter of ints) Contains the indices of atoms to get an rdf TO (e.g. for g_{AB} indicesB should contain all indices of atom type B)
			volume: (float or None) The total cell volume to assume; Default is to use the unit cell volume from first step of an input traj. Using the whole cell may not be sensible when calculating for slab geometries.

		"""
		self.distribKey = "rdf"
		self.binResObj = binResObj
		self.indicesA = indicesA
		self.indicesB = indicesB
		self.volume = volume


#TODO: Probably introduce some command objects to determine the options for these; so i can individually specify the runs but get them all combined this way
def _populateBinsWithRdfBetweenAtomGroups(inpTraj, binResObjs, indicesA, indicesB, volumes=None):
	""" Gets rdf functions for multiple binResObjs/atom groups simultaneously. This can be used to efficiently calculate multiple rdf for different element combinations (e.g. g_{OH}/g_{OO}) or the same element combo with varying bin widths
	
	Args:
		inpTraj: (TrajectoryInMemory object)
		binResObjs: (iter of BinnedResultsStandard) These contain info on bin widths etc. and are where the results are added to
		indicesA: (iter of iter of ints) Each element corresponds to one binResObj. Within that each element contains the indices of atoms to get an rdf FROM (e.g. for g_{AB} indicesA should contain all indices of atom type A)
		indicesB: (iter of iter of ints) Each element corresponds to one binResObj. Within that each element contains the indices of atoms to get an rdf TO (e.g. for g_{AB} indicesB should contain all indices of atom type B)
		volumes: (iter of floats) The total cell volume to assume for each case; Default is to use the unit cell volume. Using the whole cell may not be sensible when calculating for slab geometries.
 
	Returns
		Nothing; works in place
 
	"""
	#Sort annoying defaults
	nBins = len(binResObjs)
	volumes = _getVolumesFromTrajAndInpVolumesArg(inpTraj, nBins, volumes)

	#Create objects to handle binning for each traj step
	singleBinners = list()
	for resObj, idxListA, idxListB in it.zip_longest(binResObjs, indicesA, indicesB):
		currBinner = _RdfBinnerFixedIndices(resObj, idxListA, idxListB)
		singleBinners.append(currBinner)
	multiBinner = _MultiRdfBinnerFixedIndices(singleBinners)

	#Figure out the counts for each case(and append to the results object)
	nSteps = 0
	for trajStep in inpTraj:
		multiBinner.updateCountsFromTrajStep(trajStep)
		nSteps += 1

	#Attach the rdf
	for resObj, idxListA, idxListB,vol in it.zip_longest(binResObjs, indicesA, indicesB, volumes):
		nA, nB = len(idxListA), len(idxListB)
		_addRdfToBinValsForBinsWithCounts(resObj, vol, nA, nB, nSteps)


def _getVolumesFromTrajAndInpVolumesArg(inpTraj, nBins, volumes):
	volInit = inpTraj.trajSteps[0].unitCell.volume
	if volumes is None:
		volumes = [volInit for x in range(nBins)]
	else:
		for idx,vol in enumerate(volumes):
			if vol is None:
				volumes[idx] = volInit
	return volumes

#SOOOOOOO: We need to know both number of atoms in group A and group B to get the normalisation factor
def _addRdfToBinValsForBinsWithCounts(binRes, volTotal, nA, nB, nSteps, countKey="counts"):
	""" Adds rdf values for each bin assuming count values are present
	
	Args:
		binRes: (BinnedResultsStandard)
		volTotal: (float) The total volume of the system
		nA: (int) The number of atoms in group A for g_{AB} (this equals 1 if we're doing the rdf from a single point)
		nB: (int) The number of atoms in group B for g_{AB}
		nSteps: (int) Number of steps done
 		countKey: (str) Lets the program know where to find the "counts" in the binRes.binVals dictionary

	Returns
		 Nothing; works in place. But adds "volume" and "rdf" to binVals;volume is the volume of the bin in 3-dimensional space while rdf is obviously the radial distribution function
 
	Notes: g_{AB} (bin) = (counts(bin)/nSteps) * V_{tot} / (nA*nB*V_{bin})

	"""

	allBinEdges = binResHelp.getBinEdgePairsFromBinResObj(binRes)

	outVols, outRdf = list(), list()
	prefactor = volTotal/(nA*nB)
	for idx, currEdges in enumerate(allBinEdges):
		currCounts = binRes.binVals["counts"][idx]
		avCounts = currCounts / nSteps
		currWidth = currEdges[1]-currEdges[0]
		currCentre = currEdges[0] + (0.5*currWidth)
		currVol = 4*math.pi*(currCentre**2)*currWidth
		currRdf = avCounts*prefactor/currVol

		outVols.append(currVol), outRdf.append(currRdf)

	binRes.binVals["volume"] = outVols
	binRes.binVals["rdf"] = outRdf
	return binRes


class _MultiRdfBinnerFixedIndices():

	def __init__(self, singleRdfBinners):
		self.singleBinners = singleRdfBinners

	def updateCountsFromTrajStep(self, trajStep):
		fullDistMatrix = calcDistsHelp.calcDistanceMatrixForCell_minImageConv(trajStep.unitCell)
		for currBinner in self.singleBinners:
			currBinner.updateCountsFromDistMatrix(fullDistMatrix)


class _RdfBinnerFixedIndices():

	def __init__(self, resBins=None, indicesA=None, indicesB=None):
		""" Initializer
		
		Args:
			resBins: (BinnedResultsStandard)
			indicesA: (iter of ints)
			indicesB: (iter of ints)

		"""
		self.resBins = resBins
		self.indicesA = indicesA
		self.indicesB = indicesB

	def updateCountsFromTrajStep(self, trajStep):
		raise NotImplementedError("")

	def updateCountsFromDistMatrix(self, distMatrix):
		valsToBin = _getRadialToBinValsFromFullDistMatrix(distMatrix, indicesA=self.indicesA, indicesB=self.indicesB)
		binResHelp.binCountsFromOneDimDataSimple(valsToBin, self.resBins)

#May be able to actually make this more general than distances
def _getRadialToBinValsFromFullDistMatrix(distMatrix, indicesA=None, indicesB=None):
	""" Takes a distance matrix and extracts relevant distances to bin based on indicesA and indicesB
	
	Args:
		distMatrix: (nxm matrix) distMatrix[idxA][idxB] 
		indicesA: (iter of ints) Indices of atoms in rows to bin. Default is to use ALL
		indicesB: (iter of ints) Indices of atoms in columns to bin. Default is to use ALL
 
	Returns
		outVals: (iter of floats; len(indicesA)*len(indicesB)) Values which we will be binning later.

	Raises:
		AssertionError: If distMatrix is not square. Tests number of rows and then columns in first row to determine this.

	NOTE:
		a) The distance matrix needs to be symmetrical; this will generally mean we want the distance matrix for all atoms in the cell with each other
		b) If idxA appears in indicesA and indicesB, distMatrix[idxA][idxA] is ignored 

	"""
	indicesA = [idx for idx,unused in enumerate(distMatrix)] if indicesA is None else indicesA
	indicesB = [idx for idx,unused in enumerate(distMatrix[0])] if indicesB is None else indicesB

	#Check the matrix is square
	nRows, nCols = len(distMatrix), len(distMatrix[0])
	assert nRows==nCols, "Need a square matrix but found nRows={}, nCols={}".format(nRows,nCols)

	#Figure out all pairs of indices
	idxPairs = [ [x,y] for x,y in it.product(indicesA,indicesB) ]

	#Filter out equivalent ones ( [2,1]==[1,2] )
 	#Taken from https://stackoverflow.com/questions/38187286/find-unique-pairs-in-list-of-pairs
	ctr = collections.Counter(frozenset(x) for x in idxPairs)
	bools = [ctr[frozenset(x)]==1 for x in idxPairs]

	uniquePairs = list()
	for idx,pair in enumerate(idxPairs):
		if bools[idx]:
			uniquePairs.append(pair)
		elif  pair[0]>pair[1]:
			uniquePairs.append(pair)
		else:
			pass

	#Get the distances to bin
	outVals = list()
	for currPair in uniquePairs:
		idxA, idxB = currPair
		if idxA != idxB:
			outVals.append( distMatrix[idxA][idxB] )

	return outVals

