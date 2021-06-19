
import itertools as it

from . import binned_res as binResHelp
from . import calc_dists as calcDistsHelp
from . import calc_distrib_core as calcDistribCoreHelp


def populateInteratomicAngularDistribsFromOptionsObjs(inpTraj, optionsObjs):
	""" Populates bin result objs stored on CalcInteratomicAngularDistribOptions instances
	
	Args:
		inpTraj: (TrajectoryInMemory object)
		optionsObjs: (iter of CalcInteratomicAngularDistribOptions objects)
			 
	Returns
		Nothing; works in place on x.binResObj in optionsObjs
 
	"""
	binResObjs = [x.binResObj for x in optionsObjs]
	indices = [x.indices for x in optionsObjs]
	_populateBinsWithInteratomicAngularDistribs(inpTraj, binResObjs, indices)

#TODO: Have some functions checking the bin edges are within a tolerance of the domain
class CalcInteratomicAngularDistribOptions(calcDistribCoreHelp.CalcDistribOptionsBase):
	

	def __init__(self, binResObj, indices, checkEdges=True):
		""" Initializer
		
		Args:
			binResObj: (BinnedResultsStandard object) Note that this may get modified in place
			indices: (iter of len-3 iters of ints) Each element is 3 indices representing an angle to be calculated
			checkEdges: (Bool) If True check that the edges of bins fall within the domain (0 to 180 degrees)
 
		"""
		self.domainTol = 1e-1
		self.domain = [0,180]
		self.distribKey = "adf" #Angular distribution function
		self.binResObj = binResObj
		self.indices = indices
		if checkEdges:
			self._checkBinEdgesWithinDomain()

	def _checkBinEdgesWithinDomain(self):
		binEdges = self.binResObj.binEdges
		minEdge, maxEdge = min(binEdges), max(binEdges)
		if minEdge < self.domain[0]-abs(self.domainTol):
			raise ValueError("Bin with an edge of {} is outside domain of {}".format(minEdge, self.domain))
		if maxEdge > self.domain[1]+abs(self.domainTol):
			raise ValueError("Bin with an edge of {} is outside domain of {}".format(maxEdge, self.domain))



def _populateBinsWithInteratomicAngularDistribs(inpTraj, binResObjs, indices):
	""" Populates angular distribution functions (both probability and g(theta)) for multiple binResObjs/angular indices separately. This can be useful, for example, for quickly checking effects of varying bin widths

	Args:
		inpTraj: (TrajectoryInMemory object)
		binResObjs: (iter of BinnedResultsStandard) These contain info on bin widths etc. and are where the results are added to
		indices: (iter of iter of iter of len-3 iters) Eac element correponds to one binResObj. Within that each element contains an iter of len-3 iters; each of those contains 3 atom indices representing one angle to calculate
			 
	Returns
		Nothing; works in place
 
	"""
	#Create binners
	singleBinners = list()
	for binRes, idxList in it.zip_longest(binResObjs, indices):
		currSingleBinner = 	_InteratomicAngularBinnerFixedIndices(binRes, idxList)
		singleBinners.append(currSingleBinner)
	multiBinner = _InteratomicAngularMultiBinnerFixedIndices(singleBinners)

	#Run binners
	for tStep in inpTraj.trajSteps:
		multiBinner.updateCountsFromTrajStep(tStep)

	#Calculate the pdf and adf from the counts and number of steps
	for binResObj,idxList in it.zip_longest(binResObjs,indices):
		nAngles = len(idxList)
		_addAngularPdfAndRdfToBinObj(binResObj, nAngles)

#nSteps maybe not even needed actually
def _addAngularPdfAndRdfToBinObj(inpBinObj, nAngles):
	domain = [0,180]
	domainWidth = abs(domain[1]-domain[0])

	binEdgePairs = binResHelp.getBinEdgePairsFromBinResObj(inpBinObj)
	binWidths = [abs(x[1]-x[0]) for x in binEdgePairs]
	
	#Get pdf (probability distribution function) and adf (angular distribution function)
	totalAngles = inpBinObj.binVals
	outPdfs, outAdfs = list(), list()
	for currWidth, currCount in it.zip_longest(binWidths, inpBinObj.binVals["counts"]):
		currPdf = currCount / nAngles
		currAdf = currPdf * (domainWidth/currWidth)
		outPdfs.append(currPdf), outAdfs.append(currAdf)

	#Attach to the bins
	inpBinObj.binVals["pdf"] = outPdfs
	inpBinObj.binVals["adf"] = outAdfs


class _InteratomicAngularMultiBinnerFixedIndices():

	def __init__(self, singleBinners):
		self.singleBinners = singleBinners

	#TODO: Cache this so that we only calculate the annoying mapping stuff once
	def updateCountsFromTrajStep(self, trajStep):
		#Get all angle indices we need[Remove to cache function]
		allAngleLists =[idxList for idxList in it.chain(*[x.indices for x in self.singleBinners])]
		anglesIndicesNoDupl = [x for x in set([tuple(x) for x in allAngleLists])]
		anglesIndicesNoDupl = [list(x) for x in anglesIndicesNoDupl]

		#[remove to cache function]
		#Figure out how to map the angle lists positions to indices in each single binner
		reqdAnglesIndicesEachBinner = [list() for x in self.singleBinners]
		for binnerIdx, binnerObj in enumerate(self.singleBinners):
			for angleIdx in binnerObj.indices:
				for noDuplSingleIdx, noDuplAngleIdx in enumerate(anglesIndicesNoDupl):
					if noDuplAngleIdx == angleIdx:
						reqdAnglesIndicesEachBinner[binnerIdx].append(noDuplSingleIdx)

		#Calculate the angles for all
		allAngles = calcDistsHelp.getInterAtomicAnglesForInpGeom(trajStep.unitCell, anglesIndicesNoDupl)

		#Pass the RELEVANT calculated indices to each bin
		for binnerIdx, binnerObj in enumerate(self.singleBinners):
			currAngles = [ allAngles[idx] for idx in reqdAnglesIndicesEachBinner[binnerIdx] ]
			binnerObj.updateCountsFromCalcdAngles(currAngles)


class _InteratomicAngularBinnerFixedIndices():

	def __init__(self, resBins=None ,indices=None):
		""" Initializer
		
		Args:
			resBins: (BinnedResultsStandard)
			indices: (iter of len-3 iters) Each element is for one angle (e.g. [4,10,20] means the angle between atom indices 4,10,20)
				 
		Returns
			What Function Returns
	 
		Raises:
			Errors
		"""
		self.resBins = resBins
		self.indices = indices

	def updateCountsFromCalcdAngles(self, calcdAngles):
		binResHelp.binCountsFromOneDimDataSimple(calcdAngles, self.resBins)


