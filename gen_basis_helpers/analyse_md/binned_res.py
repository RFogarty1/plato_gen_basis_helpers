
import numpy as np

import itertools as it

class BinnedResultsStandard():
	""" Simple class for storing results of binning

	Attributes:
		binCentres: (iter of floats) Each value represents the centre of a bin
		binEdges: (iter of floats) Each value represents the edge of a bin (left->right). Same ordering as binCentres
		binVals: (dict) Keys describe a binned property (e.g. "rdf") while values are iters of floats, represnting the value of each bin

	"""

	def __init__(self, binCentres, binEdges, binVals):
		""" Initializer
		
		Args:
			binCentres: (iter of floats) Each value represents the centre of a bin
			binEdges: (iter of floats) Each value represents the edge of a bin (left->right). Same ordering as binCentres
			binVals: (dict) Keys describe a binned property (e.g. "rdf") while values are iters of floats, represnting the value of each bin
				 
		"""
		self._eqTol = 1e-5
		self.binCentres = binCentres
		self.binEdges = binEdges
		self.binVals = binVals

	@classmethod
	def fromConstantBinWidth(cls, binCentres, binWidth, binVals):
		binEdges = [x-(0.5*binWidth) for x in binCentres]
		binEdges.append( binCentres[-1] + (0.5*binWidth) )
		return cls(binCentres, binEdges, binVals)

	@classmethod
	def fromBinEdges(cls, binEdges, binVals=None):
		binVals = dict() if binVals is None else binVals
		outEdges = binEdges
		outCentres = list()
		for idx in range(len(binEdges)-1):
			currCentre = binEdges[idx] + ((binEdges[idx+1]-binEdges[idx])/2)
			outCentres.append(currCentre)
		return cls(outCentres, binEdges, binVals)

	def __eq__(self, other):
		eqTol = min(self._eqTol, other._eqTol)

		#Numerical attributes (but not the dictionary vals)
		numbAttrs = ["binCentres", "binEdges"]
		for attr in numbAttrs:
			valsA, valsB = getattr(self, attr), getattr(other,attr)
			if len( valsA ) != len( valsB ):
				return False
			diffVals = [abs(x-y) for x,y in zip(valsA,valsB)]
			if any([x>eqTol for x in diffVals]):
				return False

		#Check the binVals are equal
		binValsA, binValsB = getattr(self, "binVals"), getattr(other, "binVals")
		if binValsA.keys() != binValsB.keys():
			return False

		for key in binValsA.keys():
			valsA, valsB = binValsA[key], binValsB[key]
			if len(valsA) != len(valsB):
				return False

			try:
				unused = iter(valsA[0])
			#This will often trigger (e.g. if we're storing a single float value in each bin)
			except TypeError:
				diffVals = [abs(x-y) for x,y in zip(valsA,valsB)]
				if any([x>eqTol for x in diffVals]):
					return False

			#The else code runs if an exception wasnt triggered; this should trigger if each bin contains an iter of floats for one property
			else:
				#Check lengths are the same
				lenA = [len(x) for x in valsA]
				lenB = [len(x) for x in valsB]
				if lenA != lenB:
					return False

				#Check the values are equal
				for iterA, iterB in zip(valsA,valsB):
					diffVals = [abs(x-y) for x,y in zip(iterA,iterB)]
					if any([x>eqTol for x in diffVals]):
						return False

		return True


def binCountsFromOneDimDataSimple(inpData, binResObj, countKey="counts", raiseIfValsOutsideBins=False, initIfNeeded=True):
	""" Updates bins counts for
	
	Args:
		inpData: (iter of floats) Will update counts of bins accordingly
		binResObj: (BinnedResultsStandard)
		raiseIfValsOutsideBins: (Bool) If True raise a ValueError if we find a value which we cant bin
		initIfNeeded: (Bool) If True then it initialises the count key in the bins if its not already present as an initial step
 
	Returns
		 Nothing; works in place on binResObj. Values are placed/updated in binResObj.binVals[countKey]
 
	Raises:
		 ValueError: If "raiseIfValsOutsideBins" is True and a value is found that doesnt fit into any of the bins
	"""
	#Initialize the counts if needed
	if initIfNeeded:
		if binResObj.binVals.get(countKey,None) is None:
			binResObj.binVals[countKey] = [0 for x in range(len(binResObj.binCentres))]

	#Sort the inp data
	sortedData = sorted(inpData)

	#Sort the centres and widths
	binCentres = binResObj.binCentres
	binEdgePairs = [ [binResObj.binEdges[idx], binResObj.binEdges[idx+1]] for idx in range(len(binResObj.binCentres))]

	#Bin the data
	edgeIdx, dataIdx = 0,0
	numberPointsBinned = 0
	while (edgeIdx<len(binEdgePairs)) and (dataIdx<len(sortedData)):
		currMin, currMax = binEdgePairs[edgeIdx] 
		currVal = sortedData[dataIdx]

		if (currVal>=currMin) and (currVal<currMax):
			binResObj.binVals[countKey][edgeIdx] += 1
			numberPointsBinned += 1
			dataIdx += 1
		elif currVal>currMin:
			edgeIdx += 1
		else:
			dataIdx += 1 #This should only trigger when x is too small to fit into any bins
 

	if raiseIfValsOutsideBins:
		if numberPointsBinned < len(inpData):
			raise ValueError("Only {} values binned; out of {} present".format(numberPointsBinned, len(inpData) ))


def binResultsFromTwoDimDataSimple(inpData, binResObj, dimZeroLabel="xVals", dimOneLabel="yVals", raiseIfValsOutsideBins=True):
	""" When given x vs y, this divides the values into the bins provided in binResObj 
	
	Args:
		inpData: (iter of len-2 floats). x-values (i.e. inpData[idx][0]) determine which bin the data is placed in. y-values are some function we care about
		binResObj: (BinnedResultsStandard) Contains information on the bins
		dimZeroLabel: (str) The label to use
 
	Returns
		Nothing; works in place using binResObj. Values are placed in binResObj.binVals. If dimZeroLabel and dimOneLabel already exist then they are appended to the relevant lists; else the lists are created


	NOTES:
		a) For a bin with edges [minEdge,maxEdge] we will put x in that bin if minEdge<=x<maxEdge
		b) Overlapping bins likely wont cause any kind of error; but will lead to unpredictable results 
		c) Bins should be in order and continous

	Raises:
		 ValueError: If "raiseIfValsOutsideBins" is True and a value is found that doesnt fit into any of the bins
	"""
	#Sort the centres and widths
	binCentres = binResObj.binCentres
	sortedBinEdges = [ [binResObj.binEdges[idx], binResObj.binEdges[idx+1]] for idx in range(len(binResObj.binCentres))]

	#Sort the data
	sortedData = sorted(inpData)

	#Figure out the dict for the sorted bins
	outDicts = [ {dimZeroLabel:list(), dimOneLabel:list()} for x in range(len(sortedBinEdges)) ]

	edgeIdx, dataIdx = 0,0
	numberPointsBinned = 0
	while (edgeIdx<len(sortedBinEdges)) and (dataIdx<len(sortedData)):
		currMin, currMax = sortedBinEdges[edgeIdx]
		currDataX, currDataY = sortedData[dataIdx][0], sortedData[dataIdx][1]

		if (currDataX>=currMin) and (currDataX<currMax):
			outDicts[edgeIdx][dimZeroLabel].append( currDataX )
			outDicts[edgeIdx][dimOneLabel].append(  currDataY )
			numberPointsBinned += 1
			dataIdx += 1
		elif currDataX>currMin:
			edgeIdx += 1
		else:
			dataIdx += 1 #This should only trigger when x is too small to fit into any bins
 

	if raiseIfValsOutsideBins:
		if numberPointsBinned < len(sortedData):
			raise ValueError("Only {} values binned; out of {} present".format(numberPointsBinned, len(sortedData) ))

	#Initialise lists if needed
	if binResObj.binVals.get(dimZeroLabel,None) is None:
		binResObj.binVals[dimZeroLabel] = [list() for x in range(len(binCentres))]

	if binResObj.binVals.get(dimOneLabel,None) is None:
		binResObj.binVals[dimOneLabel] = [list() for x in range(len(binCentres))]

	#Update bins with new values
	for idx,unused in enumerate(binResObj.binCentres):
		extraXVals, extraYVals = outDicts[idx][dimZeroLabel], outDicts[idx][dimOneLabel]

		#Sort x-values out
		binResObj.binVals.get(dimZeroLabel)[idx].extend( outDicts[idx][dimZeroLabel] )

		#Sort y-values out
		binResObj.binVals.get(dimOneLabel)[idx].extend( outDicts[idx][dimOneLabel] )


def getEmptyBinResultsForValsStandard(inpVals, binWidth):
	""" Gets reasonable same-width bins for which to put inpVals in at a later point. Bins are generated such that the min value is in the centre of the first bin
	
	Args:
		inpVals: (iter of floats) Representative data which we will be binning
		binWidth: (float) The width of each bin

	Returns
		outBins: (BinnedResultsStandard object) This contains all the bins, but no data associated with them
 
	NOTE:
		Currently using np.arange as a backend; may have slight inconsistencies based on that

	"""
	minVal, maxVal = min(inpVals), max(inpVals)
	return getEmptyBinResultsFromMinMaxAndWidthStandard(minVal, maxVal, binWidth)


def getEmptyBinResultsFromMinMaxAndWidthStandard(minVal, maxVal, binWidth, extremesAtCentre=True):
	""" Gets reasonable same-width bins to put inpVals in at a later point. Bins are generated such that the bin value is in the centre of the first bin and the max value is in the centre of the final bin
	
	Args:
		minVal: (float) Expected minimum value
		maxVal: (float) Expected maximum value
		binWidth: (float) The width of the bins
		extremesAtCentre: (Bool) If True the we put minVal/maxVal at the centre of the bins (rather than the edges). 

	Returns
		outBins: (BinnedResultsStandard object) This contains all the bins, but no data associated with them
 
	NOTE:
		Currently using np.arange as a backend; may have slight inconsistencies based on that

	"""
	lowerEdge = minVal - (0.5*binWidth)
	upperEdge = maxVal + (0.5*binWidth)

	outEdges = list(np.arange(lowerEdge, upperEdge, binWidth))
	if outEdges[-1] < maxVal:
		outEdges.append( outEdges[-1]+binWidth )

	#Shift the edges to the extremes if requested
	if extremesAtCentre is False:
		nonCentredOutEdges = [x+(0.5*binWidth) for x in outEdges[:-1]]
		outEdges = nonCentredOutEdges

	return BinnedResultsStandard.fromBinEdges(outEdges)


def averageItersInEachBin(binResObj, keys=None):
	""" Takes a BinnedResultsStandard where iters of values are stored in each bin (meaninng binResObj.binVals.key is an iter of iters for that one)
	
	Args:
		binResObj: (BinnedResultsStandard) Objects containing results of binning
		keys: (iter of str) Keys in binResObj.binVals to attempt to average over. Default is to do it for all keys which are an iter of iters
			 
	Returns
		Nothing; works in place
 
	"""
	keys = binResObj.binVals.keys() if keys is None else keys

	for key in keys:
		currVals = binResObj.binVals[key]
		try:
			unused = iter(currVals[0])
		except TypeError:
			pass
		#Runs if error not triggered; in that case currVals is an iter of iters
		else:
			avVals = list()
			for valList in currVals:
				avVals.append( sum(valList)/len(valList) )
			binResObj.binVals[key] = avVals



def getBinEdgePairsFromBinResObj(binResObj):
	""" Gets [leftEdge, rightEdge] for each bin in binResObj
	
	Args:
		binResObj:(BinnedResultsStandard) Bins must be in order
			 
	Returns
		outEdges: (iter of len-2 floats) [leftEdge,rightEdge] for each bin in binResObj
 
	"""
	return [ [binResObj.binEdges[idx], binResObj.binEdges[idx+1]] for idx in range(len(binResObj.binCentres))]


#Backend checker for some distrib option classes
def _checkBinEdgesWithinDomain(binResObj, domain, domainTol):
	binEdges = binResObj.binEdges
	minEdge, maxEdge = min(binEdges), max(binEdges)
	if minEdge < domain[0]-abs(domainTol):
		raise ValueError("Bin with an edge of {} is outside domain of {}".format(minEdge, domain))
	if maxEdge > domain[1]+abs(domainTol):
		raise ValueError("Bin with an edge of {} is outside domain of {}".format(maxEdge, domain))



class NDimensionalBinnedResults():
	""" Class for storing results of binning procedures

	"""

	def __init__(self, edges, binVals=None):
		""" Initializer
		
		Args:
			edges: (iter of iter of floats) Each element is for one dimension. They should also be in order (though this wont be checked explicity)

		Returns
			What Function Returns
	 
		"""
		self.edges = edges
		self.binVals = binVals if binVals is not None else dict() #TODO: initialise? just need to make it an empty dict tbh


	#TODO: add some tolerance attrs
	def __eq__(self, other):

		#1) Check the edges 
		edgeArrayA, edgeArrayB = self.binEdgesArray, other.binEdgesArray
		if edgeArrayA.shape != edgeArrayB.shape:
			return False

		if not np.allclose(self.binEdgesArray, other.binEdgesArray):
			return False

		#2) Check the bin values
		keysA, keysB = self.binVals.keys(), other.binVals.keys()
		if keysA != keysB:
			return False

		for key in keysA:
			arrA, arrB = self.binVals[key], other.binVals[key]
			if not np.allclose(arrA,arrB):
				return False

		return True


	def _getBinEdgesList(self):
		""" Gets an iter with bin edges for each dimension. Can be (much) more useful than the array format """
		oneDimArrays = list()
		for currEdges in self.edges:
			currPairs = [ [currEdges[idx], currEdges[idx+1]] for idx in range(len(currEdges)-1)]
			oneDimArrays.append(currPairs)

		return oneDimArrays
	
	@property
	def binEdgesArray(self):
		""" Array containing all bin edges. For a 3-dim case binEdgesArray[idxA][idxB][idxC] will be a len-3 array [[minA,maxA], [minB,maxB], [minC,maxC]] which define edges for one bin in all dimensions """

		#Get edges for 1-dimensional case
		oneDimArrays = self._getBinEdgesList()

		#Get the all indices and edges
		outLens = [len(x) for x in oneDimArrays]
		oneDimIndices = list()
		for currLen in outLens:
			currIndices = [x for x in range(currLen)]
			oneDimIndices.append(currIndices)

		outIndices = [x for x in it.product(*oneDimIndices)]
		outVals = [x for x in it.product(*oneDimArrays)]

		#Combine all into an array
		outArray = np.zeros( outLens + [len(self.edges)] + [2] ) #The extra 2 is due to the min/max thing


		for indices, vals in it.zip_longest(outIndices, outVals):
			for valIdx in range(len(vals)):
				outArray[indices][valIdx] = vals[valIdx]

		return outArray
		

	@property
	def binCentresArray(self):
		""" Array contaning centres of each bin. For a 3-dim case binCentresArray[idxA][idxB][idxC] will be [centreA, centreB, centreC]"""

		edgesArray = self.binEdgesArray
		outShape = edgesArray.shape[:-1]
		outArray = np.zeros(outShape)
		allIndicesByDim = [ [x for x in range(currNumbEles)] for currNumbEles in outShape ]
		allIndices = [ x for x in it.product(*allIndicesByDim) ]

		for idxList in allIndices:
			currVal = (edgesArray[idxList][1]+edgesArray[idxList][0])/2
			outArray[idxList] = currVal

		return outArray

	def initialiseCountsMatrix(self, countKey="counts"):
		""" Initializes the matrix/tensor used to store counts to zero
		
		Args:
			countKey: (str) The matrix 
				 
		Returns
			Nothing; works in place
	 
		"""
		shapes = [len(x)-1 for x in self.edges]
		self.binVals[countKey] = np.zeros(shapes)


	#Probably gonna factor out a function that gets the indices to add to in order (so i could do more flexible things than +1)
	def addBinValuesToCounts(self, valsToBin, countKey="counts"):
		""" Updates bin values with counts
		
		Args:
			valsToBin: (iter of len-n iters) Where n is equal to the number of dimensions (length of self.edges). We add 1 to the counts for each bin
			countKey: (str) The key used for storing the counts matrix

		NOTE:
			Criteria to land in a bin with edges minEdge, maxEdge is minEdge<= x < maxEdge. Float imprecision probably mean this has little practical improtance though
	
		Returns
			Nothing; works in place
	 
		"""
		#Get pairs of edges for each bin (e.g. [[4,3],[3,2]]) . We sort each edge-pair into ascending order to make checks easier (e.g. [[4,3],[3,2]] becomes [[3,4],[2,3]])
		allEdgePairs = self._getBinEdgesList()
		sortedEdgePairs = list()
		for edgePairs in allEdgePairs:
			currPairs = [sorted(x) for x in edgePairs]
			sortedEdgePairs.append(currPairs)

		#Add each value
		for currValToBin in valsToBin:
			currIndices = [None for x in currValToBin]
			#For 3-dim case (for example) this will loop over [2,3,6] and find the right bin for 2,3,6 for the 1-dim case
			for idx,val in enumerate(currValToBin):
				for binIdx, (minEdge,maxEdge) in enumerate(sortedEdgePairs[idx]):
					if (minEdge <= val < maxEdge):
						currIndices[idx] = binIdx

			#Add a count in the relevant spot if a bin was found in EACH dimension
			if all([x is not None for x in currIndices]):
				self.binVals[countKey][tuple(currIndices)] += 1


