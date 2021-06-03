
import numpy as np

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


	#Figure out the dict for the sorted bins
	outDicts = [ {dimZeroLabel:list(), dimOneLabel:list()} for x in range(len(sortedBinEdges)) ]

	edgeIdx, dataIdx = 0,0
	numberPointsBinned = 0
	while (edgeIdx<len(sortedBinEdges)) and (dataIdx<len(inpData)):
		currMin, currMax = sortedBinEdges[edgeIdx]
		currDataX, currDataY = inpData[dataIdx][0], inpData[dataIdx][1]

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
		if numberPointsBinned < len(inpData):
			raise ValueError("Only {} values binned; out of {} present".format(numberPointsBinned, len(inpData) ))

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
	""" Gets reasonable same-width bins for which to put inpVals in at a later point. Bins are generated such that the min value in the centre of the first bin
	
	Args:
		inpVals: (iter of floats) Representative data which we will be binning
		binWidth: (float) The width of each bin

	Returns
		outBins: (BinnedResultsStandard object) This contains all the bins, but no data associated with them
 
	NOTE:
		Currently using np.arange as a backend; may have slight inconsistencies based on that

	"""
	minVal, maxVal = min(inpVals), max(inpVals)
	lowerEdge = minVal - (0.5*binWidth)
	upperEdge = maxVal + (0.5*binWidth)

	outEdges = list(np.arange(lowerEdge, upperEdge, binWidth))
	if outEdges[-1] < maxVal:
		outEdges.append( outEdges[-1]+binWidth )

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



