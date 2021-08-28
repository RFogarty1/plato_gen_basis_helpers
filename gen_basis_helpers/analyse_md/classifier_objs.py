
import numpy as np

from . import atom_combo_binval_getters as atomComboBinvalGetterHelp
from . import atom_combo_core as atomComboCoreHelp

#Note some implementations are currently in classification_binval_getters; since they were originally JUST
#backends for counting different types of atoms/molecules


class ClassifierBase():
	""" Object which returns a group of indices when given a sparseMatrixCalculator (with relevant matrices populated) """

	def classify(self, sparseMatrixCalculator):
		""" 
		
		Args:
			sparseMatrixCalculator: (_SparseMatrixCalculatorStandard) This contains various matrices (e.g. a distance matrix)
				 
		Returns
			outIndices: Generally iter of ints for the relevant group. May be iter of these; e.g. for water we will return both oxyIndices and hyIndices as separate 
	 
		"""
		raise NotImplementedError("")


def getByReferenceClassifiers(inpClassifiers, startExecCount=0):
	""" Gets "classifiers" that work by reference to the inpClassifiers. By reference means that when .classify is called, they take the values from the "reference" classifiers instead of figuring out values themselves. This is all essentially an optimization trick/hack to reduce amount of duplicate work, without some overarching co-ordinator
	
	Args:
		inpClassifiers: (iter of ClassifierBase) objects
		startExecCount: (int) The initial execution count. We track the number of times classify() is called and make sure it stays in sync with the number of calls made to the reference one [May have option to disable this check later; but for now its a useful failsafe]
 
	Returns
		byRefClassifiers: ()
 
	"""
	outVals = [_ByReferenceClassifier(inpClassifier, execCount=startExecCount) for inpClassifier in inpClassifiers]
	return outVals


class _ByReferenceClassifier(ClassifierBase):
	""" Classifies by looking at previous results of another classifier; used as an optimisation trick to not repeat same (expensive) classification operations """
	
	def __init__(self, refClassifier, execCount=0):
		self.refClassifier = refClassifier
		self.execCount = execCount

	def classify(self, sparseMatrixCalculator):
		#Check exec counts
		if self.refClassifier.execCount-1 == self.execCount:
			pass
		else:
			currArgs = [self.refClassifier.execCount, self.execCount]
			raise ValueError("Reference exec count = {}, but current object has count = {}; ref count should be 1 higher".format(*currArgs))

		outVals = self.refClassifier.storedClassifyResult
		self.execCount += 1

		return outVals


#Main difference is it returns (oxyIndices,hyIndices) rom classify; should maybe put that in doc-string somewhere
class _WaterClassifierBase():
	""" Classifies list of water indices into distinct groups of water """
	

	def classify(self, sparseMatrixCalculator):
		raise NotImplementedError("")



#Specific implementations of classifiers here
class _AtomsWithinMinDistRangeClassifier():

	def __init__(self, atomIndices, distFilterIndices, distFilterRange, minDistVal=-0.01, execCount=0):
		""" Initializer
		
		Args:
			atomIndices: (iter of ints)
			distFilterIndices: (iter of ints) Each represents an atom index. We group atomIndices by min-distance from these indices
			distFilterRange: (len-2 float iter) [minDist,maxDist] from indices in distFilterIndices for an atom to be included in the counts
			minDistVal: (float) If set to a +ve number we ignore distances smaller than it when figuring out minimum. Useful to avoid getting zeros when atomIndices and distFilterIndices overlap
			execCount: (int) Used to track how many times .classify is called; used as a safety check when using "byReference" classifiers

		"""
		self.atomIndices = atomIndices
		self.distFilterIndices = distFilterIndices
		self.distFilterRange = distFilterRange
		self.minDistVal = minDistVal
		self.execCount = 0

	def classify(self, sparseMatrixCalculator):
		#1) Create binner for minimum distances
		currArgs = [self.atomIndices, self.distFilterIndices]
		minDistBinner = atomComboBinvalGetterHelp._MinDistsGetOneDimValsToBin(*currArgs,minVal=self.minDistVal)

		#2) Get the relevant values for all distances
		minDists = minDistBinner.getValsToBin(sparseMatrixCalculator)

		#3) Get the indices which fall between distFilterRange limits
		outIndices = list()
		useFilterRange = sorted(self.distFilterRange)
		for idx, minDist in enumerate(minDists):
			if (minDist>=useFilterRange[0]) and (minDist<useFilterRange[1]):
				outIndices.append( self.atomIndices[idx] )

		self.storedClassifyResult = outIndices
		self.execCount += 1

		return outIndices




class _WaterClassifierMinDistAndNumberHBonds(_WaterClassifierBase):

	def __init__(self, oxyIndices, hyIndices, distFilterIndices, distFilterRange, nDonorFilterRange,
	             nAcceptorFilterRange, nTotalFilterRange, maxOOHBond, maxAngleHBond):
		""" Initializer
		
		Args:
			oxyIndices: (iter of ints) The oxygen indices for each water molecule
			hyIndices: (iter of len-2 ints) Same length as oxyIndices, but each contains the indices of two hydrogen indices bonded to the relevant oxygen
			distFilterIndices: (iter of ints) Each represents an atom index. We group water by distance of oxygen atoms from these indices
			distFilterRange: (len-2 float iter) [minDist,maxDist] from indices in distFilterIndices for a water to be included in this count
			nDonorFilterRange: (len-2 float iter) [minNDonor, maxNDonor] for a water.
			nAcceptorFilterRange: (len-2 float iter) [minNTotal,maxNTotal] for a water
			nTotalFilterRange: (len-2 float iter) [minNTotal,maxNTotal] for a water
			maxOOHBond: The maximum O-O distance between two hydrogen-bonded water.
			maxAngleHBond:  The maximum OA-OD-HD angle for a hydrogen bond; OA = acceptor oxygen, OD=Donor oxygen, HD=donor hydrogen
	 
		"""
		self.oxyIndices = oxyIndices
		self.hyIndices = hyIndices
		self.distFilterIndices = distFilterIndices
		self.distFilterRange = distFilterRange
		self.nDonorFilterRange = nDonorFilterRange
		self.nAcceptorFilterRange = nAcceptorFilterRange
		self.nTotalFilterRange = nTotalFilterRange
		self.maxOOHBond = maxOOHBond
		self.maxAngleHBond = maxAngleHBond

	def classify(self, sparseMatrixCalculator):
		#Step 0) Create a bunch of other binner objects to get values we need
		multiDimBinner = self._getRelevantMultiDimBinner()

		#Step 1): Get all the actual values of everything we filter by
		relValsAll = multiDimBinner.getValsToBin(sparseMatrixCalculator)

		#Step 2) Get the indices of each that are in one group
		outOxyIndices = list()
		outHyIndices = list()

		for oxyIdx,(minDist, nDonor, nAcceptor, nTotal) in enumerate(relValsAll):
			if (self.distFilterRange[0]<=minDist) and (minDist<self.distFilterRange[1]):
				if (self.nDonorFilterRange[0]<=nDonor) and (nDonor<self.nDonorFilterRange[1]):
					if (self.nAcceptorFilterRange[0]<=nAcceptor) and (nAcceptor<self.nAcceptorFilterRange[1]):
						if (self.nTotalFilterRange[0]<=nTotal) and (nTotal<self.nTotalFilterRange[1]):
							outOxyIndices.append( self.oxyIndices[oxyIdx] )
							outHyIndices.append( self.hyIndices[oxyIdx] )

		return outOxyIndices, outHyIndices


	def _getRelevantMultiDimBinner(self):
		#Deal with distances
		minDistType = "O"
		currArgs = [self.oxyIndices, self.hyIndices, self.distFilterIndices, minDistType]
		minDistBinner = atomComboBinvalGetterHelp._WaterMinDistBinValGetter(*currArgs)

		#Deal with number of donors (overpowered binner but...)
		dudFilterIndices, dudDistFilterVals = self.oxyIndices, [ [-0.1,10], [-0.1,10] ] #Should be equivalent to not filtering between groupA/groupB 
		currArgs = [self.oxyIndices, self.hyIndices, dudFilterIndices, dudDistFilterVals]
		currKwargs = {"acceptor":False, "donor":True, "maxOO":self.maxOOHBond, "maxAngle":self.maxAngleHBond}
		nDonorBinner = atomComboBinvalGetterHelp._DiscHBondCounterBetweenGroupsWithOxyDistFilterOneDimValGetter(*currArgs, **currKwargs)

		#Deal with number of acceptors
		currKwargs = {"acceptor":True, "donor":False, "maxOO":self.maxOOHBond, "maxAngle":self.maxAngleHBond}
		nAcceptorBinner = atomComboBinvalGetterHelp._DiscHBondCounterBetweenGroupsWithOxyDistFilterOneDimValGetter(*currArgs, **currKwargs)

		#Deal with total number of h-bonds
		currKwargs = {"acceptor":True, "donor":True, "maxOO":self.maxOOHBond, "maxAngle":self.maxAngleHBond}
		nTotalBinner = atomComboBinvalGetterHelp._DiscHBondCounterBetweenGroupsWithOxyDistFilterOneDimValGetter(*currArgs, **currKwargs)

		outBinner = atomComboCoreHelp._GetMultiDimValsToBinFromSparseMatrices([minDistBinner, nDonorBinner, nAcceptorBinner, nTotalBinner])

		return outBinner



class _WaterClassifierMinDistHBondsAndAdsSiteHozDists(_WaterClassifierBase):

	def __init__(self, oxyIndices, hyIndices, distFilterIndices, distFilterRange, nDonorFilterRange,
	             nAcceptorFilterRange, nTotalFilterRange, maxOOHBond, maxAngleHBond, adsSiteMinHozToOtherAdsSiteRange, execCount=0):
		""" Initializer
		
		Args:
			oxyIndices: (iter of ints) The oxygen indices for each water molecule
			hyIndices: (iter of len-2 ints) Same length as oxyIndices, but each contains the indices of two hydrogen indices bonded to the relevant oxygen
			distFilterIndices: (iter of ints) Each represents an atom index. We group water by distance of oxygen atoms from these indices
			distFilterRange: (len-2 float iter) [minDist,maxDist] from indices in distFilterIndices for a water to be included in this count
			nDonorFilterRange: (len-2 float iter) [minNDonor, maxNDonor] for a water.
			nAcceptorFilterRange: (len-2 float iter) [minNTotal,maxNTotal] for a water
			nTotalFilterRange: (len-2 float iter) [minNTotal,maxNTotal] for a water
			maxOOHBond: The maximum O-O distance between two hydrogen-bonded water.
			maxAngleHBond:  The maximum OA-OD-HD angle for a hydrogen bond; OA = acceptor oxygen, OD=Donor oxygen, HD=donor hydrogen
			adsSiteMinHozToOtherAdsSiteDistRanges: ( len-2 float iter) [minDist,maxDist] for the closest ads-site to ads-site contact using horizontal distances 
			execCount: (int) Used to track how many times .classify is called; used as a safety check when using "byReference" classifiers
	 
		"""
		self.oxyIndices = oxyIndices
		self.hyIndices = hyIndices
		self.distFilterIndices = distFilterIndices
		self.distFilterRange = distFilterRange
		self.nDonorFilterRange = nDonorFilterRange
		self.nAcceptorFilterRange = nAcceptorFilterRange
		self.nTotalFilterRange = nTotalFilterRange
		self.maxOOHBond = maxOOHBond
		self.maxAngleHBond = maxAngleHBond
		self.adsSiteMinHozToOtherAdsSiteRange = adsSiteMinHozToOtherAdsSiteRange
		self.execCount = execCount

	def classify(self, sparseMatrixCalculator):
		#1)Use another classifier to filter indices based on min-dists to adsorbates/number of h-bonds etc
		currArgs = [self.oxyIndices, self.hyIndices, self.distFilterIndices, self.distFilterRange, self.nDonorFilterRange,
		            self.nAcceptorFilterRange, self.nTotalFilterRange, self.maxOOHBond, self.maxAngleHBond]
		otherClassifier = _WaterClassifierMinDistAndNumberHBonds(*currArgs)
		oxyIndices, hyIndices = otherClassifier.classify(sparseMatrixCalculator)


		#2) Figure out the indices are adsorbed to sites with min(hozDist) in the relevant range
		adsSiteIndices = list()
		for oxyIdx in oxyIndices:
			relRow = sparseMatrixCalculator.outDict["distMatrix"][oxyIdx]
			minIdx = np.argmin(relRow[self.distFilterIndices]) #Shouldnt need the nan version here; 
			adsSiteIndices.append( self.distFilterIndices[minIdx] ) 

		#3) Get the min hoz distance for each
		minHozDistBinValGetter = atomComboBinvalGetterHelp._MinHozDistsGetValsToBin(adsSiteIndices, adsSiteIndices, minVal=0.01)
		outHozDists = minHozDistBinValGetter.getValsToBin(sparseMatrixCalculator)

		#TODO: If two are at the same site we need to count as zero min hoz-dist
		#TODO: Make sure one of the unit tests involves water filtered at THIS point
		#4) Return the relevant oxyIndices and hyIndices
		indicesInList = list()
		for idx,outHozDist in enumerate(outHozDists):
			if  (self.adsSiteMinHozToOtherAdsSiteRange[0]<=outHozDist) and (outHozDist<self.adsSiteMinHozToOtherAdsSiteRange[1]):
				indicesInList.append(idx)

		outOxyIndices = [oxyIndices[idx] for idx in indicesInList]
		outHyIndices  = [hyIndices[idx]  for idx in indicesInList]


		self.execCount += 1
		self.storedClassifyResult = (outOxyIndices, outHyIndices)

		return outOxyIndices, outHyIndices









