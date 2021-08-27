
import numpy as np

from . import atom_combo_core as atomComboCoreHelp
from . import atom_combo_binval_getters as atomComboBinvalGetterHelp



class _AtomsWithMinDistRangeCountBinvalGetter(atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase):

	def __init__(self, atomIndices, distFilterIndices, distFilterRange, minDistVal=-0.01):
		""" Initializer
		
		Args:
			atomIndices: (iter of ints)
			distFilterIndices: (iter of ints) Each represents an atom index. We group atomIndices by min-distance from these indices
			distFilterRange: (len-2 float iter) [minDist,maxDist] from indices in distFilterIndices for an atom to be included in the counts
			minDistVal: (float) If set to a +ve number we ignore distances smaller than it when figuring out minimum. Useful to avoid getting zeros when atomIndices and distFilterIndices overlap
		"""
		self.atomIndices = atomIndices
		self.distFilterIndices = distFilterIndices
		self.distFilterRange = distFilterRange
		self.minDistVal = minDistVal

	def getValsToBin(self, sparseMatrixCalculator):
		#1) Get a classifier to do the main work
		currArgs = [self.atomIndices, self.distFilterIndices, self.distFilterRange]
		classifier = _AtomsWithinMinDistRangeClassifier(*currArgs, minDistVal=self.minDistVal)
		relIndices = classifier.classify(sparseMatrixCalculator)

		#2) Just look how many indices are in the list
		return [len(relIndices)]

class _AtomsWithinMinDistRangeClassifier():

	def __init__(self, atomIndices, distFilterIndices, distFilterRange, minDistVal=-0.01):
		""" Initializer
		
		Args:
			atomIndices: (iter of ints)
			distFilterIndices: (iter of ints) Each represents an atom index. We group atomIndices by min-distance from these indices
			distFilterRange: (len-2 float iter) [minDist,maxDist] from indices in distFilterIndices for an atom to be included in the counts
			minDistVal: (float) If set to a +ve number we ignore distances smaller than it when figuring out minimum. Useful to avoid getting zeros when atomIndices and distFilterIndices overlap
		"""
		self.atomIndices = atomIndices
		self.distFilterIndices = distFilterIndices
		self.distFilterRange = distFilterRange
		self.minDistVal = minDistVal

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

		return outIndices


#Some water options below



class _WaterCountTypeBinvalGetter(atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase):
	
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

	def getValsToBin(self, sparseMatrixCalculator):
		currArgs = [self.oxyIndices, self.hyIndices, self.distFilterIndices, self.distFilterRange, self.nDonorFilterRange, self.nAcceptorFilterRange,
		            self.nTotalFilterRange, self.maxOOHBond, self.maxAngleHBond ]
		classifier = _WaterClassifierMinDistAndNumberHBonds(*currArgs)
		outOxyIndices, outHyIndices = classifier.classify(sparseMatrixCalculator)

		return [len(outOxyIndices)]


class _WaterClassifierBase():
	""" Classifies list of water indices into distinct groups of water """
	

	def classify(self, sparseMatrixCalculator):
		raise NotImplementedError("")


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
	             nAcceptorFilterRange, nTotalFilterRange, maxOOHBond, maxAngleHBond, adsSiteMinHozToOtherAdsSiteRange):
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

		return outOxyIndices, outHyIndices


#Inheriting from classifier since it needs all the same info + the classify step is where the actual work is
class _AdsorbedWaterCountTypeWithAdsSiteHozDistsBinvalGetter(atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase, _WaterClassifierMinDistHBondsAndAdsSiteHozDists):

	#This is the SAME as the classifier..... I could probably just combine the objects or similar (inherit from classifier and add binval getter here
	def __init__(self, oxyIndices, hyIndices, distFilterIndices, distFilterRange, nDonorFilterRange,
	             nAcceptorFilterRange, nTotalFilterRange, maxOOHBond, maxAngleHBond, adsSiteMinHozToOtherAdsSiteRange):
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

	def getValsToBin(self, sparseMatrixCalculator):
		oxyIndices,hyIndices = self.classify(sparseMatrixCalculator)
		return [len(oxyIndices)]

