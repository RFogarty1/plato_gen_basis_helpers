

from . import atom_combo_core as atomComboCoreHelp
from . import atom_combo_binval_getters as atomComboBinvalGetterHelp

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

