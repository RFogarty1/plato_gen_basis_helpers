
import numpy as np

from . import atom_combo_core as atomComboCoreHelp
from . import atom_combo_binval_getters as atomComboBinvalGetterHelp

from . import classifier_objs as classifierObjHelp
from . import water_derivative_species_classifiers as waterDerivClassifierHelp

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
		classifier = classifierObjHelp._AtomsWithinMinDistRangeClassifier(*currArgs, minDistVal=self.minDistVal)
		relIndices = classifier.classify(sparseMatrixCalculator)

		#2) Just look how many indices are in the list
		return [len(relIndices)]



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
		classifier = classifierObjHelp._WaterClassifierMinDistAndNumberHBonds(*currArgs)
		outOxyIndices, outHyIndices = classifier.classify(sparseMatrixCalculator)

		return [len(outOxyIndices)]


class _WaterClassifierBase():
	""" Classifies list of water indices into distinct groups of water """
	

	def classify(self, sparseMatrixCalculator):
		raise NotImplementedError("")



#Inheriting from classifier since it needs all the same info + the classify step is where the actual work is
class _AdsorbedWaterCountTypeWithAdsSiteHozDistsBinvalGetter(atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase, classifierObjHelp._WaterClassifierMinDistHBondsAndAdsSiteHozDists):

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
		self.execCount = 0 #Annoying hangover from inheriting from classifier

	def getValsToBin(self, sparseMatrixCalculator):
		oxyIndices,hyIndices = self.classify(sparseMatrixCalculator)
		return [len(oxyIndices)]



class _CountTypesBasedOnNumberHBondsToDynamicGroup(classifierObjHelp._ClassiferUsingHBondsToDynamicGroup, atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase):

	def getValsToBin(self, sparseMatrixCalculator):
		nonHyIndices, hyIndices = self.classify(sparseMatrixCalculator)
		return [len(nonHyIndices)]

#Inherit interface from "_GetOneDimValsToBinFromSparseMatricesBase" and behaviour from "_GenericNonHyAndHyClassiferUsingHBondsToGroup_simple"
class _GenericCountTypesBasedOnNumberHBondsToGroup(classifierObjHelp._GenericNonHyAndHyClassiferUsingHBondsToGroup_simple, atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase):
	pass


	def getValsToBin(self, sparseMatrixCalculator):
		nonHyIndices,hyIndices = self.classify(sparseMatrixCalculator)

		return [len(nonHyIndices)]


class _CountWaterDerivativeDistanceOnlyBinvalGetter(waterDerivClassifierHelp._WaterDerivativeDistanceOnlyClassifierGeneric, atomComboCoreHelp._GetOneDimValsToBinFromSparseMatricesBase):
	pass

	def getValsToBin(self, sparseMatrixCalculator):
		nonHyIndices, hyIndices = self.classify(sparseMatrixCalculator)
		return [len(nonHyIndices)]





