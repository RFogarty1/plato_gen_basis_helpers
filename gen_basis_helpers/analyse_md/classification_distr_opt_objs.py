
import itertools as it

from . import binned_res as binResHelp
from . import calc_distrib_core as calcDistrCoreHelp
from . import distr_opt_objs as distrOptObjHelp


class AtomClassifyBasedOnDistsFromIndicesSimpleOpts():
	""" Simple classification options for individual atoms based on their minimum distance from another group of atoms """

	def __init__(self, binResObjs, atomIndices, distFilterIndices, distFilterRanges, minDistVal=-0.01):
		""" Initializer
		
		Args:
			binResObjs: (iter of BinnedResultsStandard objects) One bin for each type of species you want to count (determined by the "Ranges" parameters)
			atomIndices: (iter of ints) The atom indices we want distributions FROM
			distFilterIndices: (iter of ints) The atom indices we want to calculate distances to in order to filter out some "atomIndices"
			distFilterRanges: (iter of len-2 float iters) The ranges used to classify the water types
			minDistVal: (float) If set to a +ve number we ignore distances smaller than it when figuring out minimum. Useful to avoid getting zeros when atomIndices and distFilterIndices overlap
 
		"""
		self.binResObjs = binResObjs
		self.atomIndices = atomIndices
		self.distFilterIndices = distFilterIndices
		self.distFilterRanges = distFilterRanges
		self.minDistVal = minDistVal

#Note: binResObj....different for each? I guess thats the least stupid way, though will make it a not-real options object
class WaterCountTypesMinDistAndHBondSimpleOpts():

	def __init__(self, binResObjs, oxyIndices, hyIndices, distFilterIndices, distFilterRanges , nDonorFilterRanges=None, nAcceptorFilterRanges=None, nTotalFilterRanges=None, maxOOHBond=3.5, maxAngleHBond=35, checkInputConsistent=True):
		""" Initializer
		
		Args:
			binResObjs: (iter of BinnedResultsStandard objects) One bin for each type of water you want to count (determined by the "Ranges" parameters)
			oxyIndices: (iter of ints) The oxygen indices for each water molecule
			hyIndices: (iter of len-2 ints) Same length as oxyIndices, but each contains the indices of two hydrogen indices bonded to the relevant oxygen
			distFilterIndices: (iter of ints) Each represents an atom index. We group water by distance of oxygen atoms from these indices
			distFilterRanges: (iter of len-2 float iters) Each contains [minDist,maxDist] from indices in distFilterIndices for a water to be included in the relevant bin
			nDonorFilterRanges: (iter of len-2 float iters) Each contains [minNDonor, maxNDonor] for a water.
			nAcceptorFilterRanges: (iter of len-2 float iters) Each contains [minNAcceptor, maxNAcceptor] for a water. Setting them to floats just below/above thresholds is sensible (e.g. [-0.1,2.1] for between 0 and two acceptors)
			nTotalFilterRanges: (iter of len-2 float iters) Each contains [minNTotal,maxNTotal] for a water
			maxOOHBond: (float) The maximum O-O distance between two hydrogen-bonded water. Angles are only calculated when this criterion is fulfilled
			maxAngleHBond: (float) The maximum OA-OD-HD angle for a hydrogen bond; OA = acceptor oxygen, OD=Donor oxygen, HD=donor hydrogen
			checkInputConsistent: (Bool) If True, run some checks at initialization to check n-dimensions consistent between the "Ranges" attributes

		Notes:
			a) For [minX,maxX] ranges we define as minX<=x<maxX
			b) nDonor/nAcceptor/nTotal filter ranges default to [-1,100] which should ALWAYS include every water molecule.

		Raises:
			ValueError:

		"""
		self._eqTol = 1e-5

		self.binResObjs = binResObjs
		self.oxyIndices = oxyIndices
		self.hyIndices = hyIndices
		self.distFilterIndices = distFilterIndices
		self.distFilterRanges = distFilterRanges
		self.nDonorFilterRanges = nDonorFilterRanges
		self.nAcceptorFilterRanges = nAcceptorFilterRanges
		self.nTotalFilterRanges = nTotalFilterRanges
		self.maxOOHBond = maxOOHBond
		self.maxAngleHBond = maxAngleHBond

		self._populateDefaultRanges()

		if checkInputConsistent:
			self._checkInputConsistent()


	@classmethod
	def fromWaterIndicesAndGeom(cls, binResObjs, waterIndices, inpGeom, distFilterIndices, distFilterRanges, nDonorFilterRanges=None, nAcceptorFilterRanges=None, nTotalFilterRanges=None, maxOOHBond=3.5, maxAngleHBond=35, checkInputConsistent=True):
		""" See initializer doc-string for args other than waterIndices and inpGeom """
		oxyIndices, hyIndices = distrOptObjHelp._getOxyAndHyIndicesFromWaterIndicesAndGeom(waterIndices, inpGeom)
		currArgs = [binResObjs, oxyIndices, hyIndices, distFilterIndices, distFilterRanges]
		currKwargs = {"nDonorFilterRanges":nDonorFilterRanges, "nAcceptorFilterRanges":nAcceptorFilterRanges, "nTotalFilterRanges":nTotalFilterRanges,
		              "maxOOHBond":maxOOHBond, "maxAngleHBond":maxAngleHBond, "checkInputConsistent":checkInputConsistent}
		return cls(*currArgs, **currKwargs)

	@classmethod
	def fromFilterObjs(cls, binResObjs, oxyIndices, hyIndices, distFilterIndices, filterObjs, maxOOHBond=3.5, maxAngleHBond=35, checkInputConsistent=True):
		""" See initializer doc-string for args other than filterObjs """
		
		#Get vals
		distFilterRanges = [x.distFilterRange for x in filterObjs]
		nDonorFilterRanges = [x.nDonorFilterRange for x in filterObjs]
		nAcceptorFilterRanges = [x.nAcceptorFilterRange for x in filterObjs]
		nTotalFilterRanges = [x.nTotalFilterRange for x in filterObjs]

		#Check for any None values; raise if in distFilterRanges but just set to wide parameters if in nDonor/nAcceptor/nTotal
		distFilterRanges      = None if all([x is None for x in distFilterRanges]) else distFilterRanges
		nDonorFilterRanges    = None if all([x is None for x in nDonorFilterRanges]) else nDonorFilterRanges
		nAcceptorFilterRanges = None if all([x is None for x in nAcceptorFilterRanges]) else nAcceptorFilterRanges
		nTotalFilterRanges    = None if all([x is None for x in nTotalFilterRanges]) else nTotalFilterRanges

		#
		currLists = [distFilterRanges, nDonorFilterRanges, nAcceptorFilterRanges, nTotalFilterRanges]
		for listIdx,currList in enumerate(currLists):
			if currList is not None:
				if any([x is None for x in currList]):
					newList = list()
					for val in currList:
						newVal = [-1,1000] if val is None else val
						newList.append(newVal)
					currLists[listIdx] = newList

		distFilterRanges, nDonorFilterRanges, nAcceptorFilterRanges, nTotalFilterRanges = currLists

		#
		currArgs = [binResObjs, oxyIndices, hyIndices, distFilterIndices, distFilterRanges]
		currKwargs = {"nDonorFilterRanges":nDonorFilterRanges, "nAcceptorFilterRanges":nAcceptorFilterRanges, "nTotalFilterRanges":nTotalFilterRanges,
		              "maxOOHBond":maxOOHBond, "maxAngleHBond":maxAngleHBond, "checkInputConsistent":checkInputConsistent}

		return cls(*currArgs, **currKwargs)



	def _populateDefaultRanges(self):
		attrs = ["nDonorFilterRanges", "nAcceptorFilterRanges", "nTotalFilterRanges"]
		nDims = len(self.distFilterRanges)
		for attr in attrs:
			if getattr(self, attr) is None:
				setattr(self, attr, [ [-1,1000] for x in range(nDims)] )


	def _checkInputConsistent(self):
		relLenAttrs = ["distFilterRanges", "nDonorFilterRanges", "nAcceptorFilterRanges", "nTotalFilterRanges"]
		relLengths = [len(getattr(self,attr)) for attr in relLenAttrs]
		if any([x>relLengths[0] for x in relLengths]):
			raise ValueError("Inconsistent lengths in relLengths")


	def __eq__(self, other):
		eqTol = min(self._eqTol, other._eqTol)

		directCmpAttrs = ["binResObjs", "oxyIndices", "hyIndices", "distFilterIndices"]
		floatCmpAttrs  = ["maxOOHBond", "maxAngleHBond"]
		rangesCmpAttrs = ["distFilterRanges", "nDonorFilterRanges", "nAcceptorFilterRanges", "nTotalFilterRanges"]

		for attr in directCmpAttrs:
			valA, valB = getattr(self,attr), getattr(other,attr)
			if valA != valB:
				return False

		for attr in floatCmpAttrs:
			valA, valB = getattr(self,attr), getattr(other,attr)
			diff = abs(valA-valB)
			if diff > eqTol:
				return False

		for attr in rangesCmpAttrs:
			valsA, valsB = getattr(self,attr), getattr(other,attr)
			if len(valsA) != len(valsB):
				return False
			for valA, valB in it.zip_longest(valsA,valsB):
				diffs = [abs(b-a) for b,a in it.zip_longest(valA,valB)]
				if any([x>eqTol for x in diffs]):
					return False

		return True



class WaterMinDistAndHBondsFilterObj():
	""" Represents some simple options to define a group of water based on hydrogen-bonds counted and min-distance from a group. Used to initialize WaterCountTypesMinDistAndHBondSimpleOpts in a different way """

	def __init__(self, distFilterRange=None, nDonorFilterRange=None, nAcceptorFilterRange=None, nTotalFilterRange=None):
		""" Initializer
		
		Args:
			distFilterRange: (minDist,maxDist) Note that this refers only to the OXYGEN atom distance to the group
			nDonorFilterRange: (minNumber, maxNumber)
			nAcceptorFilterRange: (minNumber, maxNumber)
			nTotalFilterRange: (minNumber, maxNumber)
				 
		"""
		self.distFilterRange = distFilterRange
		self.nDonorFilterRange = nDonorFilterRange
		self.nAcceptorFilterRange = nAcceptorFilterRange
		self.nTotalFilterRange = nTotalFilterRange






