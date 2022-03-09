
import itertools as it

from . import binned_res as binResHelp
from . import calc_distrib_core as calcDistrCoreHelp
from . import distr_opt_objs as distrOptObjHelp


class CompositeClassiferOptsSimple():
	""" Class for when you need to mix multiple types of classification opts; such as assigning one group based on planar distance and another on number of h-bonds """

	def __init__(self, optsObjs):
		""" Initializer
		
		Args:
			optsObjs: (iter of options objects) e.g. original purpose was to combine [ClassifyBasedOnHBondingToGroup_simple, ClassifyBasedOnHBondingToDynamicGroup]
				
		Notes:
			The optsObjs should probably lead to classifiers with the same interface; but i wont check that in any way.
 
		"""
		self.optsObjs = optsObjs




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

		#Sort input stuff
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

		#Classify attributes so we can more easily make comparisons/populate various defaults

		#Hidden parts for equality checks
		self._eqTol = 1e-5
		self._eqDirectCmpAttrs = ["binResObjs", "oxyIndices", "hyIndices", "distFilterIndices"]
		self._eqFloatCmpAttrs = ["maxOOHBond", "maxAngleHBond"]
		self._eqRangesCmpAttrs = ["distFilterRanges", "nDonorFilterRanges", "nAcceptorFilterRanges", "nTotalFilterRanges"]

		#Hidden parts for setting defaults/figuring out lengths
		#Note: Convert these into lists in the function as needed; tuples are safer here (to two dist-ranges dont share the same mem-address)
		self._defRangeVals = {"nDonorFilterRanges":(-1,1000), "nAcceptorFilterRanges":(-1,1000), "nTotalFilterRanges":(-1,1000)}
		self._rangeAttrs = ["distFilterRanges", "nDonorFilterRanges", "nAcceptorFilterRanges", "nTotalFilterRanges"] #Need to check all same length so....


		#Obv the end points really
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
		nDims = len(self.distFilterRanges)
		for attr in self._defRangeVals.keys():
			if getattr(self, attr) is None:
				currVal = self._defRangeVals[attr]
				setattr(self, attr, [ list(currVal) for x in range(nDims)] )

	def _checkInputConsistent(self):
		relLengths = [len(getattr(self,attr)) for attr in self._rangeAttrs]
		if any([x>relLengths[0] for x in relLengths]):
			raise ValueError("Inconsistent lengths in relLengths")

	def __eq__(self, other):
		eqTol = min(self._eqTol, other._eqTol)

		for attr in self._eqDirectCmpAttrs:
			valA, valB = getattr(self,attr), getattr(other,attr)
			if valA != valB:
				return False

		for attr in self._eqFloatCmpAttrs:
			valA, valB = getattr(self,attr), getattr(other,attr)
			diff = abs(valA-valB)
			if diff > eqTol:
				return False

		for attr in self._eqRangesCmpAttrs:
			valsA, valsB = getattr(self,attr), getattr(other,attr)
			if len(valsA) != len(valsB):
				return False
			for valA, valB in it.zip_longest(valsA,valsB):
				diffs = [abs(b-a) for b,a in it.zip_longest(valA,valB)]
				if any([x>eqTol for x in diffs]):
					return False

		return True




class WaterAdsorbedClassifier_usingMinHozDistsBetweenAdsorptionSitesOptsObj(WaterCountTypesMinDistAndHBondSimpleOpts):
	""" This is an options object specifically for finding certain types of adsorbed water; based on the horizontal distances betwee the atoms which they are adsorbed on. Its actually mostly inherited from WaterCountTypesMinDistAndHBondSimpleOpts though """


	def __init__(self, binResObjs, oxyIndices, hyIndices, distFilterIndices, distFilterRanges , nDonorFilterRanges=None, nAcceptorFilterRanges=None, nTotalFilterRanges=None, maxOOHBond=3.5, maxAngleHBond=35, checkInputConsistent=True, adsSiteMinHozToOtherAdsSiteRanges=None):
		""" Initializer
		
		Args:
			binResObjs: (iter of BinnedResultsStandard objects) One bin for each type of water you want to count (determined by the "Ranges" parameters)
			oxyIndices: (iter of ints) The oxygen indices for each water molecule
			hyIndices: (iter of len-2 ints) Same length as oxyIndices, but each contains the indices of two hydrogen indices bonded to the relevant oxygen
			distFilterIndices: (iter of ints) Each represents an atom index. We group water by distance of oxygen atoms from these indices (so they should be adsorption sites really)
			distFilterRanges: (iter of len-2 float iters) Each contains [minDist,maxDist] from indices in distFilterIndices for a water to be included in the relevant bin
			nDonorFilterRanges: (iter of len-2 float iters) Each contains [minNDonor, maxNDonor] for a water.
			nAcceptorFilterRanges: (iter of len-2 float iters) Each contains [minNAcceptor, maxNAcceptor] for a water. Setting them to floats just below/above thresholds is sensible (e.g. [-0.1,2.1] for between 0 and two acceptors)
			nTotalFilterRanges: (iter of len-2 float iters) Each contains [minNTotal,maxNTotal] for a water
			maxOOHBond: (float) The maximum O-O distance between two hydrogen-bonded water. Angles are only calculated when this criterion is fulfilled
			maxAngleHBond: (float) The maximum OA-OD-HD angle for a hydrogen bond; OA = acceptor oxygen, OD=Donor oxygen, HD=donor hydrogen
			checkInputConsistent: (Bool) If True, run some checks at initialization to check n-dimensions consistent between the "Ranges" attributes
			adsSiteMinHozToOtherAdsSiteDistRanges: (iter of len-2 float iters) [minDist,maxDist] for the closest ads-site to ads-site contact using horizontal distances 

		Notes:
			a) For [minX,maxX] ranges we define as minX<=x<maxX
			b) nDonor/nAcceptor/nTotal filter ranges default to [-1,100] which should ALWAYS include every water molecule.

		Raises:
			ValueError:

		"""
		#0) Initialize using the parent class
		currArgs = [binResObjs, oxyIndices, hyIndices, distFilterIndices, distFilterRanges]
		currKwargs = {"nDonorFilterRanges": nDonorFilterRanges, "nAcceptorFilterRanges": nAcceptorFilterRanges, "nTotalFilterRanges": nTotalFilterRanges,
		              "maxOOHBond": maxOOHBond, "maxAngleHBond": maxAngleHBond, "checkInputConsistent": checkInputConsistent}

		super().__init__(*currArgs, **currKwargs)

		#1) Set remaining (extended) attributes
		self.adsSiteMinHozToOtherAdsSiteRanges = adsSiteMinHozToOtherAdsSiteRanges

		#2) Extend the various hidden attributes to deal with the new ones 
		self._eqRangesCmpAttrs.append("adsSiteMinHozToOtherAdsSiteRanges") 

		self._defRangeVals["adsSiteMinHozToOtherAdsSiteRanges"] = (0.01,1000)
		self._rangeAttrs.append( "adsSiteMinHozToOtherAdsSiteRanges" )

		#Various cleanup:
		self._populateDefaultRanges()
		if checkInputConsistent:
			self._checkInputConsistent()

	@classmethod
	def fromFilterObjs(cls, binResObjs, oxyIndices, hyIndices, distFilterIndices, filterObjs, maxOOHBond=3.5, maxAngleHBond=35, checkInputConsistent=True):
		""" See initializer doc-string for args other than filterObjs """

		#Get the range attributes
		rangeAttrDefaults = {"distFilterRange": (0,1000), "nDonorFilterRange":(-1,1000), "nAcceptorFilterRange":(-1,1000),
		                      "nTotalFilterRange":(-1,1000), "minHozDistAdsSitesRange":(0.01,1000)}
		rangeAttrDict = dict()

		for key in rangeAttrDefaults.keys():
			currVals = list()
			for filterObj in filterObjs:
				currVal = getattr(filterObj, key)
				currVal = currVal if currVal is not None else rangeAttrDefaults[key]
				currVals.append(currVal)
			rangeAttrDict[key] = currVals

		#Get args/kwargs and pass them to initializer
		currArgs = [binResObjs, oxyIndices, hyIndices, distFilterIndices, rangeAttrDict["distFilterRange"]]
		currKwargs = {"nDonorFilterRanges": rangeAttrDict["nDonorFilterRange"], "nAcceptorFilterRanges": rangeAttrDict["nAcceptorFilterRange"],
		              "nTotalFilterRanges": rangeAttrDict["nTotalFilterRange"], "adsSiteMinHozToOtherAdsSiteRanges": rangeAttrDict["minHozDistAdsSitesRange"],
		              "maxOOHBond": maxOOHBond, "maxAngleHBond": maxAngleHBond, "checkInputConsistent": checkInputConsistent}

		return cls(*currArgs, **currKwargs)


class ClassifyNonHyAndHyChainedAllCommon():
	""" Classifies based on criterion being fulfilled from multiple input options objects. For example, an object for minmum number of h-bonds to one group AND min horizontal-distances to another.

	"""

	def __init__(self, classifierOpts):
		""" Initializer
		
		Args:
			classifierOpts: (iter of classifier options) These need to  lead to classifiers that return ([NonHyIndices],[HyIndices]) from their classify method
				 
		"""
		self.classifierOpts = classifierOpts

class ClassifyBasedOnHBondingToDynamicGroup():
	""" Class designed to classift molecules (e.g. water/hydroxyl) based on number of hydrogen bonds they have to a dynamically assigned group.

		Original purpose was to filter to water molecules (solCon) h-bonded to a dynamically-assinged group of water which had h-bonds to hydroxyl groups (solAds). This was for the hydroxylated Mg/water interface
	"""

	def __init__(self, dynGroupOptObj, thisGroupOptObj, mutuallyExclusive=True, checkConsistent=True, firstClassifierObjs=None):
		""" Initializer
		
		Args:
			dynGroupOptObj: (ClassifyBasedOnHBondingToGroup_simple options obj) This contains options for dynamically assigning the first group
			thisGroupOptObj: (ClassifyBasedOnHBondingToGroup_simple options obj) This contains options for assigning this group based on h-bonding with the "dynamic" group
			mututallyExclusive: (Bool) If True then a molecule cant be in both groups; it will be assigned to ONLY dynGroupOptObj
			checkConsistent: (Bool) Carries out various tests to catch possible inconsistencies between the two input opts objects (for example, .toIndices on thisGroupOptObj need to match .fromIndices on dynGroupOptObj)
			firstClassifierObj: (iter of ClassifierBase object) This is the classifier used to get the first group (dynGroupOptObj). Passing here is just a possible way of speeding things up (assuming you pass a byReference classifier; will avoid redoing painful classification work). NOTE: Need one per nTotalFilterRange element on your input groups

		Notes:
			1) Length of filter ranges on dynGroupOptObj and thisGroupOptObj need to be the same (checked with checkConsisten=True)
			2) The thisGroupObjObj toIndices need to be the same as dynGroupOptObj fromIndices (checked with checkConsistent=True)
			3) Dont RELY on checkConsistent too much; I havent really tested every way you can get input wrong.
 
		"""
		self.dynGroupOptObj = dynGroupOptObj
		self.thisGroupOptObj = thisGroupOptObj
		self.mutuallyExclusive = mutuallyExclusive
		self.firstClassifierObjs = firstClassifierObjs

		if checkConsistent:
			self._checkInputConsistent()

	def _checkInputConsistent(self):
		#Check the to/from Indices are conssitent
		if self.dynGroupOptObj.fromNonHyIndices != self.thisGroupOptObj.toNonHyIndices:
			valsA, valsB = self.dynGroupOptObj.fromNonHyIndices, self.thisGroupOptObj.toNonHyIndices
			raise ValueError("fromNonHyIndices != toNonHyIndices; values are \n {} and {}".format(valsA, valsB))

		if self.dynGroupOptObj.fromHyIndices != self.thisGroupOptObj.toHyIndices:
			valsA, valsB = self.dynGroupOptObj.fromHyIndices, self.thisGroupOptObj.toHyIndices
			raise ValueError("fromHyIndices != toHyIndices; values are \n {} and {}".format(valsA,valsB))

		#Check the filter ranges are consistent
		filterAttrs = ["nDonorFilterRanges", "nAcceptorFilterRanges", "nTotalFilterRanges"]
		for filterAttr in filterAttrs:
			lenA, lenB = len( getattr(self.dynGroupOptObj, filterAttr) ), len( getattr(self.thisGroupOptObj, filterAttr) )
			if lenA != lenB:
				currArgs = [filterAttr, lenA, lenB]
				raise ValueError("Different length filter ranges detected for {}; lengths are {} and {}".format(*currArgs))



class ClassifyNonHyAndHyBasedOnMinHozDistsToAtomGroups():
	""" Class designed to filter NonHy/Hy indices out based on their horizontal distances with other atoms
	
	"""

	def __init__(self, binResObjs, fromNonHyIndices, fromHyIndices, toNonHyIndices, toHyIndices, distFilterRanges, useIndicesFrom="nonHy", useIndicesTo="nonHy" ,minDistVal=-0.01):
		""" Initializer
		
		Args:
			binResObjs: (iter of BinnedResultsStandard objects) One bin for each type of water you want to count (determined by the "Ranges" parameters)
			fromNonHyIndices: (iter of iter of ints) The non-hydrogen indices of each molecule we're filtering
			fromHyIndices: (iter of iter of ints) Same length as nonHyFromIndices, but contain the relevant hydrogen indices
			toNonHyIndices: (iter of iter of ints) The non-hydrogen indices for all molecules we're counting hydrogen bonds TO (e.g could be hydroxyl molecules if we're filtering for water molecules with h-bonds TO hydroxyls)
			toHyIndices: (iter of iter of ints) The hydrogen indices for all molecules we're counting hydrogen bonds TO
			distFilterRanges: (iter of len-2 float iters) Each contains [minDist, maxDist] for a molecule
			useIndicesFrom: (str) - "nonHy", "hy" or "all". Which atoms we need to include for "fromNonHyIndices" and "fromHyIndices"
			useIndicesTo: (str) - "nonHy", "hy" or "all". See above
			minDistVal: (float) If set to a +ve number we ignore distances smaller than it when figuring out minimum. Useful to avoid getting zeros when atomIndices and distFilterIndices overlap

		"""
		self.binResObjs = binResObjs
		self.fromNonHyIndices = fromNonHyIndices
		self.fromHyIndices = fromHyIndices
		self.toNonHyIndices = toNonHyIndices
		self.toHyIndices = toHyIndices
		self.distFilterRanges = distFilterRanges
		self.useIndicesFrom = useIndicesFrom
		self.useIndicesTo = useIndicesTo
		self.minDistVal = minDistVal


class ClassifyBasedOnHBondingToGroup_simple():
	""" Class designed to classify molecules (e.g. water/hydroxyl) based on the number of hydrogen bonds they have to a static group. 

	Original purpose was to classify water based on number of h-bonds to hydroxyl and vice-versa
	"""

	def __init__(self, binResObjs, fromNonHyIndices, fromHyIndices, toNonHyIndices, toHyIndices, nDonorFilterRanges=None, nAcceptorFilterRanges=None, nTotalFilterRanges=None,
	             maxOOHBond=3.5, maxAngleHBond=35, checkInputConsistent=True):
		""" Initializer
		
		Args:
			binResObjs: (iter of BinnedResultsStandard objects) One bin for each type of water you want to count (determined by the "Ranges" parameters)
			fromNonHyIndices: (iter of iter of ints) The non-hydrogen indices of each molecule we're filtering
			fromHyIndices: (iter of iter of ints) Same length as nonHyFromIndices, but contain the relevant hydrogen indices
			toNonHyIndices: (iter of iter of ints) The non-hydrogen indices for all molecules we're counting hydrogen bonds TO (e.g could be hydroxyl molecules if we're filtering for water molecules with h-bonds TO hydroxyls)
			toHyIndices: (iter of iter of ints) The hydrogen indices for all molecules we're counting hydrogen bonds TO
			nDonorFilterRanges: (iter of len-2 float iters) Each contains [minNDonor, maxNDonor] for a molecule.
			nAcceptorFilterRanges: (iter of len-2 float iters) Each contains [minNAcceptor, maxNAcceptor] for a molecule. Setting them to floats just below/above thresholds is sensible (e.g. [-0.1,2.1] for between 0 and two acceptors)
			nTotalFilterRanges: (iter of len-2 float iters) Each contains [minNTotal,maxNTotal] for a molecule
			maxOOHBond: (float) The maximum O-O distance between two hydrogen-bonded water. Angles are only calculated when this criterion is fulfilled
			maxAngleHBond: (float) The maximum OA-OD-HD angle for a hydrogen bond; OA = acceptor oxygen, OD=Donor oxygen, HD=donor hydrogen [NOTE: Just replace OA and OD with relevant other atoms for non-water)
			checkInputConsistent: (Bool) If True, run some checks at initialization to check n-dimensions consistent between the "Ranges" attributes

		Notes:
			a) For [minX,maxX] ranges we define as minX<=x<maxX
			b) nDonor/nAcceptor/nTotal filter ranges default to [-1,100] which should ALWAYS include every molecule.

		Raises:
			ValueError:

		"""

		self.binResObjs = binResObjs
		self.fromNonHyIndices = fromNonHyIndices
		self.fromHyIndices = fromHyIndices
		self.toNonHyIndices = toNonHyIndices
		self.toHyIndices = toHyIndices
		self.nDonorFilterRanges = nDonorFilterRanges
		self.nAcceptorFilterRanges = nAcceptorFilterRanges
		self.nTotalFilterRanges = nTotalFilterRanges
		self.maxOOHBond = maxOOHBond
		self.maxAngleHBond = maxAngleHBond


		#Hidden parts for setting defaults/figuring out lengths
		#Note: Convert these into lists in the function as needed; tuples are safer here (to two dist-ranges dont share the same mem-address)
		self._defRangeVals = {"nDonorFilterRanges":(-1,1000), "nAcceptorFilterRanges":(-1,1000), "nTotalFilterRanges":(-1,1000)}
		self._rangeAttrs = ["nDonorFilterRanges", "nAcceptorFilterRanges", "nTotalFilterRanges"] #Need to check all same length so....

		#Various cleanup:
		self._populateDefaultRanges()
		if checkInputConsistent:
			self._checkInputConsistent()


	def _populateDefaultRanges(self):
		#Figure out number of dimensions
		rangeVals = [getattr(self,attr) for attr in self._rangeAttrs]
		if all([x is None for x in rangeVals]):
			raise ValueError("Must set at least one of the filter ranges")
		else:
			rangeVals = [x for x in rangeVals if x is not None]
			assert all([len(x)==len(rangeVals[0]) for x in rangeVals])
			nDims = len(rangeVals[0])

		#Set values
		for attr in self._defRangeVals.keys():
			if getattr(self, attr) is None:
				currVal = self._defRangeVals[attr]
				setattr(self, attr, [ list(currVal) for x in range(nDims)] )

	def _checkInputConsistent(self):
		relLengths = [len(getattr(self,attr)) for attr in self._rangeAttrs]
		if any([x>relLengths[0] for x in relLengths]):
			raise ValueError("Inconsistent lengths in relLengths")



class WaterDerivativeBasedOnDistanceClassifierOptsObj():
	""" Classify water-derivatives based on how many hydrogen neighbours each oxygen has. A hydrogen is considered a neighbour if its within a cutoff and closest to THAT oxygen out of all oxyIndices. General purpose is to get indices for water/similar from a trajectory step """

	def __init__(self, binResObjs, oxyIndices, hyIndices, maxOHDist=1.3, nNebs=None):
		""" Initializer
		
		Args:
			binResObjs: (iter of BinnedResultsStandard objects) One bin for each type of water you want to count (determined by the "Ranges" parameters)
			oxyIndices: (iter of ints) Indices of all oxygen atoms to consider
			hyIndices: (iter of ints) Indices of all hydrogen atoms to consider
			maxOHDist: (float) Maximum distance between oxygen/hydrogen for them to be considered "bonded"
			nNebs: (iter of int) Number of hydrogen neighbours needed. 0=free oxygen, 1=hydroxyl, 2=water, 3=hydronium. One entry per type you want to count

		"""
		self.binResObjs = binResObjs
		self.oxyIndices = oxyIndices
		self.hyIndices = hyIndices
		self.maxOHDist = maxOHDist
		self.nNebs = nNebs


class ClassifyByNumberNebsWithinDistanceOptsObj():
	""" Classify hydrogen ATOMS based on whether it has a neighbour within cutoff in the target groups.

	Original use case was to find free hydrogen and H2 within a simulation where water split
	"""

	def __init__(self, binResObjs, fromIndices, toIndices, maxDist=1.3, minDist=0.01, nebRanges=None):
		""" Description of function
		
		Args:
			binResObjs: (iter of BinnedResultsStandard objects) One bin for each type of water you want to count (determined by the "Ranges" parameters)
			fromIndices: (iter of ints) Indices of all hydrogen atoms to consider
			toIndices: (iter of ints) Indices of all OTHER atoms to consider (i.e. all atoms that might be considered neghbours). toIndices=fromIndices is allowed
			maxDist: (float) The maximum distance for an atom to be considered as a neighbour
			minDist: (float) The minimum distance for an atom to be considered as a neighbour; useful to have set low to avoid self-counting
			nebRanges: (iter of len-2 floats) Number of neighbours needed to fit into a classification. e.g. [ [-0.1,0.1], [0.9,1.1] ] will have one group with zero neighbours and one group with a single neighbour 
 
		"""
		self.binResObjs = binResObjs
		self.fromIndices = fromIndices
		self.toIndices = toIndices
		self.maxDist = maxDist
		self.minDist = minDist
		self.nebRanges = nebRanges

	@property
	def atomIndices(self):
		return self.fromIndices


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


class WaterMinDistHbondsAndMinHozDistAdsSiteFilterObj():

	def __init__(self, distFilterRange=None, nDonorFilterRange=None, nAcceptorFilterRange=None, nTotalFilterRange=None, minHozDistAdsSitesRange=None):
		""" Initializer
		
		Args:
			distFilterRange: (minDist,maxDist) Note that this refers only to the OXYGEN atom distance to the group
			nDonorFilterRange: (minNumber, maxNumber)
			nAcceptorFilterRange: (minNumber, maxNumber)
			nTotalFilterRange: (minNumber, maxNumber)
			minHozDistAdsSitesRange: (minDist,maxDist)
 
		"""
		self.distFilterRange = distFilterRange
		self.nDonorFilterRange = nDonorFilterRange
		self.nAcceptorFilterRange = nAcceptorFilterRange
		self.nTotalFilterRange = nTotalFilterRange
		self.minHozDistAdsSitesRange = minHozDistAdsSitesRange


