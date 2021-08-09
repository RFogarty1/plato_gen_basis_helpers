

from . import calc_distrib_core as calcDistrCoreHelp
from . import calc_radial_distrib_impl as calcRadImpl




class CalcRdfOptions(calcDistrCoreHelp.CalcRdfOptions):
	pass


class CalcPlanarDistOptions(calcRadImpl.CalcPlanarRdfOptions):
	pass


class _WaterOptsMixin():

	@property
	def primaryIndices(self):
		if self.primaryIdxType.upper() == "O":
			return self.oxyIndices
		elif self.primaryIdxType.upper() == "HA":
			return [x[0] for x in self.hyIndices]
		elif self.primaryIdxType.upper() == "HB":
			return [x[1] for x in self.hyIndices]
		else:
			raise ValueError("primaryIdxType = {} is an invalid value".format(self.primaryIdxType))

	def _getOxyAndHyIndicesFromWaterIndicesAndGeom(self, waterIndices, inpGeom):
		cartCoords = inpGeom.cartCoords
		eleList = [x[-1] for x in cartCoords]
		outOxyIndices = list()
		outHyIndices = list()
		for waterIdxList in waterIndices:
			currOxy, currHy = list(),list()
			#Sort out the indices
			for idx in waterIdxList:
				currEle = eleList[idx]
				if currEle.upper() == "O":
					currOxy.append( idx )
				elif currEle.upper() == "H":
					currHy.append( idx )
				else:
					raise ValueError("{} is in invalid element for water indices".format(currEle))

			#Add to lists
			assert len(currOxy)==1
			assert len(currHy)==2
			outOxyIndices.extend(currOxy)
			outHyIndices.append(currHy)
		return outOxyIndices, outHyIndices



class WaterMinPlanarDistOptions(calcDistrCoreHelp.CalcDistribOptionsBase, _WaterOptsMixin):

	def __init__(self, binResObj, oxyIndices, hyIndices, planeEqn=None, primaryIdxType="O", minDistType="all"):
		""" Initializer
		
		Args:
			binResObj: (BinnedResultsStandard object) Note that this may get modified in place
			oxyIndices: (iter of ints) The oxygen indices for each water molecule
			hyIndices: (iter of len-2 ints) Same length as oxyIndices, but each contains the indices of two hydrogen indices bonded to the relevant oxygen
			planeEqn: (None or ThreeDimPlaneEquation) The plane equation to calculate distrib function from
			primaryIdxType: (str) The element of the primary index. "O", "Ha" and "Hb" are the standard options
			minDistType: (str) Controls which atoms to get the minimum distance from. Current options are "all","o", and "h" (case insensitive)
				 
		"""
		self.distribKey = "rdf"
		self.binResObj = binResObj
		self.oxyIndices = oxyIndices
		self.hyIndices = hyIndices
		self.planeEqn = planeEqn
		self.primaryIdxType = primaryIdxType
		self.minDistType = minDistType

	@classmethod
	def fromWaterIndicesAndGeom(cls, binResObj, waterIndices, inpGeom, planeEqn=None, primaryIdxType="O", minDistType="all"):
		""" Alternative initializer
		
		Args:
			binResObj: (BinnedResultsStandard object) Note that this may get modified in place
			waterIndices: (iter of len-3 int iters) Each element contains the indices of a water molecule
			inpGeom: (plato_pylib UnitCell object) Used to figure out which indices correspond to Oxygen/Hydrogen
			planeEqn: (None or ThreeDimPlaneEquation) The plane equation to calculate distrib function from
			primaryIdxType: (str) The element of the primary index. "O", "Ha" and "Hb" are the standard options
			minDistType: (str)
 
		"""
		oxyIndices, hyIndices = cls._getOxyAndHyIndicesFromWaterIndicesAndGeom(None,waterIndices, inpGeom)
		currArgs = [binResObj, oxyIndices, hyIndices]
		currKwargs = {"planeEqn":planeEqn, "primaryIdxType":primaryIdxType, "minDistType":minDistType}
		return cls(*currArgs, **currKwargs)


	def __eq__(self, other):
		cmpAttrs = ["binResObj", "oxyIndices", "hyIndices", "planeEqn", "primaryIdxType", "minDistType"]
		for attr in cmpAttrs:
			valA, valB = getattr(self,attr), getattr(other,attr)
			if valA != valB:
				return False

		return True



class WaterPlanarDistOptions(calcDistrCoreHelp.CalcDistribOptionsBase, _WaterOptsMixin):

	def __init__(self, binResObj, oxyIndices, hyIndices, planeEqn=None, primaryIdxType="O"):
		""" Initializer
		
		Args:
			binResObj: (BinnedResultsStandard object) Note that this may get modified in place
			oxyIndices: (iter of ints) The oxygen indices for each water molecule
			hyIndices: (iter of len-2 ints) Same length as oxyIndices, but each contains the indices of two hydrogen indices bonded to the relevant oxygen
			planeEqn: (None or ThreeDimPlaneEquation) The plane equation to calculate distrib function from
			primaryIdxType: (str) The element of the primary index. "O", "Ha" and "Hb" are the standard options

		"""
		self.distribKey = "rdf"
		self.binResObj = binResObj
		self.oxyIndices = oxyIndices
		self.hyIndices = hyIndices
		self.planeEqn = planeEqn
		self.primaryIdxType = primaryIdxType


	@classmethod
	def fromWaterIndicesAndGeom(cls, binResObj, waterIndices, inpGeom, planeEqn=None, primaryIdxType="O"):
		""" Alternative initializer
		
		Args:
			binResObj: (BinnedResultsStandard object) Note that this may get modified in place
			waterIndices: (iter of len-3 int iters) Each element contains the indices of a water molecule
			inpGeom: (plato_pylib UnitCell object) Used to figure out which indices correspond to Oxygen/Hydrogen
			planeEqn: (None or ThreeDimPlaneEquation) The plane equation to calculate distrib function from
			primaryIdxType: (str) The element of the primary index. "O", "Ha" and "Hb" are the standard options
				 
		"""
		oxyIndices, hyIndices = cls._getOxyAndHyIndicesFromWaterIndicesAndGeom(None,waterIndices, inpGeom)
		currArgs = [binResObj, oxyIndices, hyIndices]
		currKwargs = {"planeEqn":planeEqn, "primaryIdxType":primaryIdxType}
		return cls(*currArgs, **currKwargs)


	def __eq__(self, other):
		cmpAttrs = ["binResObj", "oxyIndices", "hyIndices", "planeEqn", "primaryIdxType"]
		for attr in cmpAttrs:
			valA, valB = getattr(self,attr), getattr(other,attr)
			if valA != valB:
				return False

		return True



class DiscHBondCounterWithOxyDistFilterOptions(calcDistrCoreHelp.CalcDistribOptionsBase):
	""" Object contains options for counting discrete number of hydrogen bonds from one group to another """

	def __init__(self, binResObj, oxyIndices, hyIndices, distFilterIndices=None, distFilterVals=None, acceptor=True, donor=True, maxOO=3.5, maxAngle=35):
		""" Initializer
		
		Args:
			binResObj: (BinnedResultsStandard object) Note that this may get modified in place
			oxyIndices: (iter of ints) The oxygen indices for each water molecule
			hyIndices: (iter of len-2 ints) Same length as oxyIndices, but each contains the indices of two hydrogen indices bonded to the relevant oxygen
			distFilterIndices: (iter of ints) These are the indices we look at X-O distances for. They are used to divide water molecules into two groups. Default=None; meaning ALL water are placed in both groups (so we count hydrogen bonds between all water)
			distFilterVals: (iter) min/max distances to be in the two groups. [ [0,3], [3,5] ] would mean groupA are 0<=x<3 from distFilterIndices while groupB are 3<=x<5 from them. Default vals are [ [0,1000], [0,1000] ]
			acceptor: (Bool) If True we calculate dists/angles required to count number of groupA acceptors from groupB. Note changing the order of distFilterValues will give the reverse info (groupB acceptors from groupA)
			donor: (Bool) If True we calculate dists/angles required to count number of groupA donors to groupB. Note changing the order of distFilterValues will give the reverse info (groupB donors to groupA)
			maxOO: (float) The maximum O-O distance between two hydrogen-bonded water. Angles are only calculated when this criterion is fulfilled
			maxAngle: (float) The maximum OA-OD-HD angle for a hydrogen bond; OA = acceptor oxygen, OD=Donor oxygen, HD=donor hydrogen
 
		"""
		self.distribKey = "nHbonds"
		self.binResObj = binResObj
		self.oxyIndices = oxyIndices
		self.hyIndices = hyIndices
		self.distFilterVals = distFilterVals #Default to [[0,1000],[0,1000]
		self.distFilterIndices = distFilterIndices #Default to None; means dont filter anything
		self.donor = donor
		self.acceptor = acceptor
		self.maxOO = maxOO
		self.maxAngle = maxAngle

	@property
	def primaryIndices(self):
		return self.oxyIndices

