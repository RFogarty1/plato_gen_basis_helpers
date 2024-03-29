

from . import binned_res as binResHelp
from . import calc_distrib_core as calcDistrCoreHelp
from . import calc_radial_distrib_impl as calcRadImpl



class CalcRdfOptions(calcDistrCoreHelp.CalcRdfOptions):
	pass


class CalcPlanarDistOptions(calcRadImpl.CalcPlanarRdfOptions):
	pass


class CalcHozDistOptions(calcDistrCoreHelp.CalcDistribOptionsBase):

	def __init__(self, binResObj, indicesFrom, indicesTo, minDistAToB=False, minDistVal=0.01):
		""" Initializer
		
		Args:
			binResObj: (BinnedResultsStandard object) Note that this may get modified in place
			indicesFrom: (iter of ints) Contains the indices of atoms to get an rdf FROM (e.g. for g_{AB} indicesA should contain all indices of atom type A)
			indicesTo: (iter of ints) Contains the indices of atoms to get an rdf TO (e.g. for g_{AB} indicesB should contain all indices of atom type B)
			minDistAToB: (Bool) If False we do a normal rdf. If True, for every atom in group indicesFrom we only bin the SHORTEST distance to group indicesTo. Original use was to get fraction of oxygen atoms within a certain distance of Mg atoms.
			minDistVal: (float) The smallest distance value to use; setting > 0 allows us to use overlapping from/to indices without just finding the minimum distance is from atom X to itself

		"""
		self.binResObj = binResObj
		self.indicesFrom = indicesFrom
		self.indicesTo = indicesTo
		self.minDistAToB = minDistAToB
		self.minDistVal = minDistVal

	@property
	def primaryIndices(self):
		return self.indicesFrom


class CountNWithinDistOptions(calcDistrCoreHelp.CalcDistribOptionsBase):
	""" Options class to bin the number of atoms within certain distances (planar/horiznotal/total) 

	"""
	def __init__(self, binResObj, calcDistOpts, distRanges):
		""" Initializer
		
		Args:
			binResObj: (BinnedResultsStandard object) Note that this may get modified in place
			calcDistOpts: (CalcDistribOptionsBase) Options for how to calculate the distances (e.g. CalcHozDistOptions)
			distRanges: (len-2 float iter) We count the number of atoms, N, satisfying distRanges[0] < N <= distRanges[1] 
		
		NOTES:
			a) The bins in calcDistOpts MUST span the range distRanges, but no error is thrown if they dont

 
		"""
		self.binResObj = binResObj
		self.calcDistOpts = calcDistOpts
		self.distRanges = distRanges
		if self.distRanges[0]<=0:
			print("WARNING: lower distance range ={:.2f} may lead to self-counting".format(self.distRanges[0]))

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
		return _getOxyAndHyIndicesFromWaterIndicesAndGeom(waterIndices, inpGeom)


def _getOxyAndHyIndicesFromWaterIndicesAndGeom(waterIndices, inpGeom):
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



#Lots of parts stolen from water_rotations, which is now sorta redundant code
class WaterOrientationOptions(calcDistrCoreHelp.CalcDistribOptionsBase, _WaterOptsMixin):


	def __init__(self, binResObj, oxyIndices, hyIndices, angleType="roll", checkEdges=True, primaryIdxType="O"):
		""" Initializer
		
		Args:
			binResObj: (BinnedResultsStandard object) Note that this may get modified in place
			oxyIndices: (iter of ints) The oxygen indices for each water molecule
			hyIndices: (iter of len-2 ints) Same length as oxyIndices, but each contains the indices of two hydrogen indices bonded to the relevant oxygen
			angleType: (str - roll, pitch or azimuth) The angle type we want to bin. roll/pitch/azimuth correspond to rotations around standard x,y,z axes respectively
			checkEdges: (Bool) If True raise error if the edges of bins go beyond the domains of the angleType

		Raises:
			ValueError: If checkEdges=True and the bin edges go beyond angle domains
				 
		"""
		self.distribKey = "adf" #Angular distribution function. Not overly important what its set as
		self.domainTol = 1 #Allow bins to be up to 1 degrees outside the domain
		self.binResObj = binResObj
		self.checkEdges = checkEdges
		self.oxyIndices = oxyIndices
		self.hyIndices = hyIndices
		self._angleTypeToDomainMap = {"abs_roll":[0,90],"roll":[-90,90], "pitch":[-90,90], "azimuth":[-180,180]}
		self.angleType = angleType #IMPORTANT to do this after setting checkEdges
		self.primaryIdxType = "O"


	@classmethod
	def fromWaterIndicesAndGeom(cls, binResObj, waterIndices, inpGeom, primaryIdxType="O", checkEdges=True, angleType="roll"):
		oxyIndices, hyIndices = cls._getOxyAndHyIndicesFromWaterIndicesAndGeom(None,waterIndices, inpGeom)
		currArgs = [binResObj, oxyIndices, hyIndices]
		currKwargs = {"angleType":angleType, "checkEdges":checkEdges, "primaryIdxType":primaryIdxType}
		return cls(*currArgs, **currKwargs)


	@property
	def angleType(self):
		return self._angleType

	@angleType.setter
	def angleType(self, val):
		newDomain = self._angleTypeToDomainMap[val]
		self._checkBinEdgesWithinDomain(newDomain)
		self._angleType = val

	@property
	def domain(self):
		return self._angleTypeToDomainMap[self.angleType]

	def _checkBinEdgesWithinDomain(self, domain):
		if self.checkEdges:
			binResHelp._checkBinEdgesWithinDomain(self.binResObj, domain, self.domainTol)

	def __eq__(self,other):
		cmpAttrs = ["binResObj", "checkEdges","oxyIndices", "hyIndices", "angleType", "primaryIdxType"]
		for attr in cmpAttrs:
			valA, valB = getattr(self,attr), getattr(other,attr)
			if valA != valB:
				return False

		return True


class WaterMinDistPlusMinDistFilterOptions(calcDistrCoreHelp.CalcDistribOptionsBase, _WaterOptsMixin):
	""" Options for calculating the minimum distance between water and the "toIndices" which have minDist between filterDists for the atoms defined in filterIndices. Original purpose was to look at min(H-Mg) distance distribution for Mg atoms which did NOT have a water molecule chemisorbed """

	def __init__(self, binResObj, oxyIndices, hyIndices, toIndices, filterToIndices, filterDists, primaryIdxType="O", minDistType="all"):
		""" Initializer
		
		Args:
			binResObj: (BinnedResultsStandard object) Note that this may get modified in place
			oxyIndices: (iter of ints) The oxygen indices for each water molecule
			hyIndices: (iter of len-2 ints) Same length as oxyIndices, but each contains the indices of two hydrogen indices bonded to the relevant oxygen
			toIndices: (iter of ints) We calculate min-dist from water to these other indices
			filterIndices: (iter of ints) Indices of atoms we calculate minDist(toIndices[idxA]) from
			filterDists: (len-2 iter) [minDist, maxDist] for us to consider. Specifically we only consider values where min(filterDists) <= x < max(filterDists)
 
		"""
		self.binResObj = binResObj
		self.oxyIndices = oxyIndices
		self.hyIndices = hyIndices
		self.toIndices = toIndices
		self.filterToIndices = filterToIndices
		self.filterDists = filterDists
		self.primaryIdxType = primaryIdxType
		self.minDistType = minDistType


#TODO: Another one with an additional filter function based on minDists; that one will actually be used
class WaterMinDistOptions(calcDistrCoreHelp.CalcDistribOptionsBase, _WaterOptsMixin):

	def __init__(self, binResObj, oxyIndices, hyIndices, toIndices, primaryIdxType="O", minDistType="all", minVal=0.01):
		""" Initializer
		
		Args:
			binResObj: (BinnedResultsStandard object) Note that this may get modified in place
			oxyIndices: (iter of ints) The oxygen indices for each water molecule
			hyIndices: (iter of len-2 ints) Same length as oxyIndices, but each contains the indices of two hydrogen indices bonded to the relevant oxygen
			toIndices: (iter of ints) We calculate min-dist from water to these other indices
			primaryIdxType: (str) The element of the primary index. "O", "Ha" and "Hb" are the standard options
			minDistType: (str) Controls which atoms to get the minimum distance from. Current options are "all","o", and "h" (case insensitive)
			minVal: (float) The minimum value to take into account; this is to stop getting 0 when comparing two sets of overlapping groups
				 
		"""
		self.distribKey = "rdf"
		self.binResObj = binResObj
		self.oxyIndices = oxyIndices
		self.hyIndices = hyIndices
		self.toIndices = toIndices
		self.primaryIdxType = primaryIdxType
		self.minDistType = minDistType
		self.minVal = minVal

	@classmethod
	def fromWaterIndicesAndGeom(cls, binResObj, waterIndices, toIndices, inpGeom, primaryIdxType="O", minDistType="all"):
		""" Alternative initializer
		
		Args:
			binResObj: (BinnedResultsStandard object) Note that this may get modified in place
			waterIndices: (iter of len-3 int iters) Each element contains the indices of a water molecule
			toIndices: (iter of ints) We calculate min-dist from water to these other indices
			inpGeom: (plato_pylib UnitCell object) Used to figure out which indices correspond to Oxygen/Hydrogen
			primaryIdxType: (str) The element of the primary index. "O", "Ha" and "Hb" are the standard options
			minDistType: (str) Controls which atoms to get the minimum distance from. Current options are "all","o", and "h" (case insensitive)
 
		"""
		oxyIndices, hyIndices = cls._getOxyAndHyIndicesFromWaterIndicesAndGeom(None,waterIndices, inpGeom)
		currArgs = [binResObj, oxyIndices, hyIndices, toIndices]
		currKwargs = {"primaryIdxType":primaryIdxType, "minDistType":minDistType}
		return cls(*currArgs, **currKwargs)

	def __eq__(self,other):
		cmpAttrs = ["binResObj", "oxyIndices", "hyIndices", "toIndices", "primaryIdxType", "minDistType"]
		for attr in cmpAttrs:
			valA, valB = getattr(self,attr), getattr(other,attr)
			if valA != valB:
				return False

		return True


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
			minDistType: (str) Controls which atoms to get the minimum distance from. Current options are "all","o", and "h" (case insensitive)
 
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


class CountHBondsBetweenGenericGroupsOptions(calcDistrCoreHelp.CalcDistribOptionsBase):
	""" Options to count the number of hydrogen bonds between groups """

	def __init__(self, binResObj, fromNonHyIndices, fromHyIndices, toNonHyIndices, toHyIndices, acceptor=True, donor=True, maxOO=3.5, maxAngle=35, primaryIndices=None):
		""" Initializer
		
		Args:
			binResObj: (BinnedResultsStandard object) Note that this may get modified in place
			fromNonHyIndices: (iter of iter of ints) Each entry corresponds to an iter of indices for non-hydrogen atoms (h-bond acceptor atoms) on each molecule
			fromHyIndices: (iter of iter of ints) Each entry corresponds to an iter of indices for hydrogen atoms on each molecule
			toNonHyIndices: (iter of iter of ints) Same as fromNonHyIndices, except for molecules of the second group
			toHyIndices: (iter of iter of ints) Same as fromHyIndices, except for molecules of the second group
			acceptor: (Bool) If True we calculate dists/angles required to count number of groupA acceptors from groupB
			donor: (Bool) If True we calculate dists/angles required to count number of groupA donors to groupB
			maxOO: (float) The maximum X-X distance between two hydrogen-bonded water. For water X are the oxygen atoms; hence the variable name. Angles are only calculated when this criterion is fulfilled. NOTE: Should probably be changed almost always for non-water/hydroxyl cases
			maxAngle: (float) The maximum XA-XD-HD angle for a hydrogen bond; XA = acceptor non-hy, XD=Donor non-hy, HD=donor hydrogen
			primaryIndices: (iter of ints) Indices of atoms to associate these numbers with. Needed to slot into certain combined distributions. Default is to take [x[0] for x in self.fromNonHyIndices] (for water, this would be equivalent to taking the oxygen atom indices)

		NOTE:
			Don't have multiple NonHyIndices in one entry unless there are no hyIndices. For example, it would be fine to use both oxygen in CO2, but not for (HO)2-CO since theres no way to know which hydrogen is connected to each oxygen
	 
		"""
		self.distribKey = "nHbonds"
		self.binResObj = binResObj
		self.fromNonHyIndices = fromNonHyIndices
		self.fromHyIndices = fromHyIndices
		self.toNonHyIndices = toNonHyIndices
		self.toHyIndices = toHyIndices
		self.acceptor = acceptor
		self.donor = donor
		self.maxOO = maxOO
		self.maxAngle = maxAngle
		self._primaryIndices = primaryIndices

	#NOTE: NOT TESTED
	@property
	def primaryIndices(self):
		if self._primaryIndices is None:
			return [x[0] for x in self.fromNonHyIndices]
		else:
			return self._primaryIndices



class GetOODistsForHBondsBetweenGenericGroups(CountHBondsBetweenGenericGroupsOptions):
	""" Options for getting Oxy-Oxy (or more generally NonHy-NonHy) distances for every H-bond found;

	IMPORTANT: This bins 1 value PER H-BOND; rather than per primary index. Therefore these options cant be combined with (for example) the h-bond counting case though it wont neccesarily throw an error if you try """


	def __init__(self, *args,**kwargs):
		super().__init__(*args, **kwargs)
		self.distribKey = "hbondOODists"


class GetOOHozDistsForHBondsBetweenGenericGroups(CountHBondsBetweenGenericGroupsOptions):
	""" Options for geting Oxy-Oxy (or more generally NonHy-NonHy) horizontal distances for every H-bond found;

	IMPORTANT: This bins 1 value PER H-BOND; rather than per primary index. Therefore these options cant be combined with (for example) the h-bond counting case though it wont neccesarily throw an error if you try 
	"""

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.distribKey = "hbondOOHozDists"



class GetOOHAnglesForHBondsBetweenGenericGroups(CountHBondsBetweenGenericGroupsOptions):
	""" Options for getting O_A-O_D-H_D (or more generally NonHy_A-NonHy_D-Hy_D) angles for every H-bond found;

	IMPORTANT: This bins 1 value PER H-BOND; rather than per primary index. Therefore these options cant be combined with (for example) the h-bond counting case though it wont neccesarily throw an error if you try """

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.distribKey = "hbondOOHDists"


class GetNonHyToHyDistanceForHBonds(CountHBondsBetweenGenericGroupsOptions):

	def __init__(self, *args, **kwargs):
		super().__init__(*args,**kwargs)
		self.distribKey = "hbondNonHyToHyDists"


#
class GetAnglesForTriatomOpts(calcDistrCoreHelp.CalcDistribOptionsBase):
	""" Options for binning angles between triAtoms. Useful, for example, for looking at bond angles over a simulation

	"""

	def __init__(self, binResObj, triAtomIndices):
		""" Initializer
		
		Args:
			binResObj: (BinnedResultsStandard object) Note that this may get modified in place
			triAtomIndices: (iter of len-3 ints) Indices of atoms we want angles between. e.g. triAtomIndices = [ [1,2,3] ] will bin angles between atoms 1-2-3
				 
		"""
		self.distribKey = "triatom_angle"
		self.binResObj = binResObj
		self.triAtomIndices = triAtomIndices

	@property
	def primaryIndices(self):
		return [x[0] for x in self.triAtomIndices]

#
class GetDistsForDiatomOpts(calcDistrCoreHelp.CalcDistribOptionsBase):
	"""Options for binning distances between diatoms. Useful, for example, for looking at bondlengths

	IMPORTANT: This bins 1 value PER DIATOM rather than per primary index. Therefore these options should only really be combined with other "1-value per diatom" options, though you wont neccesarily hit an error if you dont
	"""
	
	def __init__(self, binResObj, diatomIndices):
		""" Initializer
		
		Args:
			binResObj: (BinnedResultsStandard object) Note that this may get modified in place
			diatomIndices: (iter of len-2 ints) Indices of diatoms (e.g. [ [0,1], [2,3] ])
				 
		"""
		self.distribKey = "diatom_dist"
		self.binResObj = binResObj
		self.diatomIndices = diatomIndices

	@property
	def primaryIndices(self):
		return [x[0] for x in self.diatomIndices]


class GetAngleWithGenericVectorForDiatomOpts(calcDistrCoreHelp.CalcDistribOptionsBase):
	"""Options for binning angles of diatoms with a generic vector (e.g. original use case is hydroxyl angle with surface normal)

	IMPORTANT: This bins 1 value PER DIATOM rather than per primary index. Therefore these options should only really be combined with other "1-value per diatom" options, though you wont neccesarily hit an error if you dont
	"""

	def __init__(self, binResObj, diatomIndices, inpVector, leftToRight=True):
		""" Initializer
		
		Args:
			binResObj: (BinnedResultsStandard object) Note that this may get modified in place
			diatomIndices: (iter of len-2 ints) Indices of diatoms (e.g. [ [0,1], [2,3] ])
			inpVector: (len-3 iter) Vector we take the angle with.	 
			leftToRight: (Bool) If True we generate vectors as diatom[1]-diatom[0]; if False we do diatom[0]-diatom[1]

		Notes:
			a) The angle we calculate should have domain [0,180]
			b) We calculate the angle after shifting the diatom vector to zero

		"""
		self.distribKey = "diatom_dist"
		self.binResObj = binResObj
		self.diatomIndices = diatomIndices
		self.inpVector = inpVector
		self.leftToRight = leftToRight

	@property
	def primaryIndices(self):
		return [x[0] for x in self.diatomIndices]



class GetHozDistsForDiatomOpts(calcDistrCoreHelp.CalcDistribOptionsBase):
	""" Options for binning horizontal distances between the two atoms of a diatom (e.g. between O and H in hydroxyl)

	IMPORTANT: This bins 1 value PER DIATOM rather than per primary index. Therefore these options should only really be combined with other "1-value per diatom" options, though you wont neccesarily hit an error if you dont

	"""

	def __init__(self, binResObj, diatomIndices):
		""" Initializer
		
		Args:
			binResObj: (BinnedResultsStandard object) Note that this may get modified in place
			diatomIndices: (iter of len-2 ints) Indices of diatoms (e.g. [ [0,1], [2,3] ])
				 
		"""
		self.distribKey = "diatom_dist"
		self.binResObj = binResObj
		self.diatomIndices = diatomIndices

	@property
	def primaryIndices(self):
		return [x[0] for x in self.diatomIndices]


class CountHBondsBetweenWaterGroupsOptions(calcDistrCoreHelp.CalcDistribOptionsBase):
	""" Options for counting number of hydrogen bonds between groups of water """

	def __init__(self, binResObj, fromOxyIndices, fromHyIndices, toOxyIndices, toHyIndices, acceptor=True, donor=True, maxOO=3.5, maxAngle=35):
		""" Initializer
		
		Args:
			binResObj: (BinnedResultsStandard object) Note that this may get modified in place
			fromOxyIndices: (iter of ints) The oxygen indices for each water molecule
			fromHyIndices: (iter of len-2 ints) Same length as oxyIndices, but each contains the indices of two hydrogen indices bonded to the relevant oxygen
			toOxyIndices: (iter of ints) The oxygen indices for each water molecule
			toHyIndices: (iter of len-2 ints) Same length as oxyIndices, but each contains the indices of two hydrogen indices bonded to the relevant oxygen
			acceptor: (Bool) If True we calculate dists/angles required to count number of groupA acceptors from groupB. Note changing the order of distFilterValues will give the reverse info (groupB acceptors from groupA)
			donor: (Bool) If True we calculate dists/angles required to count number of groupA donors to groupB. Note changing the order of distFilterValues will give the reverse info (groupB donors to groupA)
			maxOO: (float) The maximum O-O distance between two hydrogen-bonded water. Angles are only calculated when this criterion is fulfilled
			maxAngle: (float) The maximum OA-OD-HD angle for a hydrogen bond; OA = acceptor oxygen, OD=Donor oxygen, HD=donor hydrogen

		"""
		self.distribKey = "nHbonds"
		self.binResObj = binResObj
		self.fromOxyIndices = fromOxyIndices
		self.fromHyIndices = fromHyIndices
		self.toOxyIndices = toOxyIndices
		self.toHyIndices = toHyIndices
		self.donor = donor
		self.acceptor = acceptor
		self.maxOO = maxOO
		self.maxAngle = maxAngle

	@property
	def primaryIndices(self):
		return self.oxyIndices


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

