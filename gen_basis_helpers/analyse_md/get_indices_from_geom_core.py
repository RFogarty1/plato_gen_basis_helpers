

import copy
import itertools as it
import plato_pylib.shared.ucell_class as uCellHelp

from . import get_neb_lists as nebListHelp
from ..shared import cart_coord_utils as cartHelp
from ..shared import plane_equations as planeEqnHelp

class GetSpecialIndicesFromInpGeomTemplate():

	def __init__(self, filterFuncts):
		""" Initializer
		
		Args:
			filterFuncts: iter of f(GetSpecialIndicesFromInpGeomTemplate, inpGeom, indices)->outIndices. These functions will be called one after another. They can communicate with each other by using self.scrachSpace
				 
		Returns
			indices: (iter of int) Indices for atoms of interest
	 
		"""

		self.filterFuncts = list(filterFuncts)


	def getIndicesFromInpGeom(self,inpGeom):
		outIndices = [x for x in range(len(inpGeom.cartCoords))]
		for funct in self.filterFuncts:
			outIndices = funct(self, inpGeom, outIndices)
		return outIndices


class FilterIndicesFunction():

	def filterFunct(self, getIndicesInstance, inpGeom, inpIndices):
		raise NotImplementedError("")

	def __call__(self, getIndicesInstance, inpGeom, inpIndices):
		return self.filterFunct(getIndicesInstance, inpGeom, inpIndices)


class FilterToExcludeElesNotInList(FilterIndicesFunction):
	""" Return only indices for atoms of selected element types (Case insensitive)"""

	def __init__(self, eleKeys):
		""" Initializer
		
		Args:
			eleKeys: (iter of str) Elements to filter down to

		"""
		self.eleKeys = list(eleKeys)

	def filterFunct(self, getIndicesInstance, inpGeom, inpIndices):
		coords = inpGeom.fractCoords
		outIndices = list()
		for idx in inpIndices:
			if coords[idx][-1].upper() in [a.upper() for a in self.eleKeys]:
				outIndices.append(idx)
		return outIndices


class FilterToOuterSurfaceAtoms(FilterIndicesFunction):
	""" Return only indices for atoms on the outer surface planes """

	def __init__(self, top=True, bottom=True, distTol=1e-1, nLayers=1):
		""" Initializer
		
		Args:
			top: (Bool) Whether to include atoms on the bottom surface plane
			bottom:	 (Bool) Whether to include atoms on the top surface plane
			distTol: (float) Distance from outer surface plane an atom can be 
			nLayers: (int) The number of layers to restrict to. e.g. 2 will return indices for the first TWO layers
		"""
		self.top = top
		self.bottom = bottom
		self.distTol = distTol
		self.nLayers = nLayers

	def filterFunct(self, getIndicesInstance, inpGeom, inpIndices):
		#filter to atoms at the top/bottom surface
		outIndices = list()
		if self.top:
			topIndices = self._getIndicesForOneSide(inpGeom, inpIndices, top=True)
			outIndices.extend(topIndices)

		if self.bottom:
			botIndices = self._getIndicesForOneSide(inpGeom, inpIndices, top=False)
			outIndices.extend(botIndices)

		return sorted(list(set(outIndices)))


	def _getIndicesForOneSide(self, inpGeom, inpIndices, top=True):
		lattParams, lattAngles = inpGeom.getLattParamsList(), inpGeom.getLattAnglesList()
		tempCell = uCellHelp.UnitCell(lattParams=lattParams, lattAngles=lattAngles)

		outIndices = list()
		nextIndices = inpIndices
		for unused in range(self.nLayers):
			tempCell.fractCoords = [x for idx,x in enumerate(inpGeom.fractCoords) if idx in nextIndices]
			currPlaneEqn = cartHelp.getPlaneEqnForOuterSurfaceAtoms(tempCell, top=top) 
			currIndices = self._getFilteredIndicesForPlaneEqn(currPlaneEqn, tempCell, nextIndices)
			nextIndices = [x for x in nextIndices if x not in currIndices]
			outIndices.extend(currIndices)
		return outIndices

	def _getFilteredIndicesForPlaneEqn(self, planeEqn, inpGeom, inpIndices):
		distsFromPlane = [planeEqn.getDistanceOfPointFromPlane(x[:3]) for x in inpGeom.cartCoords]
		indicesInGeom = [idx for idx,dist in enumerate(distsFromPlane) if dist<self.distTol]
		outIndices = [inpIndices[idx] for idx in indicesInGeom]
		return outIndices


class FilterToAtomsWithinDistanceOfSurfacePlane(FilterIndicesFunction):
	""" Return only indices for atoms within a certain distance from a given surface plane """


	def __init__(self, planeEqn, maxDist, distTol=1e-1, top=True, bottom=True, minImageConv=True):
		""" Initializer
		
		Args:
			planeEqn: (ThreeDimPlaneEquation) Defines the plane to use
			maxDist: (float) Will filter to atoms CLOSER to the plane than this (e.g. maxDist=1 means we filter to distance <1)
			top: (Bool) Whether to include atoms ABOVE the plane
			bottom: (Bool) Whether to include atoms BELOW the plane
			distTol: (float) If a distance is <= to this value of the plane it is considered "on the plane" and that index will be always be returned regardless of the top/bottom parameter values (e.g. if bottom=False and an index was found at -0.5*distTol it would be returned)
			minImageConv: (Bool) Whether to take periodic boundaries into account using the minimum image convention. Setting to False is probably faster, and may be fine if your sure PBCs wont matter (e.g. your surface plane is already at the centre of the cell)
				 
		"""
		self.planeEqn = planeEqn
		self.maxDist = maxDist
		self.distTol = distTol
		self.top = top
		self.bottom = bottom
		self.minImageConv = minImageConv

	def filterFunct(self, getIndicesInstance, inpGeom, inpIndices):
		useIndices = sorted(inpIndices)
		relCoords = [x[:3] for idx,x in enumerate(inpGeom.cartCoords) if idx in useIndices]

		if self.minImageConv:
			signedDists = cartHelp.getDistancesOfAtomsFromPlaneEquation_nearestImageAware(inpGeom, self.planeEqn, useIndices, signed=True)
		else:
			signedDists = [self.planeEqn.getSignedDistanceOfPointFromPlane(x) for x in relCoords]

		#Figure out which indices meet the criterion
		outIndicesUnmapped = list()
		for idx,sDist in enumerate(signedDists):
			absDist = abs(sDist)
			if absDist<self.distTol:
				outIndicesUnmapped.append(idx)
			else:
				if self.top:
					if (sDist>0) and (absDist<self.maxDist):
						outIndicesUnmapped.append(idx)
				if self.bottom:
					if (sDist<0) and (absDist<self.maxDist):
						outIndicesUnmapped.append(idx)

		#Map back to the original indices (from inpIndices)
		outIndices = [useIndices[idx] for idx in outIndicesUnmapped]

		return outIndices


class FilterToExcludeIndicesWithoutNebsAmongstRemaning(FilterIndicesFunction):
	""" Returns only the indices for atoms which have neighbours amongst inpIndices """

	def __init__(self, maxDist, restrictToPairs=None):
		""" Initializer
		
		Args:
			maxDist: (float) Maximum distance between atoms for them to be considered neighbours
			restrictToPairs: (Optional, iter of len-2 iters) If this is set then only include these pairs of elements in the neighbour lists. E.g. if restrictToPairs = [ ["Mg","O"] ] then an Mg atom with only Mg neighbours will be filtered out
				 
		"""
		self.maxDist = maxDist
		self.restrictToPairs = restrictToPairs

	def filterFunct(self, getIndicesInstance, inpGeom, inpIndices):

		useInpIndices = sorted(inpIndices) #These NEED to be in order for us to have simple mapping between used/input indices

		if self.restrictToPairs is not None:
			sortedRestrictions = [sorted(x) for x in self.restrictToPairs]

		cellForNebs = copy.deepcopy(inpGeom)
		nebCellCartCoords = [x for idx,x in enumerate(inpGeom.cartCoords) if idx in useInpIndices]
		cellForNebs.cartCoords = nebCellCartCoords

		nebLists = nebListHelp.getNeighbourListsForInpCell_imagesMappedToCentral(cellForNebs, self.maxDist)

		outIndices = list()
		for idx,nList in enumerate(nebLists):
			if len(nList)>0:
				if self.restrictToPairs is None:
					outIndices.append( useInpIndices[idx] )
				else:
					pairLists = [ sorted([nebCellCartCoords[idx][-1], nebCellCartCoords[x][-1]]) for x in nList ]
					filteredList = [x for x in pairLists if x in sortedRestrictions]
					if len(filteredList)>0:
						outIndices.append( useInpIndices[idx] )

		return sorted(outIndices)


class FilterToExcludeIndicesFurtherOutOfPlaneThanCutoff(FilterIndicesFunction):
	""" Returns only the indices for atoms with out-of-plane distances within a cutoff of a list of input points.

	Original use case was to get indices of atoms with adsorbates above them (which meant a horizontal distance was needed rather than the total distance) """

	def __init__(self, maxDist, planeEqn, inpPoints):
		""" Initializer
		
		Args:
			maxDist: (float)
			planeEqn: (ThreeDimPlaneEquation) Defines the plane to use
			inpPoints: (iter of len-3 iters) out-of-plane distance will be determined from these points
 
		"""
		self.maxDist = maxDist
		self.planeEqn = planeEqn
		self.inpPoints = inpPoints

	def filterFunct(self, getIndicesInstance, inpGeom, inpIndices):
		relCoords = [x[:3] for idx,x in enumerate(inpGeom.cartCoords) if idx in inpIndices]

		outIndices = list()
		for idx,coord in enumerate(relCoords):
			currDists = cartHelp.getDistancesFromPointAlongPlane(coord, self.inpPoints, self.planeEqn)
			if any([x<self.maxDist for x in currDists]):
				outIndices.append( inpIndices[idx] )

		return outIndices


class FilterToExcludeIndicesBasedOnNumberOfAtomsInSurfacePlane(FilterIndicesFunction):
	""" Returns only the indices for atoms with certain number of atoms in a layer.

	Original use case was to find dissolved Mg atoms (which sit on their own away from the surface """

	def __init__(self, minAtomsInPlane, maxAtomsInPlane, planeTol=5e-1, restrictNebsToInpIndices=True):
		""" Initializer
		
		Args:
			maxAtomsInPlane: (int) Filter to include atoms which are part of a plane with <=maxAtomsInPlane
			minAtomsInPlane: (int) Filter include atoms which are part of a plane with >=minAtomsInPlane (1 would filter None out) 
			planeTol: (float) Maximum distance an atom can be from the top/botom of a plane, and still be considered as "in" that plane (can maybe be thought of as planeHeight)
			restrictNebsToInpIndices: (Bool) If True, number of atoms in a plane can only include the input indices (e.g. for an MgO system if we already filtered out the oxygen atoms and restrictNebsToInpIndices=True then only Mg atoms will contribute to minAtomsInPlane/maxAtomsInPlane)

		NOTE:
			Both maxAtomsInPlane and minAtomsInPlane are INCLUSIVE of the atom which we define the plane around. i.e. there is always AT LEAST one atom in a plane

			
		"""
		self.minAtomsInPlane = minAtomsInPlane
		self.maxAtomsInPlane = maxAtomsInPlane
		self.planeTol = planeTol
		self.restrictNebsToInpIndices = restrictNebsToInpIndices


	def filterFunct(self, getIndicesInstance, inpGeom, inpIndices):
		self._checkMinMaxConsistent()
		indicesInSamePlane = self.getIndicesInSameSurfPlaneForEachAtom(inpGeom, inpIndices, includeSelf=True)
		outIndices = list()
		for indicesInPlane in indicesInSamePlane:
			currIdx = indicesInPlane[0]
			nInPlane = len(indicesInPlane[1])
			if (nInPlane>=self.minAtomsInPlane) and (nInPlane<=self.maxAtomsInPlane):
				outIndices.append(currIdx)
		return outIndices

	def _checkMinMaxConsistent(self):
		if self.minAtomsInPlane > self.maxAtomsInPlane:
			raise ValueError(" minAtomsInPlane > maxAtomsInPlane ({}>{})".format(self.minAtomsInPlane,self.maxAtomsInPlane))

	#TODO: Dont restrict the neighbours to inpIndices; thats likely gonna be the more general use case regardless. Its easier to filter out later anyway
	def getIndicesInSameSurfPlaneForEachAtom(self, inpGeom, inpIndices, includeSelf=False):
		outVals = list()
		planeEqn = cartHelp.getPlaneEqnForOuterSurfaceAtoms(inpGeom)
		for idx in inpIndices:
			indicesInSamePlane = self._getIndicesInSameSurfPlaneOneAtomIdx(inpGeom, idx, planeEqn, includeSelf=includeSelf)
			filteredToOnlyInInpIndices = [x for x in indicesInSamePlane if x in inpIndices]
			if self.restrictNebsToInpIndices:
				outVals.append( [idx,filteredToOnlyInInpIndices] )
			else:
				outVals.append( [idx,indicesInSamePlane] )
		return outVals

	#TODO: This will likely need factoring out a some point
	def _getIndicesInSameSurfPlaneOneAtomIdx(self, inpGeom, inpIdx, planeEqn, includeSelf=False):
		#Put our current atom right at the centre of the cell; meaning we shouldnt have an issue with periodic boundaries
		tempGeom = copy.deepcopy(inpGeom)
		startFract = inpGeom.fractCoords[inpIdx][:3]
		fractTranslate = [x-y for x,y in it.zip_longest([0.5,0.5,0.5],startFract)]
		uCellHelp.applyTranslationVectorToFractionalCoords(tempGeom, fractTranslate, foldInAfter=True)

		#recentre our plane equation using the cart coords
		startCart = tempGeom.cartCoords[inpIdx]
		startDVal = planeEqn.coeffs[-1]
		useDVal = planeEqn.calcDForInpXyz(startCart[:3])
		planeEqn.coeffs = planeEqn.coeffs[:3] + [useDVal]


		#Get all the indices which are within planeTol
		outIndices = list()
		for idx,coord in enumerate(tempGeom.cartCoords):
			if idx != inpIdx:
				currDist = planeEqn.getDistanceOfPointFromPlane(coord[:3])
				if currDist < self.planeTol:
					outIndices.append(idx)
			else:
				if includeSelf:
					outIndices.append(idx)

		#Reset the planeEqn to what it started as
		planeEqn.coeffs = planeEqn.coeffs[:3] + [startDVal]

		return outIndices


class FilterToExcludeAtomsOutsideCutoffDistFromIndices(FilterIndicesFunction):
	""" Excludes indices further than a cutoff for any inpIndices. The indices of these needs to be chosen at initiation time

	"""

	def __init__(self, cutoffDist, cutoffIndices):
		""" Initializer
		
		Args:
			cutoffDist: (float) Maximum distance from inpIndices (in filterFunct) to any of atoms in cutoffIndices
			cutoffIndices: (iter of ints) Indices of atoms we're testing for being within cutoff dist
				 
		"""
		self.cutoffDist = cutoffDist
		self.cutoffIndices = sorted(cutoffIndices)

	def filterFunct(self, getIndicesInstance, inpGeom, inpIndices):

		#Figure out the indices we care about
		useIndices = sorted(inpIndices) #Safer to have them sorted
		outIndices = [x for x in useIndices if x in self.cutoffIndices]
		allIndicesForCell = self.cutoffIndices + [idx for idx in inpIndices if idx not in self.cutoffIndices]

		#Create at temporary cell with only the indices we care about + get neighbour lists (using cutoffDist)
		startCellCartCoords = inpGeom.cartCoords
		cellForNebs = copy.deepcopy(inpGeom)

		nebCellCartCoords = list()
		for idx in allIndicesForCell:
			currCoords = startCellCartCoords[idx]
			nebCellCartCoords.append(currCoords)

		cellForNebs.cartCoords = nebCellCartCoords
		nebLists = nebListHelp.getNeighbourListsForInpCell_imagesMappedToCentral(cellForNebs, self.cutoffDist)

		#Filter out inpIndices which dont have neighbours amongst self.cutoffIndices
		startIdx = len(self.cutoffIndices)
		unMappedOutIndices = list()
		for idx, nebList in enumerate(nebLists[startIdx:]):
			if len(nebList)>0:
				#Check if any of these indices matches 
				if any([x<startIdx for x in nebList]):
					unMappedOutIndices.append(idx)

		#Map back to original indices
		indicesForCellWithoutCutoffIndices = [x for x in allIndicesForCell if x not in self.cutoffIndices]
		mappedIndices = [ indicesForCellWithoutCutoffIndices[idx] for idx in unMappedOutIndices ]

		return sorted(outIndices + mappedIndices)


