

import copy
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


	def __init__(self, planeEqn, maxDist, distTol=1e-1, top=True, bottom=True):
		""" Initializer
		
		Args:
			planeEqn: (ThreeDimPlaneEquation) Defines the plane to use
			maxDist: (float) Will filter to atoms CLOSER to the plane than this (e.g. distTol=1 means we filter to distance <1)
			top: (Bool) Whether to include atoms ABOVE the plane
			bottom: (Bool) Whether to include atoms BELOW the plane
			distTol: (float) If a distance is <= to this value of the plane it is considered "on the plane" and that index will be always be returned regardless of the top/bottom parameter values (e.g. if bottom=False and an index was found at -0.5*distTol it would be returned)
				 
		"""
		self.planeEqn = planeEqn
		self.maxDist = maxDist
		self.distTol = distTol
		self.top = top
		self.bottom = bottom

	def filterFunct(self, getIndicesInstance, inpGeom, inpIndices):
		relCoords = [x[:3] for idx,x in enumerate(inpGeom.cartCoords) if idx in inpIndices]
		signedDists = [self.planeEqn.getSignedDistanceOfPointFromPlane(x) for x in relCoords]

		outIndices = list()
		for idx,sDist in enumerate(signedDists):
			absDist = abs(sDist)
			if absDist<self.distTol:
				outIndices.append(idx)
			else:
				if self.top:
					if (sDist>0) and (absDist<self.maxDist):
						outIndices.append(idx)
				if self.bottom:
					if (sDist<0) and (absDist<self.maxDist):
						outIndices.append(idx)

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

		if self.restrictToPairs is not None:
			sortedRestrictions = [sorted(x) for x in self.restrictToPairs]

		cellForNebs = copy.deepcopy(inpGeom)
		nebCellCartCoords = [x for idx,x in enumerate(inpGeom.cartCoords) if idx in inpIndices]
		cellForNebs.cartCoords = nebCellCartCoords

		nebLists = nebListHelp.getNeighbourListsForInpCell_imagesMappedToCentral(cellForNebs, self.maxDist)

		outIndices = list()
		for idx,nList in enumerate(nebLists):
			if len(nList)>0:
				if self.restrictToPairs is None:
					outIndices.append( inpIndices[idx] )
				else:
					pairLists = [ sorted([nebCellCartCoords[idx][-1], nebCellCartCoords[x][-1]]) for x in nList ]
					filteredList = [x for x in pairLists if x in sortedRestrictions]
					if len(filteredList)>0:
						outIndices.append( inpIndices[idx] )

		return outIndices


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
			currDists = [ planeEqnHelp.getOutOfPlaneDistTwoPoints(coord,x,self.planeEqn) for x in self.inpPoints]
			if any([x<self.maxDist for x in currDists]):
				outIndices.append( inpIndices[idx] )


		return outIndices


