
import plato_pylib.shared.ucell_class as uCellHelp

from ..shared import cart_coord_utils as cartHelp

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

	def __init__(self, top=True, bottom=True, distTol=1e-1):
		""" Initializer
		
		Args:
			top: (Bool) Whether to include atoms on the bottom surface plane
			bottom:	 (Bool) Whether to include atoms on the top surface plane
			distTol: (float) Distance from outer surface plane an option can be 

		"""
		self.top = top
		self.bottom = bottom
		self.distTol = distTol

	def filterFunct(self, getIndicesInstance, inpGeom, inpIndices):
		#Get a cell with filtered coordinates
		lattParams, lattAngles = inpGeom.getLattParamsList(), inpGeom.getLattAnglesList()
		tempCell = uCellHelp.UnitCell(lattParams=lattParams, lattAngles=lattAngles)
		tempCell.fractCoords = [x for idx,x in enumerate(inpGeom.fractCoords) if idx in inpIndices] 

		#filter to atoms at the top/bottom surface
		outIndices = list()
		if self.top:
			currPlaneEqn = cartHelp.getPlaneEqnForOuterSurfaceAtoms(tempCell, top=True)
			topIndices = self._getFilteredIndicesForPlaneEqn(currPlaneEqn, tempCell, inpIndices)
			outIndices.extend(topIndices)

		if self.bottom:
			currPlaneEqn = cartHelp.getPlaneEqnForOuterSurfaceAtoms(tempCell, top=False)
			bottomIndices = self._getFilteredIndicesForPlaneEqn(currPlaneEqn, tempCell, inpIndices)
			outIndices.extend(bottomIndices)

		return sorted(list(set(outIndices)))

	def _getFilteredIndicesForPlaneEqn(self, planeEqn, inpGeom, inpIndices):
		distsFromPlane = [planeEqn.getDistanceOfPointFromPlane(x[:3]) for x in inpGeom.cartCoords]
		indicesInGeom = [idx for idx,dist in enumerate(distsFromPlane) if dist<self.distTol]
		outIndices = [inpIndices[idx] for idx in indicesInGeom]
		return outIndices


