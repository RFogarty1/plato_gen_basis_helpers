
import copy
import itertools as it

import plato_pylib.shared.ucell_class as uCellHelp
import plato_pylib.shared.unit_convs as uConvHelp

from . import get_indices_from_geom_core as getIdxCore
from . import get_neb_lists as nebListHelp
from . import calc_dists as calcDistHelp

from ..shared import cart_coord_utils as cartHelp

class GetSurfaceIndicesFromGeomStandard(getIdxCore.GetSpecialIndicesFromInpGeomTemplate):

	def __init__(self, surfEles, top=True, bottom=True, distTol=1e-1):
		self.surfEles = list(surfEles)
		self.top = top
		self.bottom = bottom
		self.distTol = distTol

	@property
	def filterFuncts(self):
		eleFilter = getIdxCore.FilterToExcludeElesNotInList(self.surfEles)
		surfFilter = getIdxCore.FilterToOuterSurfaceAtoms(top=self.top, bottom=self.bottom, distTol=self.distTol)
		return [eleFilter, surfFilter]


class GetTopWaterLayerIndices():

	def __init__(self, waterDetector=None, maxHeightLayer=uConvHelp.ANG_TO_BOHR, top=True):
		""" Initializer
		
		Args:
			waterDetector: (GetWaterMoleculeIndicesFromGeomStandard object) Used to get indices related to water
			maxLayerHeight: (float) Difference in OXYGEN atom heights between top/bottom of the layer
			top: (Bool) If True then look for the O at the top of the surface, if False look for the O at the bottom of the surface
				 
		Returns
			outIndices: (iter of len-3 iters) The indices corresponding to the top(or bottom) layer of water molecules
	 
		"""
		self.waterDetector = GetWaterMoleculeIndicesFromGeomStandard() if waterDetector is None else waterDetector
		self.maxHeightLayer = maxHeightLayer
		self.top = top

	def getIndicesFromInpGeom(self, inpGeom):
		waterIndices = self.waterDetector.getIndicesFromInpGeom(inpGeom)
		allFractCoords = inpGeom.fractCoords
		allWaterAtomIndices = [x for x in it.chain(*waterIndices)]
		relFractCoords = [x for idx,x in enumerate(allFractCoords) if idx in allWaterAtomIndices and x[-1].upper()=="O"]

		#TODO: Probably have an "optimisation" option to skip this bit; the purpose of which is mainly to make sure
		#we deal with PBCs properly
		useCell = copy.deepcopy(inpGeom)
		useCell.fractCoords = relFractCoords
		uCellHelp.foldAtomicPositionsIntoCell(useCell)

		cartCoords = useCell.cartCoords

		planeEqn = cartHelp.getPlaneEqnForOuterSurfaceAtoms(useCell, top=self.top)
		distsFromPlane = [ planeEqn.getDistanceOfPointFromPlane(x[:3]) for x in cartCoords]
		outWaterIndices = list()
		for idx,val in enumerate(distsFromPlane):
			if val < self.maxHeightLayer:
				outWaterIndices.append( waterIndices[idx] )

		return outWaterIndices

class GetWaterMoleculeIndicesFromGeomStandard():

	def __init__(self, minOH=0.01, maxOH=1.2*uConvHelp.ANG_TO_BOHR, minAngle=100, maxAngle=110):
		self.minOH = minOH
		self.maxOH = maxOH
		self.minAngle = minAngle
		self.maxAngle = maxAngle

	def getIndicesFromInpGeom(self, inpGeom):
		cartCoords = inpGeom.cartCoords
		allAtomsNebLists = nebListHelp.getNeighbourListsForInpCell_imagesMappedToCentral(inpGeom, self.maxOH)
		oxyIndices = [idx for idx,unused in enumerate(cartCoords) if cartCoords[idx][-1].upper()=="O"]

		outIndices = list()
		for oIdx in oxyIndices:
			currNebList = allAtomsNebLists[oIdx]
			if len(currNebList)==2:
				eleTypes = [cartCoords[x][-1] for x in currNebList]
				if all([x.upper()=="H" for x in eleTypes]):
					oCoord, hCoordA, hCoordB = [ cartCoords[idx][:3] for idx in [oIdx]+currNebList ]
					hDistA = calcDistHelp.calcSingleDistBetweenCoords_minImageConv(inpGeom, oCoord, hCoordA)
					hDistB = calcDistHelp.calcSingleDistBetweenCoords_minImageConv(inpGeom, oCoord, hCoordB)
					angle = calcDistHelp.calcSingleAngleBetweenCoords_minImageConv(inpGeom, hCoordA, oCoord, hCoordB)
					if self._areDistsAndAnglesConsistent([hDistA,hDistB], angle):
						outIndices.append( [oIdx] + currNebList )

		return outIndices

	def _areDistsAndAnglesConsistent(self, dists, angle):
		if any([x<self.minOH for x in dists]):
			return False

		if any([x>self.maxOH for x in dists]):
			return False

		if angle<self.minAngle:
			return False

		if angle>self.maxAngle:
			return False

		return True




