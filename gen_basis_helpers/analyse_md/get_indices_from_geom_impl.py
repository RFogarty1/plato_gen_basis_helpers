
import plato_pylib.shared.unit_convs as uConvHelp

from . import get_indices_from_geom_core as getIdxCore
from . import get_neb_lists as nebListHelp
from . import calc_dists as calcDistHelp

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



class GetWaterMoleculeIndicesFromGeomStandard():

	def __init__(self, minOH=0.01, maxOH=1.2*uConvHelp.ANG_TO_BOHR, minAngle=105, maxAngle=115):
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
					angle = calcDistHelp.calcSingleAngleBetweenCoords_minImageConv(inpGeom, oCoord, hCoordA, hCoordB)
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




