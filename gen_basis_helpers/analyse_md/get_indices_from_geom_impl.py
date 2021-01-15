
from . import get_indices_from_geom_core as getIdxCore


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


