
import plato_pylib.shared.unit_convs as uConvHelp

from . import parse_from_geoms as parseFromGeomBase
from . import water_adsorbate as waterAdsHelp


class DetectOuterAdsorbedWaterLayer(parseFromGeomBase.AdsorbatesFromInpGeom):
	""" This essentially just wraps the code for finding the indices of water molecules """

	def __init__(self, waterIdxDetector):
		""" Initializer
		
		Args:
			waterIdxDetector: (GetTopWaterLayerIndices object) This detects indices of water molecules
		"""
		self.waterIdxDetector = waterIdxDetector


	def getAdsorbateObjsFromInpGeom(self, inpGeom):
		waterIndices = self.waterIdxDetector.getIndicesFromInpGeom(inpGeom)
		waterAdsObjs = waterAdsHelp.getWaterAdsorptionObjsFromInpCellAndWaterIndices(inpGeom, waterIndices)
		return waterAdsObjs

