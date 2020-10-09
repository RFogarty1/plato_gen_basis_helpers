

from . import add_adsorbates as addAdsorbsHelp
from . import build_adsorbate_layers as adsorbateLayerBuilder
from . import adsorbate_rep_objs as adsRepObjs
from . import site_occupiers as siteOccHelp

class GetSurfObjForOneFractCoverageSimpleSingleTypeAdsorbates():
	""" Callable class which generates a surface object for a specified fractional coverage of a single type adsorbate

	"""
	def __init__(self, startSurface, surfaceToSites, adsorbateObj, distance, siteOccStrat=None):
		""" Initializer
		
		Args:
			startSurface: (BaseSurface instance) The surface without adsorbates
			surfaceToSites: (BaseSurfaceToSites instance) Used to get positions of surface sites from an input surface
			adsorbateObj: (Adsorbate instance) Single object, just contains .geom for the adsorbate
			distance: (float) Distance of the adsorbate from the surface
			siteOccStrat: (Optional, BaseAdsorbateSiteOccupierStrategy) Callable used to decide which sites to occupy for partial coverages. Default is to just occupy in the order sites are listed; this is likely only useful if only one configuration is possible (e.g. full coverage or minimal coverage)

		"""
		self.startSurface = startSurface
		self.surfaceToSites = surfaceToSites
		self.adsorbateObj = adsorbateObj
		self.siteOccStrat = siteOccStrat if siteOccStrat is not None else siteOccHelp.OccupyInOrderSiteOccupier()
		self.distance = distance
	
	def _getAdsAdder(self,fractCoverage):
		args = self.adsorbateObj, fractCoverage
		kwargs = {"addStrat":self.siteOccStrat, "distance":self.distance}
		return addAdsorbsHelp.SingleTypeGetAdsorbatesForSites(*args,**kwargs)
	
	def __call__(self,fractCoverage):
		addAdsorbs = self._getAdsAdder(fractCoverage)
		layerBuilder = adsorbateLayerBuilder.AdsorbateLayerBuilder(self.surfaceToSites, addAdsorbs)
		outLayer = layerBuilder.build(self.startSurface)
		outObj = adsRepObjs.SurfaceWithAdsorbatesStandard(self.startSurface, [outLayer])
		return outObj


