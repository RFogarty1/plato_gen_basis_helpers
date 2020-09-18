
from . import adsorbate_rep_objs as adsObjs

class AdsorbateLayerBuilder():

	def __init__(self, surfaceToSites, addAdsorbatesToSites):
		""" 
		
		Args:
			surfaceToSites: (BaseSurfaceToSites obj) Helps us get the adsorbate site positions 
			addAdsorbatesToSites: (BaseGetAdsorbatesForSites obj) Helps us figure out where to put adsorbates
				 
		"""
		self.surfaceToSites = surfaceToSites
		self.addAdsorbatesToSites = addAdsorbatesToSites

	def build(self, inpSurface):
		""" Generates a AdsorbateLayer object when passed a surface object
		
		Args:
			inpSurface: (BaseSurface object) 
				 
		Returns
			outLayer: (AdsorbateLayer object) Represented a layer of adsorbed material
	 
		"""
		sitePositions = self.surfaceToSites.getSurfaceSitesFromInpSurface(inpSurface)
		surfVector = self.surfaceToSites.getOutwardsSurfaceVectorFromSurface(inpSurface)
		adsorbateList = self.addAdsorbatesToSites.getAdsorbateListForInpSites(sitePositions, surfObj=inpSurface, surfaceVector=surfVector)
		distances = self.addAdsorbatesToSites.getDistances(sitePositions, surfObj=inpSurface, surfaceVector=surfVector)
		outLayer = adsObjs.AdsorbateLayer(sitePositions, adsorbateList, distances, surfVector)
		return outLayer


