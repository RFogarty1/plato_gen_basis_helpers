
import itertools as it

class SurfaceAtomsFromInpGeom():
	""" Callable class for getting a surface object from unitCell. Main job is to filter the adsorbates out from the surface. Interface should be f(inpGeom)->surfAtoms; see self.getSurfaceObjFromInpGeom

	"""
	def getSurfaceAtomsFromInpGeom(self, inpGeom):
		""" Get co-ordinates for surface atoms 
		
		Args:
			inpGeom: (UnitCell object) Contains xyz co-ords
				 
		Returns
			surfCoords: (iter of len-4 iters) Cartesian co-ordinates ([x,y,z,symbol]) for surface atoms
	 
		"""
		raise NotImplementedError("")

	def __call__(self, inpGeom):
		return self.getSurfaceAtomsFromInpGeom(inpGeom)


class AdsorbatesFromInpGeom():
	""" Callable class for getting an iter of adsorbate objects from unitCell. Interface should be f(inpGeom)->adsorbateObjs. see self.getAdsorbateObjsFromInpGeom. 

	"""

	def getAdsorbateObjsFromInpGeom(self, inpGeom):
		""" Get an iter of adsorbate objects
		
		Args:
			inpGeom: (UnitCell object) Contains xyz co-ords
				 
		Returns
			outObjs: (iter of Adsorbate objects) 

		NOTE: This will return adsorbates where ANY atom is in the central cell, thus the number of returned values COULD be larger than the number of adsorbates per cell (and can fluctuate over a closed-system simulation)

		"""
		raise NotImplementedError("")

	def __call__(self, inpGeom, **kwargs):
		return self.getAdsorbateObjsFromInpGeom(inpGeom)



class MatchAdsObjsToSitesBasedOnPositionsStandard():
	"""Callable class to match adsorbate objects to a list of sitePositions based on applying a range of filter criteria. See self.getAdsObjsForSitePositions for the callable interface

	Attributes:
		filterFuncts: iter of functions with interface (sitePosition, surfOutwardsVector, adsObj)->Bool. Return True if adsObj consistent with being adsorbed on this site, false otherwise

	"""

	def getAdsObjsForSitePositions(self, sitePositions, surfOutwardsVector, adsObjs):
		""" Get the a list of adsorbates at each sitePosition based on filter criteria defined in this object
		
		Args:
			sitePositions: iter of iter of [x,y,z] (e.g. [ [1,2,3], [4,5,6] ]. Each represents x/y/z co-ordinates for an adsorption site
			surfOutwardsVector: len-3 iter. Vector which points outwards from a surface (one of the possible normal vectors to the surface)
			adsObjs: iter of Adsorbate objects
				 
		Returns
			adsorbedObjs: iter of iter of Adsorbate objects. 1 iter per element in sitePositions
	 
		"""
		outList = list()
		for sitePos in sitePositions:
			currAdsObjs = self._getAdsObjsForSingleSitePosition(sitePos, surfOutwardsVector, adsObjs)
			outList.append(currAdsObjs)
		return outList

	def _getAdsObjsForSingleSitePosition(self, sitePos, surfOutwardsVector, adsObjs):
		outList = list()
		for adsObj in adsObjs:
			addObj = True
			for filterFunct in self.filterFuncts:
				if filterFunct(sitePos, surfOutwardsVector, adsObj) is False:
					addObj = False
					break
			if addObj:
				outList.append(adsObj)
		return outList

	def __call__(self, sitePositions, surfOutwardsVector, adsObjs, **kwargs):
		return self.getAdsObjsForSitePositions(sitePositions, surfOutwardsVector, adsObjs)


class MatchAdsObjsFilterFunct():

	def objIsAdsorbedOnSite(self, sitePosition, surfOutwardsVector, adsObj):
		raise NotImplementedError("")

	def __call__(self, sitePosition, surfOutwardsVector, adsObj, **kwargs):
		return self.objIsAdsorbedOnSite(sitePosition, surfOutwardsVector, adsObj, **kwargs)

