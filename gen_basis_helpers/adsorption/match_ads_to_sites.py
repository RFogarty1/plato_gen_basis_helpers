
import math

from . import parse_from_geoms as parseFromGeomsBase
from ..shared import simple_vector_maths as vectHelp


class AdsObjIsH2FilterFunct(parseFromGeomsBase.MatchAdsObjsFilterFunct):

	def __init__(self, maxDist=2.0):
		self.maxDist = maxDist

	def objIsAdsorbedOnSite(self, sitePosition, surfOutwardsVector, adsObj):
		eleKeys = [x[-1].upper() for x in adsObj.geom]
		if len(eleKeys)>2:
			return False

		if not all([x=="H" for x in eleKeys]):
			return False

		dist = vectHelp.getDistTwoVectors(adsObj.geom[0][:3], adsObj.geom[1][:3])
		if dist > self.maxDist:
			return False

		return True

class AdsObjCentroidDistFromSiteFilterFunct(parseFromGeomsBase.MatchAdsObjsFilterFunct):

	def __init__(self, maxDist):
		self.maxDist = maxDist

	def objIsAdsorbedOnSite(self, sitePosition, surfOutwardsVector, adsObj):
		centroid = _getCentroidForCoordSet(adsObj.geom)
		dist = vectHelp.getDistTwoVectors(sitePosition, centroid)
		if dist>self.maxDist:
			return False

		return True

class AdsObjCentroidHorizontalDistFromSiteFilterFunct(parseFromGeomsBase.MatchAdsObjsFilterFunct):

	def __init__(self, maxDist):
		self.maxDist = maxDist

	def objIsAdsorbedOnSite(self, sitePosition, surfOutwardsVector, adsObj):
		centroid = _getCentroidForCoordSet(adsObj.geom)
		theta = vectHelp.getAngleTwoVectors( centroid, surfOutwardsVector )		
		distToAds = vectHelp.getDistTwoVectors(sitePosition, centroid)
		hozDist = distToAds*math.sin(math.radians(theta))
		if hozDist > self.maxDist:
			return False
		return True


class MatchSitesToH2AdsObjsStandard(parseFromGeomsBase.MatchAdsObjsToSitesBasedOnPositionsStandard):

	#Just setting the filter functions should be enough for this
	def __init__(self, maxDistTotal, maxDistHoz=1.0, maxHHBond=2.0):
		""" Initializer
		
		Args:
			maxDistTotal: (float) Maximum distance the CENTROID of H2 can be from a site and still be considered adsorbed
			maxDistHoz: (float) Maximum distance ORTHOGONAL to the surface vector that the H2 centroid can be while considered adsorbed
			maxHHBond: (float) Maximum distance the adsorbate H-H can be separated by

		"""
		self.checkH2Filter = AdsObjIsH2FilterFunct(maxDist=maxHHBond)
		self.maxDistFilter = AdsObjCentroidDistFromSiteFilterFunct(maxDistTotal)
		self.maxHozDistFilter = AdsObjCentroidHorizontalDistFromSiteFilterFunct(maxDistHoz)
		self.filterFuncts = [self.checkH2Filter, self.maxDistFilter, self.maxHozDistFilter]




def _getCentroidForCoordSet(inpCoords):
	centroid = [0,0,0]
	for coords in inpCoords:
		centroid = [x+c for x,c in zip(centroid,coords)]
	centroid = [x/len(inpCoords) for x in centroid]
	return centroid

