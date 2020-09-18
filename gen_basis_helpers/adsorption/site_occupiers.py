
import itertools as it

from ..shared import simple_vector_maths as vectHelp

class BaseAdsorbateSiteOccupierStrategy():

	def matchAdsorbatesToSites(self, sitePositions, adsorbates, surfaceVector=None, distances=None, surfCell=None):
		""" (Base class docstring) Gets a list of adsorbates corresponding to list of sitePositions. None is used to denote an empty site
		
		Args:
			sitePositions: (iter of len 3 iters) Each element contains x,y,z co-ordinates for a site
			adsorbates: (iter of AdsorbateBase objects) Adsorbates in the order they are to be added (order may or may not matter dependending on strategt used)
			surfaceVector (Optional, some strats may use): (len-3 iter) vector points out from the surface
			distances (Optional, some strats may use): (iter of len-3 iters) The distance of each adsorbate from the relevant site
			surfCell (Optional, some strats may use): (BaseSurface obj) 
			
		Returns:
			occSites: (iter of adsorbates) Same length as sitePositions. Each element contains either None (unoccupied site) or an adsorbate object. The position of each index (idx) corresponds to sitePositions[idx]
	 
		"""
		raise NotImplementedError("")

	def __call__(self, *args, **kwargs):
		return self.matchAdsorbatesToSites(*args, **kwargs)


class OccupyInOrderSiteOccupier(BaseAdsorbateSiteOccupierStrategy):

	def __init__(self):
		pass

	def matchAdsorbatesToSites(self, sitePositions, adsorbates, **kwargs):
		outAdsorbates = [None for x in sitePositions]
		for idx, obj in enumerate(adsorbates):
			outAdsorbates[idx] = obj	
		return outAdsorbates


class OccupyClosestIgnoringAdsorbateGeomsSiteOccupier(BaseAdsorbateSiteOccupierStrategy):
	
	def __init__(self):
		pass

	def _findIdxOfNearestUnoccPosition(self, sitePositions, outAdsorbates):
		occPositions = [pos for pos,ads in it.zip_longest(sitePositions,outAdsorbates) if ads is not None]
		unOccIndices =  [idx for idx,(pos,ads) in enumerate(it.zip_longest(sitePositions,outAdsorbates)) if ads is None]
		nearestIdx = unOccIndices[0]
		nearestDist = min( [  vectHelp.getDistTwoVectors(x,sitePositions[nearestIdx]) for x in occPositions] )

		for idx in unOccIndices[1:]:
			currNearest = min( [  vectHelp.getDistTwoVectors(x,sitePositions[idx]) for x in occPositions] )
			if (currNearest < nearestDist):
				nearestDist, nearestIdx = currNearest, idx

		return nearestIdx


	def matchAdsorbatesToSites(self, sitePositions, adsorbates, **kwargs):
		outAdsorbates = [None for x in sitePositions]
	
		for idx,obj in enumerate(adsorbates):
			if idx==0:
				outAdsorbates[0] = obj #First adsorbate is added to first site always
			else:
				toOccIdx = self._findIdxOfNearestUnoccPosition(sitePositions, outAdsorbates)
				outAdsorbates[toOccIdx] = obj

		return outAdsorbates


class MinimizeAverageOfNonPeriodicOccSiteDistsStandardSiteOccupier(BaseAdsorbateSiteOccupierStrategy):

	def __init__(self):
		pass

	def _getIdxForNextAdsorptionSite(self, sitePositions, outAdsorbates):
		occPositions = [pos for pos,ads in it.zip_longest(sitePositions,outAdsorbates) if ads is not None]
		unOccIndices =  [idx for idx,(pos,ads) in enumerate(it.zip_longest(sitePositions,outAdsorbates)) if ads is None]

		outIdx = unOccIndices[0]
		minAvgDist = self._getAvDistanceOfPointAToOtherPoints( sitePositions[outIdx], occPositions )

		for idx in unOccIndices[1:]:
			currVal = self._getAvDistanceOfPointAToOtherPoints( sitePositions[idx], occPositions )
			if currVal < minAvgDist:
				outIdx=idx
				minAvgDist = currVal

		return outIdx

	def _getAvDistanceOfPointAToOtherPoints(self, pointA, otherPoints):
		allDists = [vectHelp.getDistTwoVectors(pointA, x) for x in otherPoints]
		avDist = sum(allDists)/len(allDists)
		return avDist

	#NOTE: Factor out the duplication from OccupyClosestIgnoringAdsorbateGeomsSiteOccupier
	def matchAdsorbatesToSites(self, sitePositions, adsorbates, **kwargs):
		outAdsorbates = [None for x in sitePositions]

		for idx, obj in enumerate(adsorbates):
			if idx==0:
				outAdsorbates[0] = obj #Add the first adsorbate to first site
			else:
				toOccIdx = self._getIdxForNextAdsorptionSite(sitePositions, outAdsorbates)
				outAdsorbates[toOccIdx] = obj

		return outAdsorbates
