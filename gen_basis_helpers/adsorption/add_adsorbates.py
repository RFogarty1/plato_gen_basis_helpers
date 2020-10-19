
import copy
import numpy as np

from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra

from . import site_occupiers as addStratHelp
from ..shared import cart_coord_utils as cartHelp
from ..shared import simple_vector_maths as vectHelp

class BaseGetAdsorbatesForSites():
	""" (Base docstring) Callable class, when passed a list of positions for adsorption sites it figures out which sites are occupied by which adsorbate. See getAdsorbateListForInpSites for callable interface

	"""
	def getAdsorbateListForInpSites(self, sitePositions, surfObj=None, surfaceVector=None):
		""" When given a list of site positions, this returns a mixed list of adsorbate objects (in indices corresponding to occupied sites) and None (for indices corresponding to unoccupied sites) 
		
		Args:
			sitePositions: (iter of len-3 iters) Each element is [x,y,z] of one adsorption site
			surfaceVector (Optional, some strats may use): (len-3 iter) vector points out from the surface
			surfCell (Optional, some strats may use): (BaseSurface obj)

		Returns
			occSites: (iter of adsorbates) Same length as sitePositions. Each element contains either None (unoccupied site) or an adsorbate object. The position of each index (idx) corresponds to sitePositions[idx]
	 
		"""
		raise NotImplementedError("")

	def getDistances(self, sitePositions, surfObj=None, surfaceVector=None):
		""" List of distances of adsorbates from adsorption site. Follow same logic as return value of __call__ (None for unocc sites)
		"""
		raise NotImplementedError("")

	def __call__(self, inpSites, surfObj=None, surfaceVector=None):
		return self.getAdsorbateListForInpSites(inpSites, surfObj=surfObj, surfaceVector=surfaceVector)


class SingleTypeGetAdsorbatesForSites(BaseGetAdsorbatesForSites):
	""" Callable class, takes a list of positions and figures out which sites are to be occupied by the single type of adsorbate held. See getAdsorbateListForInpSites for callable interface

	"""
	def __init__(self, adsorbateObj, fractCoverage, addStrat=None, distance=0, fractTol=1e-2, nAds=None, useFractCover=True):
		""" Initializer
		
		Args:
			adsorbateObj: (Adsorbate obj) contains the geometry for the single adsorbate
			fractCoverage: (float) Fractional coverage for the surface; if inconsistent with the number of sites then it will throw an error later (upon call, not initialisation)
			addStrat: (BaseAdsorbateSiteOccupierStrategy, Optional) function that determines how we choose which sites to occupy. By default it will just occupy the first n-sites its passed (can be thought of as random, only really useful for full coverage)
			distance: (float, Optional) Distance from the adsorption site (NOT neccesarily a surface atom) to place the adsorbate objects. Some addStrat may use this info
			fractTol: (float,optional) Sites will be occupied such that actual fract coverage equals fractCoverage+-fractTol, if this is not possible an error will be thrown

			nAds: (int, optional) Specificy the integer number of adsorbates to use. This only takes effect if useFractCover=False
			useFractCover:(Bool, optional, default=True) Whether to use fractional coverage to determine number of adsorbates. If False will use nAds to determine

		"""
		self.adsorbateObj = adsorbateObj
		self.fractCoverage = fractCoverage
		self.addStrat = addStrat if addStrat is not None else addStratHelp.OccupyInOrderSiteOccupier()
		self.fractTol = fractTol
		self.distance = distance
		self.nAds = nAds
		self.useFractCover = useFractCover

	def getDistances(self, sitePositions, surfObj=None, surfaceVector=None):
		adsorbateSites = self.getAdsorbateListForInpSites(sitePositions, surfObj=surfObj, surfaceVector=surfaceVector)
		outDistances = list()
		for x in adsorbateSites:
			if x is None:
				outDistances.append(None)
			else:
				outDistances.append(self.distance)
		return outDistances

	def getAdsorbateListForInpSites(self, sitePositions, surfObj=None, surfaceVector=None):
		adsorbateObjs = self.getInputAdsorbateObjs( len(sitePositions) )
		distances = [self.distance for x in adsorbateObjs]
		outAdsorbates = self.addStrat(sitePositions, adsorbateObjs, surfaceVector=surfaceVector, surfObj=surfObj, distances=distances)
		return outAdsorbates

	def getInputAdsorbateObjs(self, nSites):
		numbAdsorbates = self.getNumberOfAdsorbateObjs(nSites)
		return [copy.deepcopy(self.adsorbateObj) for x in range(numbAdsorbates)]

	def getNumberOfAdsorbateObjs(self,nSites):
		if self.useFractCover is False:
			return self.nAds

		integerFractCoverages = [ x/nSites for x in range(nSites+1) ]
		diffs = [abs(x-self.fractCoverage) for x in integerFractCoverages]
		nearestIdx = min( [(idx,x) for idx,x in enumerate(diffs)], key=lambda x:x[1] )[0]
		outFractCoverage = integerFractCoverages[nearestIdx]
		if abs(outFractCoverage-self.fractCoverage)>self.fractTol:
			currArgs = [self.fractCoverage, nSites, outFractCoverage]
			raise ValueError("{} is an invalid fractional coverage for {} sites; the nearest sensible value is {}".format(*currArgs))

		return nearestIdx #+1 to account for the zero-based indexing


class AddWaterAdsorbatesToBilayerSitesStandard(BaseGetAdsorbatesForSites):

	def __init__(self, waterA, waterB, distA, distB=None):
		""" Initializer
		
		Args:
			waterA: (Adsorbate obj) Contains the geometry for the single adsorbate
			waterB: (Adsorbate obj) Contains the geometry for the single adsorbate
			distA: (float) Distance of the waterA adsorption from the adsorption
			distB: (float, Optional) If None then waterB is put at the same distance as waterA
				 
		"""
		self.waterA = waterA
		self.waterB = waterB
		self.distA = distA
		self.distB = distB if distB is not None else distA


	def getDistances(self, sitePositions, surfObj=None, surfaceVector=None):
		waterTypes = self._getWaterTypeForEachSite(sitePositions)
		outDists = list()
		for x in waterTypes:
			if x=="A":
				outDists.append( self.distA )
			else:
				outDists.append( self.distB )
		return outDists

	#Current implementation forms a graph of nearest neighbours, with distance set to "1" between
	# each node. We then get a distance matrix linking each node. If we assume the source point is
	# "waterA" then any a distance of 1 must be "waterB" and distance 2 must be "waterA" etc.
	def getAdsorbateListForInpSites(self, sitePositions, surfObj=None, surfaceVector=None):
		sites = self._getWaterTypeForEachSite(sitePositions)
		outSites = list()
		for x in sites:
			if x=="A":
				outSites.append( copy.deepcopy(self.waterA) )
			else:
				outSites.append( copy.deepcopy(self.waterB) )

		return outSites

	def _getWaterTypeForEachSite(self, sitePositions):
		adjacencyMatrix = self._getAdjacencyMatrix(sitePositions)
		spFormat = csr_matrix(adjacencyMatrix)
		distMatrix = dijkstra(spFormat) 
		outSites = list()
		for val in distMatrix[0]:
			if (round(val)%2 == 0):
				outSites.append( "A" )
			else:
				outSites.append( "B" )
		return outSites

	def _getAdjacencyMatrix(self, sitePositions):
		nearestNebDist = cartHelp.getClosestDistanceBetweenTwoPoints(sitePositions)
		nPos = len(sitePositions)
		outMatrix = np.zeros( (nPos,nPos), dtype=np.int32 )

		distTol = 1e-1
		for rowIdx in range(nPos):
			for colIdx in range(nPos):
				if rowIdx==colIdx:
					outMatrix[rowIdx][colIdx] = 0
				else:
					posA, posB = sitePositions[rowIdx], sitePositions[colIdx]
					currDist = vectHelp.getDistTwoVectors(posA,posB)
					if abs(currDist-nearestNebDist) < distTol:
						outMatrix[rowIdx][colIdx] = 1

		return outMatrix


