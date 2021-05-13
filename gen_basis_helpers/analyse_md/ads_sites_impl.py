
import itertools as it
from . import ads_sites_core as coreAdsSiteHelp
from . import calc_dists as distsHelp

class TopStandard(coreAdsSiteHelp.FixedIndicesAdsSiteBase):

	def __init__(self, atomIdx, siteName="top"):
		""" Initializer. Represents top site for a surface (where adsorbate sites directly above a top-plane atom
		
		Args:
			atomIdx: (int) The index of the atom defining the top site
				 
		"""
		self.atomIdx = atomIdx
		self.siteName = siteName

	#May need this vectorised later? Which would likely mean a composite class
	def positionFromGeom(self, inpGeom, inpCartCoords=None, inpIdx=None):
		""" Returns position of this site when given an input geometry
		
		Args:
			inpGeom: (plato_pylib UnitCell object)
			inpCartCoords: (inpGeom.cartCoords) Passing these directly MIGHT be faster than getting them from inpGeom
			inpIdx: (int, Optional) The index of the atom defining the top site; this is stored on self but passing as an input argument will overwrite that
 
		Returns
			outCoord: (len-3 iter) [x,y,z] position of the site
	 
		"""
		cartCoords = inpGeom.cartCoords if inpCartCoords is None else inpCartCoords
		atomIdx = self.atomIdx if inpIdx is None else inpIdx
		return [x for x in cartCoords[atomIdx][:3]]


class BridgeStandard(coreAdsSiteHelp.FixedIndicesAdsSiteBase):

	def __init__(self, atomIndices, siteName="bridge"):
		""" Initializer. Represents bridge site; where adsorbate is equidistant between 2 surface atoms
		
		Args:
			atomIndices: (len-2 iter of ints) The indices of the two adsorbate sites
				 
		"""
		self.atomIndices = atomIndices
		self.siteName = siteName

	def positionFromGeom(self, inpGeom, inpCartCoords=None):
		cartCoords = inpGeom.cartCoords if inpCartCoords is None else inpCartCoords
		assert len(self.atomIndices)==2
		coordA = [x for x in cartCoords[self.atomIndices[0]][:3]]
		coordB = [x for x in cartCoords[self.atomIndices[1]][:3]]
		nearestImageCoordB = distsHelp.getNearestImageNebCoordsBasic(inpGeom, coordA, coordB)
		outCoord = [(a+b)/2 for a,b in it.zip_longest(coordA, nearestImageCoordB)]
		return outCoord


class HollowStandard(coreAdsSiteHelp.FixedIndicesAdsSiteBase):

	def __init__(self, atomIndices, siteName="hollow"):
		""" Initializer. Represents a hollow site, which will be at the centroid of 3 surface atoms
		
		Args:
			atomIndices: (len-3 iter of ints) The indices of the two adsorbate sites
				 
		"""
		self.atomIndices = atomIndices
		self.siteName = siteName

	def positionFromGeom(self, inpGeom, inpCartCoords=None):
		cartCoords = inpGeom.cartCoords if inpCartCoords is None else inpCartCoords
		coordsToAverage = list()
		firstCoord = cartCoords[ self.atomIndices[0] ][:3]
		secondCoord = distsHelp.getNearestImageNebCoordsBasic(inpGeom, firstCoord, cartCoords[self.atomIndices[1]][:3])
		thirdCoord  = distsHelp.getNearestImageNebCoordsBasic(inpGeom, firstCoord, cartCoords[self.atomIndices[2]][:3])

		coordsToAverage = [firstCoord, secondCoord, thirdCoord]

		nCoords = len(coordsToAverage)
		outX = sum([ x[0] for x in coordsToAverage ])/ nCoords
		outY = sum([ x[1] for x in coordsToAverage ])/ nCoords
		outZ = sum([ x[2] for x in coordsToAverage ])/ nCoords

		return [outX, outY, outZ]





