
import copy
import itertools as it
import math

import numpy as np

import plato_pylib.shared.ucell_class as uCellHelp
import plato_pylib.utils.supercell as supCellHelp

from . import surf_to_sites_shared as surfToSiteHelp
from ..shared import cart_coord_utils as cartHelp
from ..shared import plane_equations as planeEqn
from ..shared import simple_vector_maths as vectHelp




class Hcp0001SurfaceToSitesSharedMixin():
	
	def getOutwardsSurfaceVectorFromSurface(self, inpSurface):
		inpGeom = inpSurface.unitCell
		surfPlaneEqn = self.getSurfacePlaneEqn(inpGeom)
		return vectHelp.getUnitVectorFromInpVector(surfPlaneEqn.coeffs[:3])

	#TODO: I should likely be grouping close-planes together (using a tolerance criterion)
	def getSurfacePlaneEqn(self, inpCell):
		abPlaneEqn = cartHelp.getABPlaneEqnWithNormVectorSameDirAsC(inpCell.lattVects)
		return cartHelp.getPlaneEqnForOuterSurfaceAtoms(inpCell, top=self.top)

	def getFirstLayerAtomPositions(self, inpCell):
		surfPlane = self.getSurfacePlaneEqn(inpCell)
		cartCoords = [x[:3] for x in inpCell.cartCoords]
		indicesInSurfPlane = cartHelp.getFilteredIndicesForCoordsInInputPlane(cartCoords, surfPlane, planeTolerance=self.minInterPlaneDist)
		return [cartCoords[x] for x in indicesInSurfPlane]


class HcpSurfaceToHcpSites(Hcp0001SurfaceToSitesSharedMixin, surfToSiteHelp.BaseSurfaceToSites):

	def __init__(self, top=True):
		self.top = top
		self.minInterPlaneDist = 0.5

	def getSurfaceSitesFromInpSurface(self, inpSurface):
		inpGeom = inpSurface.unitCell
		surfPlaneEqn = self.getSurfacePlaneEqn(inpGeom)
		secondLayerPlaneEqn = self.getSecondLayerPlaneEquation(inpGeom)
		secondLayerAtomPositions = self.getSecondLayerAtomPositions(inpGeom)

		#Find the distance between the parralel planes
		interPlaneDist = getDistanceTwoParralelPlanes(surfPlaneEqn, secondLayerPlaneEqn)
		unitOutwardsSurfVector = vectHelp.getUnitVectorFromInpVector(surfPlaneEqn.coeffs[:3])
		outwardsSurfVector = [interPlaneDist*x for x in unitOutwardsSurfVector]

		#Get the adsorbate site positions
		outPositions = list()
		for aPos in secondLayerAtomPositions:
			atomPos = list(aPos)
			sitePos = [ x+y for x,y in it.zip_longest(atomPos,outwardsSurfVector)]
			outPositions.append(sitePos)

		return outPositions
	
	def getSecondLayerPlaneEquation(self,inpCell):
		surfPlane = self.getSurfacePlaneEqn(inpCell)
		cartCoords = [x[:3] for x in inpCell.cartCoords]
		nearestAtomContainingPlane = getNearestAtomContainingPlane(cartCoords, surfPlane, minDistance=self.minInterPlaneDist)
		return nearestAtomContainingPlane

	def getSecondLayerAtomPositions(self, inpCell):
		surfPlane = self.getSurfacePlaneEqn(inpCell)
		cartCoords = [x[:3] for x in inpCell.cartCoords]
		nearestAtomContainingPlane = getNearestAtomContainingPlane(cartCoords, surfPlane, minDistance=self.minInterPlaneDist)
		indicesInSecondPlane = cartHelp.getFilteredIndicesForCoordsInInputPlane(cartCoords, nearestAtomContainingPlane, planeTolerance=self.minInterPlaneDist)
		return [cartCoords[x] for x in indicesInSecondPlane]

	def __call__(self, inpSurface):
		return self.getSurfaceSitesFromInpSurface(inpSurface)


class HcpSurfaceToWaterBilayerSites(Hcp0001SurfaceToSitesSharedMixin, surfToSiteHelp.BaseSurfaceToSites):

	def __init__(self, top=True, firstUnoccStrat=None):
		""" Description of function
		
		Args:
			top: (Optional, Bool) If True
			firstUnoccStrat: (Optional, f(sitePositions)) Function to determine the first site to be unoccupied

		Returns
			What Function Returns
	 
		Raises:
			Errors
		"""
		self.top = top
		self.minInterPlaneDist = 0.5
		self.firstUnoccStrat = lambda x:x[0] if firstUnoccStrat is None else firstUnoccStrat


	def getSurfaceSitesFromInpSurface(self, inpSurface):
		atopGetter = HcpSurfaceToAtopSites(top=self.top)
		atopSites = atopGetter(inpSurface)
		surfVector = self.getOutwardsSurfaceVectorFromSurface(inpSurface)
		outIndices = self._getBilayerOccupiedSiteIndicesFromAtopSites(atopSites, surfVector)
		outSites = [x for idx,x in enumerate(atopSites) if idx in outIndices]
		return outSites

	def _getBilayerOccupiedSiteIndicesFromAtopSites(self, atopPositions, surfVector):
		firstPosition = self.firstUnoccStrat(atopPositions)
		vectA, vectB = self._getBilayerCellVectorsFromAtopPositions(atopPositions, surfVector)

		unoccIndices = list()
		for idx,pos in enumerate(atopPositions):
			if self._siteIsReachable(pos, firstPosition, vectA, vectB, surfVector):
				unoccIndices.append(idx)

		occIndices = [idx for idx,x in enumerate(atopPositions) if idx not in unoccIndices]

		fractCoverage = len(occIndices) / ( len(occIndices)+len(unoccIndices) )
		if abs(fractCoverage-0.66)>1e-2:
			raise ValueError("Bilayer sites/atop Sites = {}, this value should be 0.66.".format(fractCoverage))

		return occIndices

	def _siteIsReachable(self, testPos, startPos, vectA, vectB, surfVector, intTol=1e-1):
		diffVector = np.array([x-y for x,y in it.zip_longest(testPos, startPos)])
		vectorMatrix = np.array( [np.array(vectA), np.array(vectB), np.array(surfVector)] ).transpose()
		outCoeffs = np.dot(np.linalg.inv(vectorMatrix), diffVector)
		isReachable=False
		if all([abs(round(x)-x)<intTol for x in outCoeffs]):
			isReachable = True
		return isReachable

	def _getBilayerCellVectorsFromAtopPositions(self, atopPositions, surfVector):
		firstSite = self.firstUnoccStrat(atopPositions)
		idxFirstSite = self._getIdxOfSiteInList(firstSite,atopPositions)
		otherSitePositions = [x for idx,x in enumerate(atopPositions) if idx!=idxFirstSite]
		nearestSiteCoords = cartHelp.getNearestCoordToInputPoint(firstSite, otherSitePositions)

		#Get the three-vector basis
		rotatePlus120Matrix = vectHelp.getRotationMatrixAroundAxis(surfVector,120)
		vectA = [x-y for x,y in it.zip_longest(nearestSiteCoords,firstSite)]
		vectB = np.dot( rotatePlus120Matrix, vectA )
		vectC = np.dot( rotatePlus120Matrix, vectB )

		outVectA = [x-y for x,y in it.zip_longest(vectA, vectB)]
		outVectB = [x-y for x,y in it.zip_longest(vectA, vectC)]
		return outVectA, outVectB

	def _getIdxOfSiteInList(self, sitePos, allPositions):
		distsFromSite = list()
		for currPos in allPositions:
			currDist = vectHelp.getDistTwoVectors(currPos,sitePos)
			distsFromSite.append(currDist)

		indices = [idx for idx,dist in enumerate(distsFromSite) if dist <1e-4]
		assert len(indices) == 1
		return indices[0]

class HcpSurfaceToAtopSites(Hcp0001SurfaceToSitesSharedMixin, surfToSiteHelp.BaseSurfaceToSites):

	def __init__(self, top=True):
		self.top = top
		self.minInterPlaneDist = 0.5

	def getSurfaceSitesFromInpSurface(self, inpSurface):
		return self.getFirstLayerAtomPositions(inpSurface.unitCell)



class HcpSurfaceToFccHollowSites(Hcp0001SurfaceToSitesSharedMixin,surfToSiteHelp.BaseSurfaceToSites):

	def __init__(self, top=True, foldCoordsIntoCell=True):
		self.top = top
		self.minInterPlaneDist = 0.5
		self.foldCoordsIntoCell = foldCoordsIntoCell

	def getSurfaceSitesFromInpSurface(self, inpSurface):
		#Step 1 - Get the hcp sites (remember to pass top=self.top)
		#Step 2 - Get the triangle of top-layer atoms around those sites
		#Step 3 - You can THEN get the neighbouring triangle in a consistent way (maybe ALWAYS fix origin closest to actual origin + always move along x when possible)
		#Step 3.5 - Note if you do B + (C-A) then the new triangle we care about is B,C,D (hollow site is in the middle of that)
		#Step 4 - The middle of this new triangle is where you then put the new site
		hcpSiteGenerator = HcpSurfaceToHcpSites(top=self.top)
		hcpSites = hcpSiteGenerator(inpSurface)
		surfacePlaneEqn = self.getSurfacePlaneEqn(inpSurface.unitCell)
		geomWithImages = supCellHelp.getUnitCellSurroundedByNeighbourCells(inpSurface.unitCell, alongC=False)
		outFccSites = [self._getFccSiteFromHcpSite(geomWithImages, surfacePlaneEqn, x) for x in hcpSites]
		if self.foldCoordsIntoCell:
			outFccSites = self._getOutSitesFoldedIntoCell(outFccSites, inpSurface.unitCell)
		return outFccSites

	def _getFccSiteFromHcpSite(self, geomWithImages, surfacePlaneEqn, hcpXyz):
		#Firstly we get the three in-plane atom co-ordinates which are around the hcp site
		inpCoords = [x[:3] for x in geomWithImages.cartCoords]
		nearestThreeInPlaneIndices = cartHelp.getIndicesNearestInPlanePointsUpToN(3, hcpXyz, inpCoords, surfacePlaneEqn)

		#Now we get a neighbouring triangle, which contains the fcc-hollow site ~ at the centre
		coordA, coordB, coordC = _getCoordsSortedByDistanceFromOrigin([inpCoords[x] for x in nearestThreeInPlaneIndices])
		baVect = [x-y for x,y in it.zip_longest(coordB,coordA)]
		coordD = [x+y for x,y in it.zip_longest(coordC,baVect)]
		
		#Get the fcc hollow site position
		outSite = [(a+b+c)/3 for a,b,c in it.zip_longest(coordB, coordC, coordD)]

		return outSite

	def _getOutSitesFoldedIntoCell(self, outSites, inpCell):
		cartCoords = list(inpCell.cartCoords)
		newCoords = [x+["fake_ele"] for x in outSites]
		inpCell.cartCoords = cartCoords + newCoords
		uCellHelp.foldAtomicPositionsIntoCell(inpCell)
		outCartCoords = inpCell.cartCoords
		outSites = list()
		for x in outCartCoords:
			if x[3] == "fake_ele":
				outSites.append(x[:3])
		return outSites


class HcpSurfaceToBridgeSites(Hcp0001SurfaceToSitesSharedMixin, surfToSiteHelp.BaseSurfaceToSites):

	def __init__(self, top=True, alongA=True, alongB=True):
		""" Initialiser
		
		Args:
			top: (Optional, Bool) If True then sites are added to the top-layer (defined as having most +ve components along c)
			alongA: (Optional, Bool) If True include bridge sites along the first lattice vector
			alongB: (Optional, Bool) If True include bridge sites along the second lattice vector
 
		"""
		self.minInterPlaneDist = 0.5
		self.top = top
		self.alongA = alongA
		self.alongB = alongB

	def getSurfaceSitesFromInpSurface(self, inpSurface):
		startCell = inpSurface.unitCell
		uVectA, uVectB = [vectHelp.getUnitVectorFromInpVector(x) for x in startCell.lattVects[:2]]
		topSites = self.getFirstLayerAtomPositions(startCell)
		surfPlaneEqn = self.getSurfacePlaneEqn(startCell)
		cartCoords = [x[:3] for x in startCell.cartCoords]
		planeAtomIndices = cartHelp.getFilteredIndicesForCoordsInInputPlane(cartCoords, surfPlaneEqn, planeTolerance=self.minInterPlaneDist)
		
		#NOTE: This assumes ~perfect surface; namely that each top-layer atom has 3-equidistant neighbours
		aSites, bSites = list(), list()
		for idx,site in enumerate(topSites):
			currIdxInFullCell = planeAtomIndices[idx]
			nearestNebDist = cartHelp.getNearestInPlaneDistanceGivenInpCellAndAtomIdx(startCell, currIdxInFullCell, surfPlaneEqn)
			currASite = [x+(t*nearestNebDist*0.5) for x,t in it.zip_longest(site,uVectA)]
			currBSite = [x+(t*nearestNebDist*0.5) for x,t in it.zip_longest(site,uVectB)]
			aSites.append(currASite), bSites.append(currBSite)

		#Filter based on whether we're populating all sites or not
		outSites = list()
		if self.alongA:
			outSites.extend( aSites )
		if self.alongB:
			outSites.extend( bSites )

		return outSites


def _getCoordsSortedByDistanceFromOrigin(inpCoords):
	idxVsCoord = [x for x in enumerate(inpCoords)]
	distsFromOrigin = [vectHelp.getDistTwoVectors([0,0,0],x) for x in inpCoords]
	idxVsDist = [x for x in enumerate(distsFromOrigin)]
	return [inpCoords[idx] for idx,dist in sorted(idxVsDist, key=lambda x:[0])]



def getNearestAtomContainingPlane(inputCoords, planeEquation, minDistance=1e-1):
	allDists = list()
	for x in inputCoords:
		currDist = planeEquation.getDistanceOfPointFromPlane(x)
		allDists.append(currDist)

	labelledDists = [ [idx,x] for idx,x in enumerate(allDists) ]
	filteredDists = [x for x in labelledDists if abs(x[1])>minDistance]
	minDistIdx = min (filteredDists, key=lambda x:x[1]) [0]
	minDistPoint = inputCoords[minDistIdx]
	outDVal = planeEquation.calcDForInpXyz(minDistPoint)
	outCoeffs = list(planeEquation.coeffs[:3]) + [outDVal]
	return planeEqn.ThreeDimPlaneEquation(*outCoeffs)


def getDistanceTwoParralelPlanes(planeA,planeB):
	lenPlaneVectorA, lenPlaneVectorB = [vectHelp.getLenOneVector(x) for x in [planeA.coeffs[:3], planeB.coeffs[:3]]]
	assert abs(lenPlaneVectorA-lenPlaneVectorB)<1e-1 
	interPlaneDist = abs( planeA.coeffs[-1] - planeB.coeffs[-1] )/lenPlaneVectorA
	return interPlaneDist


def getDistFromFccHollowFromExpBondlengths(surfToAdsL, surfNearestNebL):
	""" Gets the expected distance between adsorbate and fcc hollow site based on expected bond lengths (should also work ok for hcp site)
	
	Args:
		surfToAdsL: (float) Distance between surface atom and adsorbate
		surfNearestNebL: (float) Distance between surface nearest neighbour atoms, lattice parameter a in the primitive cell

	Returns
		outDist: (float) Expected distance between fcc-site and the adsorbate
	
	Raises:
		ValueError: If surfToAdsL is too small a domain error will occur; this is because if distance is too short it becomes impossible for it to be in the centre of 3 surface atoms

	"""
	sideA = (0.5*surfNearestNebL) / math.cos(math.radians(30))
	return math.sqrt( (surfToAdsL**2) - (sideA**2) )



