
import itertools as it

import plato_pylib.shared.ucell_class as uCellHelp
import plato_pylib.utils.supercell as supCellHelp

from ..shared import cart_coord_utils as cartHelp
from ..shared import plane_equations as planeEqn
from ..shared import simple_vector_maths as vectHelp


class BaseSurfaceToSites():
	""" Callable class that takes a surface object and returns positions of adsorbate sites. Also has a function that returns a vector pointing out from the surface (which allows us to place adsorbates at varying distances from the adsorption sites

	"""

	def getOutwardsSurfaceVectorFromSurface(self, inpSurface):
		raise NotImplementedError("")

	def getSurfaceSitesFromInpSurface(self, inpSurface):
		raise NotImplementedError("")

	def __call__(self, inpSurface):
		return self.getSurfaceSitesFromInpSurface(inpSurface)


class Hcp0001SurfaceToSitesSharedMixin():
	
	def getOutwardsSurfaceVectorFromSurface(self, inpSurface):
		inpGeom = inpSurface.unitCell
		surfPlaneEqn = self.getSurfacePlaneEqn(inpGeom)
		return vectHelp.getUnitVectorFromInpVector(surfPlaneEqn.coeffs[:3])

	#TODO: I should likely be grouping close-planes together (using a tolerance criterion)
	def getSurfacePlaneEqn(self, inpCell):
		abPlaneEqn = cartHelp.getABPlaneEqnWithNormVectorSameDirAsC(inpCell.lattVects)

		#Find the plane-equations for the top or bottom plane [ASSUMES PERFECTLY FLAT SURFACE for now]
		#(Stolen from the generic surface code for now; TODO: Factor this out)
		allDVals = list()
		for x in inpCell.cartCoords:
			allDVals.append( abPlaneEqn.calcDForInpXyz(x[:3]) )
		maxD, minD = max(allDVals), min(allDVals)
		
		if self.top:
			outPlaneEquation = planeEqn.ThreeDimPlaneEquation( *(abPlaneEqn.coeffs[:3] + [maxD]) )
		else:
			outPlaneEquation = planeEqn.ThreeDimPlaneEquation( *(abPlaneEqn.coeffs[:3] + [minD]) )
			outPlaneEquation.coeffs = [x*-1 for x in outPlaneEquation.coeffs]
		return outPlaneEquation


class HcpTopSurfaceToSites(Hcp0001SurfaceToSitesSharedMixin, BaseSurfaceToSites):

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


class HcpSurfaceToFccHollowSites(Hcp0001SurfaceToSitesSharedMixin,BaseSurfaceToSites):

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
		hcpSiteGenerator = HcpTopSurfaceToSites(top=self.top)
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


