
import itertools as it

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
		return getSurfaceSitesFromInpSurface(inpSurface)


#TODO: Make this much neater even just for the "top-only" implementation before moving on
class HcpTopSurfaceToSites(BaseSurfaceToSites):

	def __init__(self, top=True):
		self.top = top
		self.minInterPlaneDist = 0.5

	def getOutwardsSurfaceVectorFromSurface(self, inpSurface):
		inpGeom = inpSurface.unitCell
		surfPlaneEqn = self.getSurfacePlaneEqn(inpGeom)
		return vectHelp.getUnitVectorFromInpVector(surfPlaneEqn.coeffs[:3])

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

	#TODO: This can likely be shared between all surface site types
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


def groupCoordsBasedOnSignedDistanceFromPlane(inputCoords, planeEquation, planeTol=1e-1):
	""" Takes a list of co-ordinates an returns indices grouped based on signed distances from an input plane equation
	
	Args:
		inputCoords: (iter of len-3 iters) Each contains [x,y,z] co-ordinates
		planeEquation: (ThreeDimPlaneEquation object) Represents a plane (usually surface plane used)
		planeTol: (float) The maximum planar separation between two co-ordinates for them to be considered in the same plane
 
	Returns
		indices: (iter of iter) Each contains the indices of the co-ordinates in a plane (the same signed distance from input planeEquation)
		avDistances: (iter of iter) Each index (idx) contains the average signed distance of indices[idx] from the input plane
 
	"""
	#TODO:Its VERY difficult to know how to group into planes, e.g. imagine tolerance of 0.1 with distances and -0.8,-0.9,-1.0; theres no UNIQUE grouping here
	#TODO: Test the grouping function separately with various test cases
	#STEP 1: get a set of distances
	#STEP 2: 
	pass



