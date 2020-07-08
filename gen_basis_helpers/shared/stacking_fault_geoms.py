
import copy
import itertools as it
from . import cart_coord_utils as cartCoordHelp
from . import plane_equations as planeEqnHelp
from . import simple_vector_maths as vectHelp

class BaseStackingFaultGeomGenerator():
	""" Base class for mapping a perfect cell geometry (generally a surface) and a fractional displacement into a structure with a stacking fault 

	"""

	def getGeomForGivenDisplacement(inpGeom, displacement, centralIdx=None, planeTolerance=None):
		""" Return a geometry for the stacking fault displacement for a given input geometry
		
		Args:
			inpGeom: plato_pylib UnitCell object, contains the geometry for the perfect structure [this likely really needs to be perfect aswell, meaning a supercell based directly on a 2-atom cell]
			displacement: float, Value (usually 0 to 1) describing the magnitude of the stacking fault
			centralIdx: int (Optional), The index of an atom in the plane containing the stacking fault
			planeTolerance: float (Optional), The maximum separation between two atoms z-planes for them to be considered in the same plane
			
			NOTE: The optional attributes should default to values stored within the class
 
		Returns
			outGeom: plato_pylib UnitCell object containing the input geometry with the required stacking fault
	 
		"""
		raise NotImplementedError("")



#TODO: Can likely factor a lot of this out into a standard class (Template pattern)
class HcpI2StackingFaultGeomGenerator(BaseStackingFaultGeomGenerator):

	def __init__(self, centralIdx=None, planeTolerance=5e-2):
		self.centralIdx = centralIdx
		self.planeTolerance = planeTolerance

	def _getCentralAtomIdx(self, inpCell, planeTolerance=None):
		planeTolerance = self.planeTolerance if planeTolerance is None else planeTolerance

		if self.centralIdx is not None:
			return self.centralIdx

		surfacePlane = cartCoordHelp.getABPlaneEqnWithNormVectorSameDirAsC_uCellInterface(inpCell)
		idxVsDist = list()

		#Get the distances of each atom from the bottom/top plane of the surface
		for idx,coord in enumerate(inpCell.cartCoords):
			currDist = surfacePlane.getSignedDistanceOfPointFromPlane(coord[:3])
			idxVsDist.append( (idx,currDist) )

		#group atoms into their planes
		allDists = [ idxVsDist[0][1] ]
		uniquePlaneDists = [ idxVsDist[0][1] ]
		uniquePlaneAtomIndices = [ [idxVsDist[0][0]] ]
		for idx,dist in idxVsDist[1:]:
			sortedIntoPlane = False
			#Check if this plane has already been found
			for pIdx, pDist in enumerate(uniquePlaneDists):
				if abs(pDist-dist) < planeTolerance:
					uniquePlaneAtomIndices[pIdx].append(idx)
					sortedIntoPlane = True

			if not sortedIntoPlane:
				uniquePlaneDists.append( dist )
				uniquePlaneAtomIndices.append( [idx] ) 

		#Figure out the number of atoms in each plane; throw if not all equal
		nAtomsPerPlane = len(uniquePlaneAtomIndices[0]) 
		maxInPlane = max([len(x) for x in uniquePlaneAtomIndices])
		minInPlane = min([len(x) for x in uniquePlaneAtomIndices])
		assert all([len(x)==nAtomsPerPlane for x in uniquePlaneAtomIndices]), "Each ab plane should have the same number of atoms in it; but between {} and {} found in current system".format(minInPlane,maxInPlane)

		#Find the central plane index
		#Number of planes will always be even, and since we displace above we need to take the bottom of the 2 central as our pivot
		# (else we go B-A-B -> B-A-C instead of A-B-A -> B-C-A)
		planeIdxVsDist = [(idx,dist) for idx,dist in enumerate(uniquePlaneDists)]
		sortedIndicesVsDists = sorted(planeIdxVsDist, key=lambda x:x[1], reverse=True) #sorts in terms of high to low (higest plane is first index)
		nElements = len(sortedIndicesVsDists)
		assert nElements%2==0, "hcp 0001 Cells should have an even number of atomic layers"
		if int(nElements/2)%2 == 0:
			idxNeeded = int(nElements/2) #Even number of AB layers; if 4 we'd take the 5th plane from the bottom (3rd A in AB-AB-AB-AB listing high to low)
		else:
			idxNeeded = int(nElements/2)-1 #Odd number of AB layers; +1 or -1 would likely both be fine here

		planeIdx = sortedIndicesVsDists[ idxNeeded ][0]
		centralIdx = uniquePlaneAtomIndices[planeIdx][0]

		return centralIdx

	def getGeomForGivenDisplacement(self,inpGeom, displacement, centralIdx=None, planeTolerance=None):
		outCell = copy.deepcopy(inpGeom)
		self._displaceCellInPlace(outCell, displacement, centralIdx, planeTolerance)
		return outCell

	def _getDisplacementVectorForDispParamEqualsOne(self, inpGeom):
		atomIdx = self._getCentralAtomIdx(inpGeom)
		lattVects = inpGeom.lattVects
		surfacePlane = planeEqnHelp.ThreeDimPlaneEquation.fromTwoPositionVectors(lattVects[0],lattVects[1])
		nearestInPlaneNebDistance = cartCoordHelp.getNearestInPlaneDistanceGivenInpCellAndAtomIdx(inpGeom, atomIdx, surfacePlane) #Works if we assume all are the same; which they should be

		#I think ANY vector along the ab plane would work, hence we use a for simplicity (could make the class configurable if this becomes an issue)
		dispUnitVector = vectHelp.getUnitVectorFromInpVector(lattVects[0])
		dispVectorMagnitude = (1/3)*nearestInPlaneNebDistance

		return [x*dispVectorMagnitude for x in dispUnitVector]

	def _applyDisplacementVectorToRelevantAtomsInCell(self, inpGeom, displaceVector, centralIdx=None, planeTolerance=None):
		#Sort out default args
		centralIdx = self._getCentralAtomIdx(inpGeom) if centralIdx is None else centralIdx
		planeTolerance = self.planeTolerance if planeTolerance is None else planeTolerance

		lattVects = inpGeom.lattVects
		surfacePlane = cartCoordHelp.getABPlaneEqnWithNormVectorSameDirAsC_uCellInterface(inpGeom)

		relAtomDVal = surfacePlane.calcDForInpXyz( inpGeom.cartCoords[centralIdx][:3] )
		surfacePlane.d = relAtomDVal #Now specifically the surface-plane we're applying the shift to

		#Figure out which atoms are below the surface plane
		indicesOfAtomsBelowOrOnPlane = list()
		for idx,coords in enumerate(inpGeom.cartCoords):
			currSignedDist = surfacePlane.getSignedDistanceOfPointFromPlane(coords[:3])
			if (currSignedDist - planeTolerance <= 0):
				indicesOfAtomsBelowOrOnPlane.append(idx)

		#Apply the displacement
		outCoords = copy.deepcopy(inpGeom.cartCoords)
		for idx, coords in enumerate(outCoords):
			if idx in indicesOfAtomsBelowOrOnPlane:
				currCoords = [x+d for x,d in it.zip_longest(coords[:3],displaceVector)] + [coords[-1]]
				outCoords[idx] = currCoords

		inpGeom.cartCoords = outCoords

	def _displaceCellInPlace(self, inpGeom, displacement, centralIdx=None, planeTolerance=None):
		_checkAnglesConsistentWithHcp0001(inpGeom)
		displaceFactorOneVector = self._getDisplacementVectorForDispParamEqualsOne(inpGeom)
		displaceVector = [x*displacement for x in displaceFactorOneVector]
		self._applyDisplacementVectorToRelevantAtomsInCell(inpGeom, displaceVector, centralIdx, planeTolerance)



def _checkAnglesConsistentWithHcp0001(inpCell, errorTol=5e-1):
	expAngles = [90,90,120]
	actAngles = inpCell.getLattAnglesList()
	angleDiffs = [abs(exp-act) for exp,act in it.zip_longest(expAngles,actAngles)]
	if any([x>errorTol for x in angleDiffs]):
		raise ValueError("Angles need to be {}, not {}".format(expAngles,actAngles))


