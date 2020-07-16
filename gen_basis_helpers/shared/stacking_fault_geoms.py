
import copy
import itertools as it
import math

import numpy as np

from . import cart_coord_utils as cartCoordHelp
from . import geom_constraints as geomConstraintHelp
from . import plane_equations as planeEqnHelp
from . import simple_vector_maths as vectHelp

class BaseStackingFaultGeomGenerator():
	""" Base class for mapping a perfect cell geometry (generally a surface) and a fractional displacement into a structure with a stacking fault 

	"""

	def getGeomForGivenDisplacement(self, inpGeom, displacement, centralIdx=None, planeTolerance=None):
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

	def getGeomConstraints(self, inpGeom):
		""" Get a geometry constraints object for the input geometry
		
		Args:
			inpGeom: (plato_pylib UnitCell object)
				 
		Returns
			geoConstraints: (GeomConstraints object) Contains a representation of the atoms that need to be fixed when running geometry opts
	 
		"""
		raise NotImplementedError("")


class HcpStackingFaultGeomGeneratorTemplate(BaseStackingFaultGeomGenerator):


	def getGeomForGivenDisplacement(self,inpGeom, displacement, centralIdx=None, planeTolerance=None):
		outCell = copy.deepcopy(inpGeom)
		self._applyPerfectToDispZeroDisplacementsToCell(outCell,planeTolerance)
		self._displaceCellInPlace(outCell, displacement, centralIdx, planeTolerance)
		return outCell


	def _applyPerfectToDispZeroDisplacementsToCell(self, inpCell, planeTolerance=5e-2):
		raise NotImplementedError("")

	def _getCentralAtomIdx(self, inpCell, planeTolerance=None):
		planeTolerance = self.planeTolerance if planeTolerance is None else planeTolerance

		if self.centralIdx is not None:
			return self.centralIdx

		surfacePlane = cartCoordHelp.getABPlaneEqnWithNormVectorSameDirAsC_uCellInterface(inpCell)

		#group atoms into their planes
		uniquePlaneDists, uniquePlaneAtomIndices = _getUniquePlaneDistsAndAtomIndicesFromSurfacePlaneAndCartCoords(surfacePlane, inpCell.cartCoords, planeTolerance)

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


def _getUniquePlaneDistsAndAtomIndicesFromSurfacePlaneAndCartCoords(inpPlaneEqn, cartCoords, planeTolerance):
	idxVsDist = _getAtomIndicesVsSignedDistsFromPlane(inpPlaneEqn, cartCoords)
	return _getUniquePlaneDistsAndAtomIndicesFromIdxVsDistList(idxVsDist, planeTolerance)

def _getAtomIndicesVsSignedDistsFromPlane(inpPlaneEqn, cartCoords):
	idxVsDist = list()
	#Get the distances of each atom from the bottom/top plane of the surface
	for idx,coord in enumerate(cartCoords):
		currDist = inpPlaneEqn.getSignedDistanceOfPointFromPlane(coord[:3])
		idxVsDist.append( (idx,currDist) )

	return idxVsDist

def _getUniquePlaneDistsAndAtomIndicesFromIdxVsDistList(idxVsDist, planeTolerance):
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

	return uniquePlaneDists, uniquePlaneAtomIndices


class HcpI1StackingFaultGeomGenerator(HcpStackingFaultGeomGeneratorTemplate):

	def __init__(self, centralIdx=None, planeTolerance=5e-2):
		self.centralIdx = centralIdx
		self.planeTolerance = planeTolerance

	#TODO: We want to make it so centralAtomIdx is basically unchanged w.r.t doing the next stacking fault.
	#Using ANY centralIdx in the "new" plane is probably good enough though
	def _applyPerfectToDispZeroDisplacementsToCell(self, inpCell, planeTolerance=5e-2):
		planeTolerance = self.planeTolerance if planeTolerance is None else planeTolerance

		#Get index/coorodinates of one atom in our "A" plane
		startCartCoords = inpCell.cartCoords
		outCartCoords = copy.deepcopy(inpCell.cartCoords)
		centralIdx = self._getCentralAtomIdx(inpCell, planeTolerance)
		centralCoords = inpCell.cartCoords[centralIdx][:3]
		surfacePlaneEqn = cartCoordHelp.getABPlaneEqnWithNormVectorSameDirAsC_uCellInterface(inpCell)
		centralPlaneDist = surfacePlaneEqn.getSignedDistanceOfPointFromPlane( centralCoords )

		#Step 1 - we need to group all the AB levels below centralIdx; start by getting all and then filtering
		uniquePlaneDists, uniquePlaneAtomIndices = _getUniquePlaneDistsAndAtomIndicesFromSurfacePlaneAndCartCoords(surfacePlaneEqn, inpCell.cartCoords, planeTolerance)
		filteredPlaneDists, filteredPlaneAtomIndices = list(), list()
		for pDist, atomIndices in it.zip_longest(uniquePlaneDists, uniquePlaneAtomIndices):
			if pDist-centralPlaneDist > planeTolerance: #Only want planes below the central atom plane
				filteredPlaneDists.append(pDist)
				filteredPlaneAtomIndices.append(atomIndices)
		assert len(filteredPlaneDists)%2==0, "Need to have an even number of planes to apply shifts to"

		#Step 2 - swap z co-ordinates of all the relevant A planes wih the B planes below [we've forced the plane to lie in xy due to use of the hcp0001 surface]
		planeIdxVsDist = [(idx,dist) for idx,dist in enumerate(filteredPlaneDists)]
		sortedIndicesVsDists = sorted(planeIdxVsDist, key=lambda x:x[1], reverse=True) #sorts in terms of high to low (higest plane is first index)

		for idx in range(0,len(sortedIndicesVsDists),2):
			aPlaneIdx, bPlaneIdx = sortedIndicesVsDists[idx][0], sortedIndicesVsDists[idx+1][0]
			aPlaneZ = startCartCoords[filteredPlaneAtomIndices[aPlaneIdx]][2]
			bPlaneZ = startCartCoords[filteredPlaneAtomIndices[bPlaneIdx]][2]
			#Sort out plane A->B
			for atomIdx in filteredPlaneAtomIndices[aPlaneIdx]:
				outCartCoords[atomIdx][2] = bPlaneZ

			#Sort out plane B->A
			for atomIdx in filteredPlaneAtomIndices[bPlaneIdx]:
				outCartCoords[atomIdx][2] = aPlaneZ

#		for idx, (pIdx,pDist) in enumerate(sortedIndicesVsDists):
#			sortedIndicesVsDists[idx]



	def _displaceCellInPlace(self, inpGeom, displacement, centralIdx=None, planeTolerance=None):
		pass

#TODO: Can likely factor a lot of this out into a standard class (Template pattern)
class HcpI2StackingFaultGeomGenerator(HcpStackingFaultGeomGeneratorTemplate):

	def __init__(self, centralIdx=None, planeTolerance=5e-2, fraction_10m10=1, fraction_m2110=0):
		self.centralIdx = centralIdx
		self.planeTolerance = planeTolerance
		self.fraction_10m10 = fraction_10m10
		self.fraction_m2110 = fraction_m2110

	#The zero displacement case is the perfect cell in this instance
	def _applyPerfectToDispZeroDisplacementsToCell(self, inpCell, planeTolerance=5e-2):
		pass


	def _displaceCellInPlace(self, inpGeom, displacement, centralIdx=None, planeTolerance=None):
		_checkAnglesConsistentWithHcp0001(inpGeom)
		displaceFactorOneVector = self._getDisplacementVectorForDispParamEqualsOne(inpGeom)
		displaceVector = [x*displacement for x in displaceFactorOneVector]
		self._applyDisplacementVectorToRelevantAtomsInCell(inpGeom, displaceVector, centralIdx, planeTolerance)


	def _getUnitDisplacementVector(self, inpGeom):
		aVect,bVect = inpGeom.lattVects[0], inpGeom.lattVects[1]
		assert aVect[-1]==0
		assert bVect[-1]==0
		a1Vect, a2Vect = vectHelp.getUnitVectorFromInpVector(aVect), vectHelp.getUnitVectorFromInpVector(bVect)

		#Find the third unit vector; this is 120 degrees from the both the first two (2 equations) and its z component is 0
		cos120 = math.cos(math.radians(120))
		twoDimVectA3 = np.linalg.inv( np.array( (aVect[:2],bVect[:2]) ) ) @ np.array( [cos120, cos120] )
		a3Vect = vectHelp.getUnitVectorFromInpVector( [x for x in twoDimVectA3] + [0] )

		#Get a displacement vector
		vect10m10 = [a-b for a,b in it.zip_longest(a1Vect,a3Vect)]
		vectm2110 = [(-2*a)+b+c for a,b,c in it.zip_longest(a1Vect,a2Vect,a3Vect)]
		outVect = vectHelp.getUnitVectorFromInpVector( [(a*self.fraction_10m10)+ (b*self.fraction_m2110) for a,b in it.zip_longest(vect10m10,vectm2110)] )

		return outVect

	def _getDisplacementVectorForDispParamEqualsOne(self, inpGeom):
		#Start by assuming the m2110 is ALWAYS true
		dispVectorMagnitude = self._getDispMagnitudeForFactorEqualsOne(inpGeom)
		dispUnitVector = self._getUnitDisplacementVector(inpGeom)
		outVector = [x*dispVectorMagnitude for x in dispUnitVector]
		return outVector


	def _getDispMagnitudeForFactorEqualsOne(self, inpGeom):
		atomIdx = self._getCentralAtomIdx(inpGeom)
		lattVects = inpGeom.lattVects
		surfacePlane = planeEqnHelp.ThreeDimPlaneEquation.fromTwoPositionVectors(lattVects[0],lattVects[1])
		nearestInPlaneNebDistance = cartCoordHelp.getNearestInPlaneDistanceGivenInpCellAndAtomIdx(inpGeom, atomIdx, surfacePlane) #Works if we assume all are the same; which they should be
		dispVectorMagnitude = (1/3)*nearestInPlaneNebDistance
		return dispVectorMagnitude

	#TODO: Factor most of this out into a function
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


	def getGeomConstraints(self, inpGeom):
		coords = inpGeom.cartCoords
		cellConstraints = geomConstraintHelp.CellConstraints([True,True,True],[True,True,True])
		cartPosConstraints = list()
		for idx,unused in enumerate(coords):
			currConstraint = geomConstraintHelp.AtomicCartesianConstraint(idx, fixX=True)
			cartPosConstraints.append(currConstraint)
		atomicPosConstraints = geomConstraintHelp.AtomicPositionConstraints( atomicCartConstraints=cartPosConstraints)
		#Code above is actually old, wrong code. fixY=True would make it technically more correct, but not allow full relaxation
		#Need to code up ability to constrain atoms in planes to get proper constraints for general stacking faults
		raise NotImplementedError("") 
		return geomConstraintHelp.GeomConstraints(atomicPosConstraints,cellConstraints)




def _checkAnglesConsistentWithHcp0001(inpCell, errorTol=5e-1):
	expAngles = [90,90,120]
	actAngles = inpCell.getLattAnglesList()
	angleDiffs = [abs(exp-act) for exp,act in it.zip_longest(expAngles,actAngles)]
	if any([x>errorTol for x in angleDiffs]):
		raise ValueError("Angles need to be {}, not {}".format(expAngles,actAngles))


