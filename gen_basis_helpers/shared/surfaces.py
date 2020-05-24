
import copy
import itertools as it
import math


import numpy as np

import plato_pylib.shared.ucell_class as uCell
import plato_pylib.utils.supercell as supCell

import gen_basis_helpers.shared.base_surface as baseSurface



class GenericSurface(baseSurface.BaseSurface):

	def __init__(self, singleLayerCell, nLayers, lenVac=None, lenAbsoluteVacuum=None):
		""" Initialiser
		
		Args:
			singleLayerCell: A plato_pylib UnitCell object representing a single layer of the surface. The surface needs to be defined along the ab axes
			nLayers: The number of surface layers to use (1 layer = 1 singleLayerCell)
			EXACTLY ONE of these following parameters needs to be set to something other than None
			lenvac: (float, optional-ish) The amount of vacuum TO ADD between surface images
			lenAbsoluteVacuum: (float, optional-ish) The amount of vacuum between surfaces along the surface normal vector (i.e. the minimum distance betweeen surface planes when applying periodic boundary conditions)
	
		"""

		self._singleLayerCell = singleLayerCell
		self._nLayers = nLayers
		self._setLenVacFromInitializerInputArgs(lenVac, lenAbsoluteVacuum)

	@property
	def lenVac(self):
		return self._lenVac

	def _setLenVacFromInitializerInputArgs(self, lenVac, lenAbsoluteVacuum):
		if (lenVac is None) and (lenAbsoluteVacuum is None):
			raise AttributeError("lenVac or lenAbsoluteVacuum need setting")

		if (lenVac is not None) and (lenAbsoluteVacuum is not None):
			raise AttributeError("Only ONE of lenVac and lenAbsoluteVacuum can be set; not both")

		if lenVac is not None:
			self._lenVac = lenVac
		elif lenAbsoluteVacuum is not None:
			self._lenVac = 0 #Needs to be defined as SOMETHING for the setter to work properly
			self.lenAbsoluteVacuum = lenAbsoluteVacuum 
		else:
			raise ValueError("Error: Code hit a section it never should")
	
	@property
	def unitCell(self):
		superCell = supCell.superCellFromUCell(self._singleLayerCell, [1,1,self._nLayers])
		vacToAdd = self._getAmountOfVacuumToAddAlongC()
		addVacuumToUnitCellAlongC(superCell, vacToAdd)
		return superCell

	@property
	def lenAbsoluteVacuum(self):
		""" Length of vacuum between surface planes in two periodic images (bottom plane in one cell and top plane in another). Put another way, this is the minimum interaction distance between atoms in different cells along c.
		"""
		lattVects = self._singleLayerCell.lattVects
		abPlaneEquation = ThreeDimPlaneEquation.fromTwoPositionVectors(lattVects[0],lattVects[1],normaliseCoeffs=True)

		#Find the plane-equations for the top plane
		allDVals = list()
		for x in self._singleLayerCell.cartCoords:
			allDVals.append( abPlaneEquation.calcDForInpXyz(x[:3]) )
		maxD = max(allDVals)
		topPlaneEquation = ThreeDimPlaneEquation( *(abPlaneEquation.coeffs[:3] + [maxD]) )

		#Find a position lying in the bottom plane for our single cell
		unused, bottomAtomIndex = min( (idx,val) for val,idx in enumerate(allDVals) ) #Any atom in this plane is fine
		bottomAtomPosition = self._singleLayerCell.cartCoords[bottomAtomIndex][:3]
		bottomPositionInCellAbove = [x1+x2 for x1,x2 in it.zip_longest(bottomAtomPosition,lattVects[2])]

		#Get the vacuum sepataion when we dont add any ADDITIONAL vacuum on
		outSepNoAdditionalVacuum = topPlaneEquation.getDistanceOfPointFromPlane(bottomPositionInCellAbove)

		#Get the vacuum separation
		outVacSep = outSepNoAdditionalVacuum + self._lenVac 

		return outVacSep

	@lenAbsoluteVacuum.setter
	def lenAbsoluteVacuum(self,val):
		currVal = self.lenAbsoluteVacuum
		diff = val-currVal
		self._lenVac += diff

	#When alpha or beta do not equal 90 degrees, we need to add a larger amount of vacuum to get the same
	#separation between surface planes
	def _getAmountOfVacuumToAddAlongC(self):
		lattVects = self._singleLayerCell.lattVects
		surfPlaneEquation = ThreeDimPlaneEquation.fromTwoPositionVectors(lattVects[0],lattVects[1],normaliseCoeffs=True)
		unitSurfaceNormal = surfPlaneEquation.coeffs[:3]
		unitCVector = [x/self._singleLayerCell.lattParams["c"] for x in lattVects[2]]
		componentAlongC = np.dot( np.array(unitSurfaceNormal), np.array(unitCVector) )
		vacToAdd = self._lenVac/componentAlongC
		return vacToAdd	

	@property
	def surfaceArea(self):
		gamma = self._singleLayerCell.lattAngles["gamma"]
		a,b,c = self._singleLayerCell.getLattParamsList()
		area = a*b*math.sin(math.radians(gamma)) #Area of a parralelogram: A=bh; h=a sin(theta)
		return area

class Hcp0001Surface(GenericSurface):


	def __init__(self, singleLayerCell, nLayers, lenVac=None, lenAbsoluteVacuum=None):
		""" Initialiser
		
		Args:
			singleLayerCell:  A plato_pylib UnitCell object representing the bulk structure. This should be whats needed
			           to reperesent a single layer, so supercells with [n,m,1] dimensions should be used
			nLayers: The number of surface layers to use (1 layer = 1 singleLayerCell)
			EXACTLY ONE of these following parameters needs to be set to something other than None
			lenvac: (float, optional-ish) The amount of vacuum TO ADD between surface images
			lenAbsoluteVacuum: (float, optional-ish) The amount of vacuum between surfaces along the surface normal vector (i.e. the minimum distance betweeen surface planes when applying periodic boundary conditions)
	
		"""

		self._singleLayerCell = singleLayerCell
		self._nLayers = nLayers
		self._setLenVacFromInitializerInputArgs(lenVac, lenAbsoluteVacuum)
		self._checkLattAnglesConsistent()

	def _checkLattAnglesConsistent(self):
		errorTol = 0.1
		expLattAngles = {"alpha":90.0, "beta":90.0, "gamma":120.0}
		expAlternativeLattAngles = {"alpha":90.0, "beta":90.0, "gamma":60.0}
		actLattAngles = self._singleLayerCell.lattAngles
		for key in expLattAngles.keys():
			if (abs(expLattAngles[key] - actLattAngles[key]) > errorTol) and (abs(expAlternativeLattAngles[key] - actLattAngles[key]) > errorTol):
				raise baseSurface.InvalidSurfaceError("Hcp surfaces need to be built from unit-cells ith angles of 90/90/120, but input has angles of {}".format(self._singleLayerCell.lattAngles))






class Rocksalt001Surface(GenericSurface):

	def __init__(self, singleLayerCell, nLayers, lenVac=None, lenAbsoluteVacuum=None):
		""" Initialiser
		
		Args:
			singleLayerCell: A plato_pylib UnitCell object representing a single layer of the surface. This can be built from the primitive cell by using getSingleLayerRocksalt001FromPrimitiveCell; extra space can be added in x and y directions by passing the relevant supercell
			nLayers: The number of surface layers to use. Note each layer is two atoms deep
			lenvac: (float, optional-ish) The amount of vacuum TO ADD between surface images
			lenAbsoluteVacuum: (float, optional-ish) The amount of vacuum between surfaces along the surface normal vector (i.e. the minimum distance betweeen surface planes when applying periodic boundary conditions)

		
		"""
		super().__init__(singleLayerCell, nLayers, lenVac=lenVac, lenAbsoluteVacuum=lenAbsoluteVacuum)


def addVacuumToUnitCellAlongC(inpCell, lenVac):
	#Need to add the vacuum, half at top and half at bottom
	lattParams = copy.deepcopy( inpCell.lattParams )
	lattParams["c"] += lenVac
	startLattVects = copy.deepcopy( inpCell.lattVects )

	cartCoords = copy.deepcopy(inpCell.cartCoords)

	for row in cartCoords:
		row[2] += 0.5*lenVac

	inpCell.lattParams = lattParams #NOTE: We MUST change lattParams first, since this changes the cart-coords of atoms
	inpCell.cartCoords = cartCoords


def getSingleLayerHcp1010FromPrimitiveCell(primCell):
	_checkInpCellIsHcpPrimitive(primCell)

	#Trying diff approach
	uVectA = [0,0,1] #this is for the 0001 surface
	uVectB = [1,0,0] #This should be for the 10m10 surface
	lattVectTransformMatrix = _getRotationMatrixLinkingTwoUnitVectors(uVectA, uVectB)

	#Cell vectors part
	outCell = copy.deepcopy(primCell)
	outCell.putCAlongZ = True
	finalLattVects = lattVectTransformMatrix.dot( np.array(outCell.lattVects) )
	outCell.lattVects = finalLattVects

	#Fract coords part
	fractCoords = _getResultFromTransformMatrixToFractCoordsInclAtomicSymbols(outCell.fractCoords, lattVectTransformMatrix)
	outCell.fractCoords = fractCoords
	_centreCFractCoordsForInpCell(outCell)

	return outCell

#https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
def _getRotationMatrixLinkingTwoUnitVectors(uVectA, uVectB):
	assert abs(1 - _getLengthOfVector(uVectA))<1e-3
	assert abs(1 - _getLengthOfVector(uVectB))<1e-3
	assert len(uVectA)==3
	assert len(uVectB)==3

	crossProd = np.cross(uVectA, uVectB)
	v1,v2,v3 = crossProd
	vx = np.array( [ [0    , -1*v3,  1*v2],
	                 [v3   ,  0   , -1*v1],
	                 [-1*v2, v1   , 0    ] ] )

	vxSquared = vx.dot(vx)
	cosTheta = np.dot(uVectA,uVectB)
	angularFactor = 1 / (1 + cosTheta)
	vxSquaredWithAngular = angularFactor*vxSquared
	rotMatrix = np.identity(3) + vx + vxSquaredWithAngular

#	np.identity(3)
	return rotMatrix

def _getLengthOfVector(inpVector):
	sqrSum = sum( [x**2 for x in inpVector] )
	return math.sqrt(sqrSum)


def _getResultFromTransformMatrixToFractCoordsInclAtomicSymbols(inpFractCoords, transformMatrix):
	coordPart = np.array( [x[:3] for x in inpFractCoords] )
	atomListPart = [x[-1] for x in inpFractCoords]

	outList = list()
	for row,atomSymb in zip(coordPart,atomListPart):
		currCoords = (transformMatrix.dot(row)).tolist()
		currCoords += [atomSymb]
		outList.append(currCoords)
	
	return outList


def _centreCFractCoordsForInpCell(inpCell):
	fractCoords =  inpCell.fractCoords
	zComps = [x[2] for x in fractCoords]
	highestVal, lowestVal = max(zComps), min(zComps)
	distFromTop, distFromBottom = 1-highestVal, lowestVal
	shiftValue = 0.5*(distFromTop - distFromBottom)
	outZComps = [x+shiftValue for x in zComps]
	for x in fractCoords:
		x[2] += shiftValue
	inpCell.fractCoords = fractCoords


def _checkInpCellIsHcpPrimitive(inpCell):
	a,b,c = inpCell.getLattParamsList()
	alpha, beta, gamma = inpCell.getLattAnglesList()
	assert abs(a-b) < 1e-4
	assert abs(c-b) > 1e-4
	assert len(inpCell.fractCoords) == 2
	assert abs(alpha-90) < 1e-1
	assert abs(beta-90) < 1e-1 
	assert abs(gamma-120) < 1e-1

def getSingleLayerRocksalt001FromPrimitiveCell(primCell):
	""" Description of function
	
	Args:
		primCell: (plato_pylib UnitCell object) This must be a primitive rock salt cell. Meaning 2 atoms, with each lattice angle=60 degrees and a=b=c
			 
	Returns
		 outCell: (plato_pylib UnitCell object) Contains a 4-atom unit-cell which forms the basis for forming rock-salt 001 surfaces (it is one layer)
 
	Raises:
		 AssertionError: If the input primitive cell is not a rock salt cell. Probably not garanteed to be raised if an incompatible cell is passed, but at least catches some errors.
	"""

	assert _uCellIsRockSaltPrimitive(primCell, printError=True), "UnitCell {} is not a primitive rock salt cell".format(primCell)

	
	#From here onwards i assume primCell is in the correct format
	inpLattParam = primCell.lattParams["a"]
	lattParam = primCell.lattParams["a"] * (1/math.sqrt(0.5))
	_tempA = 1/math.sqrt(2) #if a unit vector is halfway between x/y axes, then each component hsa this value
	lattVects = [ [_tempA*inpLattParam,     _tempA*inpLattParam, 0],
	              [-1*_tempA*inpLattParam, _tempA*inpLattParam, 0],
	              [0.0                   , 0.0,                lattParam] ]
	eleA, eleB = [x[-1] for x in primCell.fractCoords]
	fractCoords = [ [0.0,0.0,0.0,eleA], [0.0,0.0,0.5,eleB], [0.5,0.5,0.0,eleB], [0.5,0.5,0.5,eleA] ]



	return uCell.UnitCell.fromLattVects(lattVects, fractCoords)


def _uCellIsRockSaltPrimitive(inpCell, printError=True, angleTol=1e-1, lattParamTol=1e-1):
	isPrim = True
	printMsg = ""

	#Check right number of atoms
	if len(inpCell.fractCoords) != 2:
		printMsg = "2 atoms expected for rocksalt primitive cell, found {}".format( len(inpCell.fractCoords) )
		isPrim = False

	#Check correct angle
	angles = [x for x in inpCell.lattAngles.values()]
	expAngle = 60.0
	if not all([ abs((x-expAngle))<angleTol for x in angles]):
		printMsg = "Angles should also be 60.0 degrees for rocksalt primitive cell, found angles of {}".format(angles)
		isPrim=False

	#Check correct lattice paramter relationships
	lattParams = [x for x in inpCell.lattParams.values()]
	if not all([ abs(x-lattParams[0]) < lattParamTol for x in lattParams]):
		printMsg = "Lattice parameters should all be equal for rocksalt primitive cell, found parameters of {}".format(lattParams)
		isPrim=False

	if printError:
		print(printMsg)

	return isPrim



#TODO: Move this to its own module
class ThreeDimPlaneEquation():
	""" Class representing a 3-dimension plane equation ax + by +cz = d

	"""
	def __init__(self, a, b, c, d):
		""" Initializer
		
		Args:
			a,b,c,d (floats): Parameters for the plane equation ax + by + cz = d
				 
		"""
		self.a = a
		self.b = b
		self.c = c
		self.d = d

	@classmethod
	def fromTwoPositionVectors(cls, inpVectA, inpVectB, normaliseCoeffs=True):
		""" Alternative initializer. Creates the object using two input POSITION vectors (i.e. both MUST pass through origin)
		
		Args:
			inpVectA: (len 3 float iter) Position vector
			inpVectB: (len 3 float iter) Position vector
			normaliseCoeffs: (Bool, Optional) If True then always return coefficients for the normalised (i.e. unit) vector normal to the plane. Default=True
				 
		"""
		vectA, vectB = np.array(inpVectA), np.array(inpVectB)
		normVect = np.cross(vectA,vectB)
		if normaliseCoeffs:
			lenNormVect = math.sqrt( sum([x**2 for x in normVect]) )
			normVect = np.array( [x/lenNormVect for x in normVect] )

		#D is always zero, since the position vectors are also points on the plane
		outCoeffs = [x for x in normVect] + [0] 

		return cls(*outCoeffs)


	def getDistanceOfPointFromPlane(self, inpXyz):
		""" Calculates the distance of a point from the plane
		
		Args:
			inpXyz: (len 3 float iter) [x,y,z] co-ordinates
				 
		Returns
			outDist: (float) The distance between input point and the nearest point on this plane
	 
		"""
		#Step 1 = Find the d value for the parralel plane this lies on; the displacement vector between point and plane is then the vector normal to this plane
		#Step 2  = use the fact the normal vector points the same way as the displacement vecctor to get a length
		dValForThisPoint = self.calcDForInpXyz(inpXyz)
		diffInDVals = abs( dValForThisPoint - self.d )
		lenNormalVectorToThisPlane = math.sqrt( (self.a**2) + (self.b**2) + (self.c**2) )
		outDist = diffInDVals / lenNormalVectorToThisPlane
		return outDist


	def calcDForInpXyz(self, inpXyz):
		""" For a given xyz calculate d. If d=self.d then the point lies on this plane; else it lies on a parralel plane with d being the output to this function
		
		Args:
			inpXyz: (len 3 float iter) co-ordinates for a point
				 
		"""
		assert len(inpXyz) == 3
		return sum( [param*coord for param,coord in it.zip_longest([self.a,self.b,self.c],inpXyz)] )


	@property
	def coeffs(self):
		""" Returns a,b,c,d coeffs (in ax + by +cz = d) as a len-4 iter
		"""
		return [self.a,self.b,self.c,self.d]

