
import copy
import itertools as it
import math


import numpy as np

import plato_pylib.shared.ucell_class as uCell
import plato_pylib.utils.supercell as supCell

import gen_basis_helpers.shared.base_surface as baseSurface



class GenericSurface(baseSurface.BaseSurface):

	def __init__(self, singleLayerCell, nLayers, lenVac):
		""" Initialiser
		
		Args:
			singleLayerCell: A plato_pylib UnitCell object representing a single layer of the surface. 
			nLayers: The number of surface layers to use (1 layer = 1 singleLayerCell)
			lenvac: The amount of vacuum required between surface images
	
		"""

		self._singleLayerCell = singleLayerCell
		self._nLayers = nLayers
		self._lenVac = lenVac
	
	@property
	def unitCell(self):
		superCell = supCell.superCellFromUCell(self._singleLayerCell, [1,1,self._nLayers])
		addVacuumToUnitCellAlongC(superCell, self._lenVac)
		return superCell

	@property
	def surfaceArea(self):
		return (self._singleLayerCell.volume / self._singleLayerCell.lattParams["c"])


class Hcp0001Surface(baseSurface.BaseSurface):

	def __init__(self, bulkUCell, nLayers, lenVac):
		""" Initialiser for Hcp001Surface object
		
		Args:
			bulkUCell: A plato_pylib UnitCell object representing the bulk structure. This should be whats needed
			           to reperesent a single layer, so supercells with [n,m,1] dimensions should be used
			nLayers: The number of surface layers to use
			lenvac: The amount of vacuum required between surface images
		
		Raises:
			InvalidSurfaceError(ValueError): If the angles are wrong for the hcp surface 
		"""
		self._bulkCell = bulkUCell
		self._nLayers = nLayers
		self._lenVac = lenVac

		self._checkLattAnglesConsistent()

	def _checkLattAnglesConsistent(self):
		errorTol = 0.1
		expLattAngles = {"alpha":90.0, "beta":90.0, "gamma":120.0}
		expAlternativeLattAngles = {"alpha":90.0, "beta":90.0, "gamma":60.0}
		actLattAngles = self._bulkCell.lattAngles
		for key in expLattAngles.keys():
			if (abs(expLattAngles[key] - actLattAngles[key]) > errorTol) and (abs(expAlternativeLattAngles[key] - actLattAngles[key]) > errorTol):
				raise baseSurface.InvalidSurfaceError("Hcp surfaces need to be built from unit-cells ith angles of 90/90/120, but input has angles of {}".format(self._bulkCell.lattAngles))


	@property
	def unitCell(self):
		startCell = copy.deepcopy( self._bulkCell )
		outBulkCell = supCell.superCellFromUCell(self._bulkCell,[1,1,self._nLayers]) 
		addVacuumToUnitCellAlongC(outBulkCell, self._lenVac)
		return outBulkCell

	@property
	def surfaceArea(self):
		methA = self._calcSurfaceAreaUsingVolumeDivByC()
		methB = self._calcSurfaceAreaUsingTrig()
		allowedDiff = 0.01*max([methA,methB])
		assert abs(methB-methA) < allowedDiff, "Varying estimated for surface area, methA={}, methB={}".format(methA,methB)
		return self._calcSurfaceAreaUsingTrig()

	def _calcSurfaceAreaUsingVolumeDivByC(self):
		return (self._bulkCell.volume / self._bulkCell.lattParams["c"])

	def _calcSurfaceAreaUsingTrig(self):
		#latt params a/b define two sides of a triangle, with an angle of 120 between them. Hence cosine rule can be used to get
		#other info needed to get us a surface area
		lattParams = self._bulkCell.lattParams
		baseSqr = (lattParams["a"]**2) + (lattParams["b"]**2) - (2*lattParams["a"]*lattParams["b"]*math.cos(math.radians(120)))
		base = math.sqrt(baseSqr)
		#Use herons formula to get the area of a triangle from 3 sides
		p = (lattParams["a"] + lattParams["b"] + base) / 2
		areaSqr = p*(p-lattParams["a"])*(p-lattParams["b"])*(p-base)
		area = math.sqrt(areaSqr)
		return area*2 #The cell surface is made of two of these trigangles (its a parallelogram in general)






class Rocksalt001Surface(GenericSurface):

	def __init__(self, singleLayerCell, nLayers, lenVac):
		""" Initialiser
		
		Args:
			singleLayerCell: A plato_pylib UnitCell object representing a single layer of the surface. This can be built from the primitive cell by using getSingleLayerRocksalt001FromPrimitiveCell; extra space can be added in x and y directions by passing the relevant supercell
			nLayers: The number of surface layers to use. Note each layer is two atoms deep
			lenvac: The amount of vacuum required between surface images
		
		"""
		super().__init__(singleLayerCell, nLayers, lenVac)


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



