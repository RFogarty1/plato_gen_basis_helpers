
import copy
import itertools as it
import math


import numpy as np

import plato_pylib.shared.ucell_class as uCell
import plato_pylib.utils.supercell as supCell

from . import plane_equations as planeEqn

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
		abPlaneEquation = planeEqn.ThreeDimPlaneEquation.fromTwoPositionVectors(lattVects[0],lattVects[1],normaliseCoeffs=True)

		#Find the plane-equations for the top plane
		allDVals = list()
		for x in self._singleLayerCell.cartCoords:
			allDVals.append( abPlaneEquation.calcDForInpXyz(x[:3]) )
		maxD = max(allDVals)
		topPlaneEquation = planeEqn.ThreeDimPlaneEquation( *(abPlaneEquation.coeffs[:3] + [maxD]) )

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
		surfPlaneEquation = planeEqn.ThreeDimPlaneEquation.fromTwoPositionVectors(lattVects[0],lattVects[1],normaliseCoeffs=True)
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


def _getSingleLayerHcp1010FromPrimitiveCell_longTermination(primCell):
	""" Gets the hcp(10-10) primtive cell with long termination; meaning the first interlayer spacing is larger than the second
		WARNING: THIS FUNCTION ISNT WELL TESTED and wasnt written to neccesarily be super-general. Be VERY careful using it """
	startPrimCell = getSingleLayerHcp1010FromPrimitiveCell(primCell)

	assert len(startPrimCell.cartCoords)==2

	#We need to get the image atom for our current top layer
	allCart = startPrimCell.cartCoords
	atomACart, atomBCart = [np.array(x[:3]) for x in allCart]
	lattVectC = np.array( startPrimCell.lattVects[2] )

	atomDistsFromAbPlane = _getDistancesFromAtomsToAbPlaneForInpUCell(startPrimCell)


	if atomDistsFromAbPlane[0] > atomDistsFromAbPlane[1]:
		topAtomIdx, topAtomXyz = 0, allCart[0][:3]
	else:
		topAtomIdx, topAtomXyz = 1, allCart[1][:3]

	#Now apply a translation vector. I THINK this will ALWAYS be the negative direction based on the distance from ab plane
	#Might not be true if we allow starting fract coords that have any values > 1 or < 0
	cVect = startPrimCell.lattVects[2] 
	newBottomAtomXyz = [atomCoord-translation for atomCoord,translation in it.zip_longest(topAtomXyz,cVect)]

	#Now modify the cartesian coordinates
	allCart[topAtomIdx][:3] = newBottomAtomXyz
	startPrimCell.cartCoords = allCart
	return startPrimCell


#NOT TESTED PROPERLY
def _getDistancesFromAtomsToAbPlaneForInpUCell(inpStruct):
    lattVects = inpStruct.lattVects
    bottomPlane = planeEqn.ThreeDimPlaneEquation.fromTwoPositionVectors(lattVects[0], lattVects[1])
    cartCoords = [x[:3] for x in inpStruct.cartCoords]
    distsFromPlane = list()
    for x in cartCoords:
        currDist= bottomPlane.getDistanceOfPointFromPlane(x)
        distsFromPlane.append(currDist)
    return distsFromPlane



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
	angles = inpCell.getLattAnglesList()
	expAngle = 60.0
	if not all([ abs((x-expAngle))<angleTol for x in angles]):
		printMsg = "Angles should also be 60.0 degrees for rocksalt primitive cell, found angles of {}".format(angles)
		isPrim=False

	#Check correct lattice paramter relationships
	lattParams = inpCell.getLattParamsList()
	if not all([ abs(x-lattParams[0]) < lattParamTol for x in lattParams]):
		printMsg = "Lattice parameters should all be equal for rocksalt primitive cell, found parameters of {}".format(lattParams)
		isPrim=False

	if printError:
		print(printMsg)

	return isPrim




def getSingleLayerRocksalt011FromPrimitiveCell(primCell):
	""" 
	
	Args:
		primCell: (plato_pylib UnitCell object) This must be a primitive rock salt cell. Meaning 2 atoms, with each lattice angle=60 degrees and a=b=c
			 
	Returns
		 outCell: (plato_pylib UnitCell object) Contains a 4-atom unit-cell which forms the basis for forming rock-salt 011 surfaces (it is one layer)
 
	Raises:
		 AssertionError: If the input primitive cell is not a rock salt cell. Probably not garanteed to be raised if an incompatible cell is passed, but at least catches some errors.
	"""
	assert _uCellIsRockSaltPrimitive(primCell, printError=True), "UnitCell {} is not a primitive rock salt cell".format(primCell)
	atomA, atomB = primCell.fractCoords[0][-1], primCell.fractCoords[1][-1]
	conventLattParam = (2*primCell.lattParams["a"]) / math.sqrt(2)
	otherLattParams = primCell.lattParams["a"]
	outputLattParams = [conventLattParam, otherLattParams, otherLattParams]
	outputLattAngles = [90,90,90]
	fractCoords = [ [0.0, 0.0, 0.0, atomA],
	                [0.0, 0.5, 0.5, atomB],
	                [0.5, 0.0, 0.0, atomB],
	                [0.5, 0.5, 0.5, atomA] ]
	outCell = uCell.UnitCell(lattParams=outputLattParams, lattAngles=outputLattAngles)
	outCell.fractCoords = fractCoords
	return outCell


#I WAS going to make the surface via this and some general implementation; but gave up basically
def _getConventionalRocksaltCellFromPrimitiveCell(primCell):
	assert _uCellIsRockSaltPrimitive(primCell, printError=True), "UnitCell {} is not a primitive rock salt cell".format(primCell)
	atomA = primCell.fractCoords[0][-1]
	atomB = primCell.fractCoords[1][-1]
	lattParamAll = (2*primCell.lattParams["a"]) / math.sqrt(2)
	outLattParams = [lattParamAll for x in range(3)]
	outLattAngles = [90,90,90]
	fractCoords = [ [0.0, 0.0, 0.0, atomA],
	                [0.5, 0.0, 0.0, atomB],
	                [0.0, 0.5, 0.0, atomB],
	                [0.5, 0.5, 0.0, atomA],
	                [0.0, 0.0, 0.5, atomB],
	                [0.5, 0.0, 0.5, atomA],
	                [0.0, 0.5, 0.5, atomA],
	                [0.5, 0.5, 0.5, atomB] ]
	outCell = uCell.UnitCell(lattParams=outLattParams, lattAngles=outLattAngles)
	outCell.fractCoords = fractCoords
	return outCell



def getSingleLayerBrucite0001FromPrimitiveCell(primCell):
	""" Get a single layer of brucite 0001 (with OH terminations) from a primitive cell (which may well have Mg/OH terminations)
	
	Args:
		primCell: (plato_pylib UnitCell object) This must be a primitive brucite hexagonal cell. Meaning Mg(OH)2 in fract coords; and a=b!=c, alpha=beta=90, gamma=120
			 
	Returns
		 outCell: (plato_pylib UnitCell object) Contains a unit-cell which forms the basis for forming brucite 0001 surfaces (it is one layer)
 
	Raises:
		 AssertionError: If the input primitive cell is not a brucite hexagonal cell. Probably not garanteed to be raised if an incompatible cell is passed, but at least catches some errors.
	"""
	assert _uCellIsBrucitePrimitive(primCell), "Input cell is not a brucite primitive cell"

	#Step 1 is to find the top OH positions
	outCell = copy.deepcopy(primCell)
	topHPos,topHIdx = _findLowestOrHighestZPosAndIdxForElementInFractCoords(outCell.fractCoords, "H", lowOrHigh="high")
	topOPos,topOIdx = _findLowestOrHighestZPosAndIdxForElementInFractCoords(outCell.fractCoords, "O", lowOrHigh="high")
	mgPos,mgIdx = _findLowestOrHighestZPosAndIdxForElementInFractCoords(outCell.fractCoords, "Mg", lowOrHigh="high")
	bottomHPos, bottomHIdx = _findLowestOrHighestZPosAndIdxForElementInFractCoords(outCell.fractCoords, "H", lowOrHigh="low")
	bottomOPos, bottomOIdx = _findLowestOrHighestZPosAndIdxForElementInFractCoords(outCell.fractCoords, "O", lowOrHigh="low")

	#If Mg is the the middle of the cell then do nothing
	if (mgPos < topHPos) and (mgPos < topOPos) and (mgPos > bottomHPos) and (mgPos > bottomOPos):
		return outCell

	#If Mg is at the top of the cell; push it to the bottom of the cell
	if (mgPos > topOPos):
		translationVector = outCell.lattVects[-1]
		cartCoords = copy.deepcopy(outCell.cartCoords) #Copy should be unnecesary really
		cartCoords[mgIdx] = [x-t for x,t in it.zip_longest(cartCoords[mgIdx][:3],translationVector)] + ["Mg"]
		outCell.cartCoords = cartCoords

	#Next step is to displace by the c-vector (i.e. get their images). This SHOULD be the same as subtracting c from the z-coord
	#TODO: The translation vector is not garanteed to be POSITIVE; need some way of testing if i should add or subtract
	translationVector = outCell.lattVects[-1]
	cartCoords = copy.deepcopy(outCell.cartCoords) #Copy should be unnecesary really
	cartCoords[topHIdx] = [x-t for x,t in it.zip_longest(cartCoords[topHIdx][:3],translationVector)] + ["H"]
	cartCoords[topOIdx] = [x-t for x,t in it.zip_longest(cartCoords[topOIdx][:3],translationVector)] + ["O"]
	outCell.cartCoords = cartCoords

	_centreCFractCoordsForInpCell(outCell)
	
	return outCell


def _findLowestOrHighestZPosAndIdxForElementInFractCoords(fCoords, element, lowOrHigh="low"):
	eleList = [fCoord[-1] for fCoord in fCoords]
	zPositions = [fCoord[-2] for fCoord in fCoords]
	topZ, bottomZ = min(zPositions)-0.1, max(zPositions)+0.1 #Start values must allow ALL atoms to be < or > than as required

	for idx,(zPos,ele) in enumerate( zip(zPositions,eleList) ):
		if lowOrHigh.lower()=="low":
			if (ele.lower() == element.lower()) and (zPos<bottomZ):
				bottomZ, outIdx = zPos, idx
		elif lowOrHigh.lower()=="high":
			if (ele.lower() == element.lower()) and (zPos>topZ):
				topZ, outIdx = zPos, idx
		else:
			raise ValueError("{} is an invalid option for lowOrHigh".format(lowOrHigh))
	
	outZ = bottomZ if lowOrHigh.lower()=="low" else topZ
	return outZ,outIdx



def _uCellIsBrucitePrimitive(inpCell, printError=True,  angleTol=1e-1, lattParamTol=1e-1):
	isPrim=True
	printMsg = ""

	#Check angles 
	angles = [x for x in inpCell.getLattAnglesList()]
	expAngles = [90,90,120]
	for exp,act in it.zip_longest(expAngles, angles):
		if abs(exp-act)>angleTol:
			printMsg = "Angles should be {} for Brucite primitive cell; but found angles of {}".format(expAngles,angles)
			isPrim=False

	#Check lattice parameters (a=b!=c always true for the primitive cell)
	a,b,c = [x for x in inpCell.getLattParamsList()]
	if (abs(b-a)>lattParamTol) or (c<=b):
		printMsg = "Lattice parameters {} are inconsistent with a hexagonal cell where a=b!=c".format([a,b,c])
		isPrim=False

	#Check correct elements are present
	expElements = sorted( ["H","H","O","O","Mg"] )
	actElements = sorted( [x[-1] for x in inpCell.fractCoords] )
	if expElements != actElements:
		printMsg = "Elements {} are incorrect for brucite".format(actElements)
		isPrim=False

	if (printError is True) and (isPrim is False):
		print(printMsg)

	return isPrim




