
import copy
import itertools as it
import math

import numpy as np

import plato_pylib.shared.unit_convs as uConvHelp

from . import adsorbate_rep_objs as adsObjHelp

from ..analyse_md import calc_dists as distsHelp
from ..shared import plane_equations as planeEqnHelp
from ..shared import simple_vector_maths as vectHelp


def getStandardWaterAdsObjFromXyz(inpXyz, surfNormal=None):
	""" Creates a WaterAdsorbateStandard object from input xyz co-ords
	
	Args:
		inpXyz: (len-3 iter of len-4 iters) Cartesian co-ords for water. E.g. [ [0,0,0,"O"], [0,0.7,0.7,"H"], [0,0.7,-0.7,"H"] ]
		surfNormal: (len-3 iter) Normal vector for the surface, used to define the azimuthal axis. Only ever tested with default of [0,0,1]
 
	Returns
		 outObj: (WaterAdsorbateStandard) 
 
	"""
	surfNormal = [0,0,1] if surfNormal is None else surfNormal

	dists, angle = getWaterInternalCoordsFromXYZ(inpXyz)
	surfNormal = [0,0,1]
	roll, pitch, azimuth = getStandardAxisRotationsFromXyz(inpXyz, surfNormal=surfNormal)
	refPosGetter = WaterRefPosGetterStandard()
	oIndices = [idx for idx,x in enumerate(inpXyz) if x[-1].upper()=="O"]
	assert len(oIndices)==1
	translationVector = inpXyz[oIndices[0]][:3]

	kwargDict = {"refPosGetter":refPosGetter, "pitch":pitch, "roll":roll,
	             "azimuthal":azimuth, "translationVector": translationVector}
	outObj = WaterAdsorbateStandard(dists, angle, **kwargDict)

	return outObj


class WaterAdsorbateStandard(adsObjHelp.Adsorbate):

	def __init__(self, ohDists, angle, refPosGetter=None,
	             pitch=None, roll=None, azimuthal=None, translationVector=None):
		""" Initialiser
		
		Args:
			ohDists: (len-2 float iter) Distances between the two O-H atoms
			angle: (float) The H-O-H angle
			refPosGetter: f(ohDists,angle) Returns a set of x/y/z coords from internal coordinates. Also needs properties defining rollAxis, pitchAxis and azimuthalAxis
			pitch: (Float, Optional, default=0) Rotation around an axis orthogonal to the roll/azimuthal axes, essentially used to tilt O-H up or down. Generally +ve values will tilt OH upwards ([0,-1,0] pitch axis)
			roll : (Float, Optional, default=0) Rotation around the bisect of the two OH/OH bonds in the reference geometry. Generally expected to be ~0 or ~180
			azimuthal: (Float, Optional, default=0) Rotation about the z-axis
			translationVector: (Float, len-3 iter, default=[0,0,0]): Apply this translation to all atoms 
		"""
		self._eqTol = 1e-5
		self.ohDists = list(ohDists)
		self.angle = angle
		self.refPosGetter = refPosGetter if refPosGetter is not None else WaterRefPosGetterStandard()
		self.pitch = pitch if pitch is not None else 0
		self.roll = roll if roll is not None else 0
		self.azimuthal = azimuthal if azimuthal is not None else 0
		self.translationVector = list(translationVector) if translationVector is not None else [0.0,0.0,0.0]

	@property
	def geom(self):
		outGeom = self.refPosGetter(self.ohDists, self.angle)
		self._applyRollToOutgeom(outGeom)
		self._applyPitchToOutGeom(outGeom)
		self._applyAzimuthalToOutGeom(outGeom)
		self._applyTranslationToOutGeom(outGeom)
		return outGeom

	#Roll is along the x-axis by default, but more generally its along the bisect angle really.
	#REALLY we can have a general function for this
	def _applyRollToOutgeom(self, geom):
		self._applyRotationToGeom(geom, self.refPosGetter.rollAxis, self.roll)

	def _applyPitchToOutGeom(self, geom):
		self._applyRotationToGeom(geom, self.refPosGetter.pitchAxis, self.pitch)

	def _applyAzimuthalToOutGeom(self, geom):
		self._applyRotationToGeom(geom, self.refPosGetter.azimuthalAxis, self.azimuthal)

	def _applyRotationToGeom(self, geom, axis, angle):
		outGeom = getGeomRotatedAroundAxis(geom, axis, angle)
		self._moveNewCoordsToGeomInPlace(outGeom, geom)

	def _applyTranslationToOutGeom(self, geom):
		newCoords = list()
		for oldCoords in geom:
			currCoords = [x+t for x,t in it.zip_longest(oldCoords[:3],self.translationVector)]
			newCoords.append(currCoords)
		self._moveNewCoordsToGeomInPlace(newCoords, geom)

	def _moveNewCoordsToGeomInPlace(self, newCoords, geom):
		for idx, coords in enumerate(newCoords):
			geom[idx][:3] = coords[:3]

	def __eq__(self, other):
		eqTol = min(self._eqTol, other._eqTol)
		floatAttrs = ["angle", "pitch", "roll", "azimuthal"]
		floatListAttrs = ["ohDists", "translationVector"]

		if self.refPosGetter != other.refPosGetter:
			return False

		for attr in floatAttrs:
			valA, valB = getattr(self, attr), getattr(other,attr)
			if abs(valA-valB) > eqTol:
				return False

		for attr in floatListAttrs:
			valsA, valsB = getattr(self,attr), getattr(other,attr)
			diffs = [abs(x-y) for x,y in it.zip_longest(valsA,valsB)]
			if any([x>eqTol for x in diffs]):
				return False

		return True


class WaterRefPosGetterStandard():
	""" Takes internal co-ordinates for water and converts into a set of cartesian co-ordinates. The callable interface is as follows:

	Args:
		ohDists: (len-2 iter) O-H distances
		angle: (float) H-O-H angle in degrees

	"""
	def __init__(self):
		self._eqTol = 1e-5

	@property
	def rollAxis(self):
		return [1,0,0]

	@property
	def pitchAxis(self):
		return [0,-1,0]

	@property
	def azimuthalAxis(self):
		return [0,0,1]

	def __call__(self, ohDists, angle):
		outDists = sorted(ohDists)
		outSymbols = ["O","H","H"]
		usedAngle = angle / 2
		xVals = [math.cos(math.radians(usedAngle))*d for d in outDists]
		yVals = [math.sin(math.radians(usedAngle))*d for d in outDists]
		yVals[-1] *= -1
		outVals = [ [0.0,0.0,0.0,"O"] ]
		for x,y,symbol in it.zip_longest(xVals,yVals,outSymbols[1:]):
			currVals = [x,y,0,symbol]
			outVals.append(currVals)
		return outVals		

	def __eq__(self, other):
		if type(self) != type(other):
			return False

		eqTol = min(self._eqTol, other._eqTol)
		props = ["rollAxis", "pitchAxis", "azimuthalAxis"]
		for prop in props:
			diffs = [abs(x-y) for x,y in zip(getattr(self,prop), getattr(other,prop))]
			if any([x>eqTol for x in diffs]):
				return False
		return True


def getWaterInternalCoordsFromXYZ(inpCoords):
	""" Get two OH-distances and the H-O-H angle when given a set of water monomer cartesian co-ordinates
	
	Args:
		inpCoords: (n length iter of len-4 iters) Each entry is [x,y,z,symbol] where symbol is the chemical symbol
			 
	Returns
		outDists: (Len-2 float iter) The two O-H bondlengths
 		angle: (float) The H-O-H angle

	Raises:
		 ValueError: If input co-ordinates dont contain one oxygen and two hydrogen
	"""
	_checkInputCoordsAreH2O(inpCoords)
	oCoord = [x[:3] for x in inpCoords if x[-1].upper()=="O"][0]
	hCoords = [x[:3] for x in inpCoords if x[-1].upper()=="H"]
	vectA = [b-a for b,a in it.zip_longest(hCoords[0],oCoord)]
	vectB = [b-a for b,a in it.zip_longest(hCoords[1],oCoord)]

	outDists = [vectHelp.getDistTwoVectors(oCoord, x) for x in hCoords]
	outAngle = vectHelp.getAngleTwoVectors(vectA,vectB)
	return outDists, outAngle

def _checkInputCoordsAreH2O(inpCoords):
	""" raises ValueError if input coords dont correspond to a H2O chemical formula
	
	Args:
		inpCoords: (n length iter of len-4 iters) Each entry is [x,y,z,symbol] where symbol is the chemical symbol

	Raises:
		 ValueError: If inpCoords do not correspond to h2o
	"""
	if len(inpCoords) != 3:
		raise ValueError("water coords are length {}, should be length 3".format( len(inpCoords) ))

	numbH = len( [x[-1] for x in inpCoords if x[-1].upper()=="H"] )
	numbO = len( [x[-1] for x in inpCoords if x[-1].upper()=="O"] )
	if ( (numbH != 2) or (numbO != 1) ):
		raise ValueError("H2O expected, but got H{}O{}".format(numbH,numbO))



#simply from wikipedia.
def getGeomRotatedAroundAxis(inpGeom, axis, angle):
	""" Note rotations appear counter-clockwise when looking from end of axis to origin using right-hand rule co-ordinate system
	
	"""
	rotationMatrix = vectHelp.getRotationMatrixAroundAxis(axis,angle)
	geomToRotate = np.array( [x[:3] for x in inpGeom] ).transpose()
	rotated = np.dot(rotationMatrix, geomToRotate)
	outGeom = [list(x[:3]) for x in rotated.transpose()]
	return outGeom



def getStandardAxisRotationsFromXyz(inpCoords, surfNormal=None):
	""" Gets roll, pitch and azimuthal angles that relate a standard water orientation to the input co-ordinates
	
	Args:
		inpCoords: (n length iter of len-4 iters) Each entry is [x,y,z,symbol] where symbol is the chemical symbol
		surfNorm: (len 3-vector) Normal to the surface; needed to define azimuthal axis uniquely (it ALWAYS points the same direction as the surface normal OR orthogonal to it), Default [0,0,1] 
 
	Returns (All values in degrees with domain -90 to +90)
		roll: (float) 
		pitch: (float)
		azimuth: (float)
	Raises:
		ValueError: If any angles magnitude is close to 90 degrees; the function becomes badly behaved in this case
 
	"""
	refPosGetter =  WaterRefPosGetterStandard()
	#1) Need to get a "standard" orientation
	ohDists, angle = getWaterInternalCoordsFromXYZ(inpCoords)
	standardOrientation = refPosGetter(ohDists, angle)

	#2) Get the input geometry in terms of O,H,H order
	reorderedCoords = list()
	for x in inpCoords:
		if x[-1] == "O":
			reorderedCoords.append(x)
	for x in inpCoords:
		if x[-1] == "H":
				reorderedCoords.append(x)

	#3) Check the order of elements is the same in standardOrientation and out coords
	eleListA, eleListB = [x[-1] for x in reorderedCoords], [x[-1] for x in standardOrientation]
	assert eleListA==eleListB

	#4) Get the new rotated co-ordinate system
	rollVector = vectHelp.getUnitVectorFromInpVector( _getBisectVectorForInpWaterCoords(inpCoords) )
	azimuthVector = _getAzimuthVectorForInpWaterCoords(inpCoords, surfVector=surfNormal)
	pitchVector = getGeomRotatedAroundAxis([rollVector], azimuthVector, -90)[0]


	#5) Get the rotation matrix linking the two co-ordinate systems
	startVectors = np.array([refPosGetter.rollAxis, refPosGetter.pitchAxis, refPosGetter.azimuthalAxis])
	finalVectors = np.array([x[:3] for x in [rollVector,pitchVector,azimuthVector]]).transpose()
	rotMatrix = np.dot( finalVectors, np.linalg.inv(startVectors) )

	#6) Get the roll, pitch, azimuth from the rotatation matrix (which we know is R_az@R_pi@R_roll, hence we have an analyitcal form for each matrix element)
	roll, pitch, azimuth = 0,0,0
	roll = math.degrees( math.atan( rotMatrix[2,1]/rotMatrix[2,2] ) )
	pitch = math.degrees( math.asin( rotMatrix[2,0] ) )	
	azimuth = math.degrees( math.atan( rotMatrix[1,0]/rotMatrix[0,0] ) ) 

	#Cant deal with values near 90 at the moment; may try and sort later but for now just throw
	for angle in roll,pitch,azimuth:
		if abs( abs(angle)-90 )  < 1e-2:
			raise ValueError("[roll,pitch,azimuth] = {}. One of these is too close to 90 degrees".format([roll,pitch,azimuth]))


	return roll,pitch,azimuth

def _getAzimuthVectorForInpWaterCoords(inpCoords, surfVector=None):
	surfVector = [0,0,1] if surfVector is None else surfVector
	inPlaneVects = _getOHVectorsFromInpCoords(inpCoords)
	planeEqn = planeEqnHelp.ThreeDimPlaneEquation.fromTwoPositionVectors(*inPlaneVects)
	normToMolecularPlane = planeEqn.coeffs[:3]
	if vectHelp.getDotProductTwoVectors(surfVector,normToMolecularPlane) < 0:
		normToMolecularPlane = [-1*x for x in normToMolecularPlane]
	return normToMolecularPlane

def _getBisectVectorForInpWaterCoords(inpCoords):
	vectA, vectB = _getOHVectorsFromInpCoords(inpCoords)
	vectA, vectB = [vectHelp.getUnitVectorFromInpVector(x) for x in [vectA,vectB]]
	outVect = [a+b for a,b in it.zip_longest(vectA,vectB)]
	return outVect

def _getOHVectorsFromInpCoords(inpCoords):
	oCoord = [x[:3] for x in inpCoords if x[-1].upper()=="O"][0]
	hCoords = [x[:3] for x in inpCoords if x[-1].upper()=="H"]
	vectA = [b-a for b,a in it.zip_longest(hCoords[0],oCoord)]
	vectB = [b-a for b,a in it.zip_longest(hCoords[1],oCoord)]
	return [vectA, vectB]


def getWaterAdsorptionObjsFromInpCellAndWaterIndices(inpCell, indices):
	""" Gets a list of WaterAdsorbateStandard from an inputCell and indices of water molecules
	
	Args:
		inpCell: (plato_pylib UnitCell object)
		indices: (iter of len-3 iters) Each contains indices for a single water molecule
			 
	Returns
		 outObjs: (iter of WaterAdsorbateStandard objs) The cartesian co-ords for the Oxygen will be unchanged, the Hydrogen will be the images closest to the oxygen 
 
	"""
	lattVects = inpCell.lattVects
	cartCoords = inpCell.cartCoords

	outObjs = list()
	for currIndices in indices:
		currStartCoords = [cartCoords[idx] for idx in currIndices]
		currXYZ = _getOCentredWaterCoords_minImageConv(inpCell, currStartCoords)
		currObj = getStandardWaterAdsObjFromXyz(currXYZ)
		outObjs.append(currObj)

	return outObjs

def _getOCentredWaterCoords_minImageConv(inpCell, waterCoords):
	#Get indices for O/H; so that we can centre on O
	oIndices = [idx for idx,x in enumerate(waterCoords) if x[-1].upper()=="O"]
	hIndices = [idx for idx,x in enumerate(waterCoords) if x[-1].upper()=="H"]
	assert len(oIndices)==1
	assert len(hIndices)==2
	
	#Get co-ords by minimum image convention (i think....)
	oCoords = waterCoords[oIndices[0]]
	outCoords = [ oCoords ]
	for idx in hIndices:
		currStartCoords = waterCoords[idx]
		currOutCoords = distsHelp.getNearestImageNebCoordsBasic( inpCell, oCoords[:3], waterCoords[idx][:3] ) + ["H"]
		outCoords.append(currOutCoords)

	return outCoords


#TODO: Add topLayerAdsObj=None. Allow user to pass an adsorbate obj (ignoring translation vector) that will replace the 
def	getAdsorbateObjsForNextWaterBilayerBasic(inpAds, bilayerSpacing, layerTolerance=0.2*uConvHelp.ANG_TO_BOHR, top=True, surfNormal=None, topLayerTemplate=None):
	""" Gets a list of adsorbate objects which form the next water bilayer for a structure
	
	Args:
		inpAds: (iter of WaterAdsorbateStandard objs) Each represents one water molecule in the top bilayer (we add either above/below this)
		bilayerSpacing: (float) The closest gap between OXYGEN atoms between two bilayers
		layerTolerance: (float) The maximum distance between two oxygens for them to be considered in the same single layer of a bilayer [should be set pretty small]
		surfNormal: [len-3 iter] [x,y,z] for the surface normal vector. Default = [0,0,1]
		topLayerTemplate: (WaterAdsorbateStandard obj) By default the code will create the new layers by transforming the two types of water adsorbates (flat/Hup or Hdown) in the inpAds objects. If you set topLayerTemplate the code will act as if THIS adsorbate was found at the top of the inpAds objs bilayer structure. This is useful for creating extra adsorbate layers on HDown structures but using Hup configurations for the next layers
	
	NOTE:
		All inpsAds objects must be using the default WaterRefPosGetterStandard, which has the oxygen at origin in the "standard" orientation

	Returns
		outAdsObjs: (iter of WaterAdsorbateStandard objs) 
 
	Raises:
		 ValueError: If inpAds span more than 2 layers. Layers are defined using layerTolerance
	"""

	surfNormal = [0,0,1] if surfNormal is None else vectHelp.getUnitVectorFromInpVector(surfNormal)

	#1) Figure out the average height between each "monolayer" making up the bilayer
	oCoords = list()
	for adsObj in inpAds:
		currCoords = [x for x in adsObj.geom if x[-1].upper()=="O"]
		assert len(currCoords)==1
		oCoords.append(currCoords[0])

	planeEqn = planeEqnHelp.ThreeDimPlaneEquation(*surfNormal,0)
	signedODists = [planeEqn.getSignedDistanceOfPointFromPlane(x[:3]) for x in oCoords]
	botDist, topDist = min(signedODists), max(signedODists)

	#1.5) Divide the bilayer into a "top" and "bottom" layer
	botIndices, topIndices = list(), list()
	for idx,dist in enumerate(signedODists):
		distFromBot, distFromTop = abs(dist-botDist), abs(dist-topDist)
		if distFromBot<layerTolerance:
			botIndices.append(idx)
		if distFromTop<layerTolerance:
			topIndices.append(idx)
		if (distFromBot<layerTolerance) and (distFromTop<layerTolerance):
			raise ValueError("Oxygen atom couldnt be assigned to either \"top\" or \"bottom\" of the bilayer. This likely means layerTolerance is set too high")

	#TODO: Change to a ValueError
	totalIndices = len(botIndices) + len(topIndices)
	if totalIndices != len(signedODists):
		raise ValueError("Assigned {} indices despite {} water present; This is likely due to either problems with the number of layers OR layerTolerance is set wrong".format(totalIndices,len(signedODists)))

	#1.9) NOW we figure out the average gap between "top" and "bottom" layers
	avTopLayerSignedDist = sum([signedODists[idx] for idx in topIndices]) / len(topIndices)
	avBotLayerSignedDist = sum([signedODists[idx] for idx in botIndices]) / len(botIndices)
	layerHeight = abs(avTopLayerSignedDist-avBotLayerSignedDist)

	#2.5) Get the typical "top" and "bottom" orientations. They should all be pretty similar regardless.
	#This is obviously a very imperfect way to do this; User supplying the "bottom" and "top" layer orientations will probably be the normal 
	#way regardless though
	topLayerTemplate = copy.deepcopy( inpAds[topIndices[0]] ) if topLayerTemplate is None else topLayerTemplate
	botLayerTemplate = copy.deepcopy( inpAds[botIndices[0]] )

	#3) Apply the relevant translations
	if top:
		firstLayerTVect = [bilayerSpacing*x for x in surfNormal]
		secondLayerTVect = [(bilayerSpacing+(2*layerHeight))*x for x in surfNormal]
		outAdsObjs = list()
		#First sort out the closest layer
		for idx in topIndices:
			startTVect = inpAds[idx].translationVector
			currAdsObj = copy.deepcopy(botLayerTemplate) #NEarest to top layer, so orientation flips to "bot-layer like"
			currAdsObj.translationVector = [x+y for x,y in zip(startTVect, firstLayerTVect)]
			currAdsObj.azimuthal += 180
			outAdsObjs.append(currAdsObj)

		for idx in botIndices:
			startTVect = inpAds[idx].translationVector
			currAdsObj = copy.deepcopy(topLayerTemplate)
			currAdsObj.translationVector = [x+y for x,y in zip(startTVect, secondLayerTVect)]
			currAdsObj.azimuthal += 180
			outAdsObjs.append(currAdsObj)

	else:
		raise NotImplementedError("")


	return outAdsObjs
