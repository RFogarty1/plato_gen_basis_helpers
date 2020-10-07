
import itertools as it
import math

import numpy as np

from . import adsorbate_rep_objs as adsObjHelp
from ..shared import simple_vector_maths as vectHelp

class WaterAdsorbateStandard(adsObjHelp.Adsorbate):

	def __init__(self, ohDists, angle, refPosGetter=None,
	             pitch=None, roll=None, azimuthal=None):
		""" Initialiser
		
		Args:
			ohDists: (len-2 float iter) Distances between the two O-H atoms
			angle: (float) The H-O-H angle
			refPosGetter: f(ohDists,angle) Returns a set of x/y/z coords from internal coordinates. Also needs properties defining rollAxis, pitchAxis and azimuthalAxis
			pitch: (Float, Optional, default=0) Rotation around an axis orthogonal to the roll/azimuthal axes, essentially used to tilt O-H up or down. Generally +ve values will tilt OH upwards ([0,-1,0] pitch axis)
			roll : (Float, Optional, default=0) Rotation around the bisect of the two OH/OH bonds in the reference geometry. Generally expected to be ~0 or ~180
			azimuthal: (Float, Optional, default=0) Rotation about the z-axis
 
		"""
		#TODO
		#	translationVector: (Float, len-3 iter): Apply this translation to the water
		self.ohDists = list(ohDists)
		self.angle = angle
		self.refPosGetter = refPosGetter if refPosGetter is not None else WaterRefPosGetterStandard()
		self.pitch = pitch if pitch is not None else 0
		self.roll = roll if roll is not None else 0
		self.azimuthal = azimuthal if azimuthal is not None else 0
#		self.translationVector = list(translationVector) if translationVector is not None else [0.0,0.0,0.0]

	@property
	def geom(self):
		outGeom = self.refPosGetter(self.ohDists, self.angle)
		self._applyRollToOutgeom(outGeom)
		self._applyPitchToOutGeom(outGeom)
		self._applyAzimuthalToOutGeom(outGeom)
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

	def _moveNewCoordsToGeomInPlace(self, newCoords, geom):
		for idx, coords in enumerate(newCoords):
			geom[idx][:3] = coords[:3]

class WaterRefPosGetterStandard():
	""" Takes internal co-ordinates for water and converts into a set of cartesian co-ordinates. The callable interface is as follows:

	Args:
		ohDists: (len-2 iter) O-H distances
		angle: (float) H-O-H angle in degrees

	"""
	def __init__(self):
		pass

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
	ux, uy, uz = vectHelp.getUnitVectorFromInpVector( axis )
	rotationMatrix =  np.zeros( [3,3] )

	cosTheta = math.cos(math.radians(angle))
	sinTheta = math.sin(math.radians(angle))
	rotationMatrix[0][0] = cosTheta + ( (ux**2) * (1-cosTheta) )
	rotationMatrix[0][1] = (ux*uy)*(1-cosTheta) - (uz*sinTheta)
	rotationMatrix[0][2] = (ux*uz)*(1-cosTheta) + (uy*sinTheta)
	rotationMatrix[1][0] = (ux*uy*(1-cosTheta)) + (uz*sinTheta)
	rotationMatrix[1][1] = cosTheta + ((uy**2)*(1-cosTheta))
	rotationMatrix[1][2] = (uy*uz)*(1-cosTheta) - (ux*sinTheta)
	rotationMatrix[2][0] = (uz*ux)*(1-cosTheta) - (uy*sinTheta)
	rotationMatrix[2][1] = (uz*uy)*(1-cosTheta) + (ux*sinTheta)
	rotationMatrix[2][2] = cosTheta + ((uz*uz)*(1-cosTheta))

	geomToRotate = np.array( [x[:3] for x in inpGeom] ).transpose()
	rotated = np.dot(rotationMatrix, geomToRotate)
	outGeom = [list(x[:3]) for x in rotated.transpose()]

	return outGeom





