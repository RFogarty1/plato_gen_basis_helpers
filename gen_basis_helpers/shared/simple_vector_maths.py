
import itertools as it
import math

import numpy as np

def getUnitVectorFromInpVector(inpVector):
	lenVect = getLenOneVector(inpVector)
	return [x/lenVect for x in inpVector]

def getLenOneVector(vectA):
	return math.sqrt( sum([x**2 for x in vectA]) )

def getDistTwoVectors(vectA,vectB):
	sqrDiff = [ (a-b)**2 for a,b in it.zip_longest(vectA,vectB) ]
	return math.sqrt( sum(sqrDiff) )

def getAngleTwoVectors(vectA,vectB):
	normFactorA = math.sqrt( sum( [x**2 for x in vectA] ) )
	normFactorB = math.sqrt( sum( [x**2 for x in vectB] ) )

	normA = [x/normFactorA for x in vectA]
	normB = [x/normFactorB for x in vectB]

	dotProd = getDotProductTwoVectors(normA,normB)
	return math.degrees( math.acos(dotProd) )

def getDotProductTwoVectors(vectA,vectB):
	return sum( [a*b for a,b in it.zip_longest(vectA,vectB)] )


#https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
def getRotationMatrixLinkingTwoUnitVectors(uVectA, uVectB):
	assert abs(1 - getLenOneVector(uVectA))<1e-3
	assert abs(1 - getLenOneVector(uVectB))<1e-3
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

	return rotMatrix


def getMatrixForMultipleRotations(axes, angles):
	""" Gets a single rotation matrix representing the input rotations along axes and angles. Note the order rotations are applied are the same as they are input (the first element in axes/angles is the first rotation applied)
	
	Args:
		axes: (iter of len-3 iters) Each element is an axis of rotation
		angles: (iter of floats) Each element is a float
			 
	Returns
		outMatrix: (3x3 "matrix"; really iter of iters) Rotation matrix for the combined rotations requested
 
	"""
	outMatrix = np.identity(3)
	for axis,angle in it.zip_longest(axes,angles):
		outMatrix = np.array( getRotationMatrixAroundAxis(axis,angle) ) @ outMatrix

	return outMatrix

#simply from wikipedia.
def getRotationMatrixAroundAxis(axis, angle):
	""" Gets matrix representation for a rotation around a given axis. Note rotations appear counter-clockwise when looking from end of axis to origin using right-hand rule co-ordinate system
	
	Args:
		axis: (len-3 float iter) Vector defining axis to rotate around, e.g [1,0,0] for x-axis
		angle: (float) Angle (in degrees) to rotate around

	Returns
		 rotMatrix: (3x3 np array), use rotMatrix@geom to get the geometry rotated around this axis (assuming geom has 1 co-ord per column)
 
	"""

	ux, uy, uz = getUnitVectorFromInpVector( axis )
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
	return rotationMatrix


def applyMultipleRotationsAroundAxisToInpCoords(inpCoords, axes, angles, inPlace=True):
	""" Applies multiple rotations to a matrix
	
	Args:
		inpCoords: (nx3 "matrix"; can be an iter of iters)
		axes: (iter of len-3 iters) Each represents an axis to rotate around
		angles: (iter of floats) Angle to rotate around this axis
		inPlace: (Bool) Whether to modify the inpCoords in place
 
	Returns
		outMatrix: (nx3 "matrix; really iter of iters) Contains the rotated co-ordinates
 
	"""
	outMatrix = inpCoords
	for axis, angle in it.zip_longest(axes,angles):
		outMatrix = applyRotationAroundAxisToInpCoords(outMatrix, axis, angle, inPlace=False)

	if inPlace:
		for rowIdx in range(len(outMatrix)):
			for colIdx in range(len(outMatrix[rowIdx])):
				inpCoords[rowIdx][colIdx] = outMatrix[rowIdx][colIdx]

	return outMatrix

def applyRotationAroundAxisToInpCoords(inpCoords, axis, angle, inPlace=True):
	""" Apply a rotation to a set of co-ordinates; either in place or return the output matrix
	
	Args:
		inpCoords: (nx3 "matrix"; can be an iter of iters)
		axis: (len-3 iter) The axis of rotation
		angle: (float) The angle to rotate around this
		inPlace: (Bool) Whether to modify the inpCoords in place
			 
	Returns
		 outMatrix: (nx3 "matrix; really iter of iters) Contains the rotated co-ordinates
 
	"""
	rotationMatrix = getMatrixForMultipleRotations([axis],[angle])
	outMatrixTranspose = rotationMatrix @ np.array(inpCoords).transpose()
	outMatrix = outMatrixTranspose.transpose().tolist()

	if inPlace:
		for rowIdx in range(len(outMatrix)):
			for colIdx in range(len(outMatrix[rowIdx])):
				inpCoords[rowIdx][colIdx] = outMatrix[rowIdx][colIdx]

	return outMatrix


def getStandardRotationAnglesFromRotationMatrix(rotMatrix, sinThetaTol=1e-2):
	""" Gets three angles representing rotations of rotMatrix. This involves refactoring the rotMatrix into R = R_zR_yR_x where we define axes x=[1,0,0], y=[0,-1,0], z=[0,0,1]
	
	Args:
		rotMatrix: (3x3 matrix) Rotation matrix. Numpy array or iter of iters should work.
		sinThetaTol: (float) The maximum the sin(theta) matrix element can be above 1 or below -1. This is essentially to account for numerical errors; we will assume the matrix element is exactly 1 (i.e. theta=90 degrees) in this case
 
	Returns
		outAngles: (len-3 iter of floats) [theta_x, theta_y, theta_z]; These may also be called [roll, pitch, azimuthal]

	Edge Cases/ Restrictions:
		a) ThetaY has domain [-90,90] inclusive (arcsin used to determine). ThetaX and ThetaZ have domain -180<theta<=180 (atan2 used to determine)
		b) We use arctan2 to get theta_x and theta_z. Unlike normal arctan this works for theta=\pm 90 degrees. We pass both sin(theta) and cos(theta) to this function; if sin(theta) is +ve and cos(theta) is 0 we get +90 degrees, if sin(theta) is -ve and cos(theta) is 0 we get -90 degrees. NOTE: Normal math.tan actually works here anyway if rotMatrix is a np array (i think)
 
	"""

#		c) If theta_y is \pm 90 degrees then we cant uniquely solve for theta_x and theta_z; only their sum or difference can be determined. Hence in this case (specifically when the matrix element which should be sin(theta) is >=1 or <=-1) we set theta_z to zero and solve for theta_x. This is somewhat arbitrary but this edge case shouldnt occur often regardless

	#Get the thetaY case; checking for domain
	sinTheta = rotMatrix[2][0]
	if sinTheta > 1:
		if abs(sinTheta - 1) < sinThetaTol:
			thetaY = 90
		else:
			raise ValueError("sinTheta = {} is not allowed".format(sinTheta))

	elif sinTheta < -1:
		if abs(sinTheta + 1) < sinThetaTol:
			thetaY = -90
		else:
			raise ValueError("sinTheta = {} is not allowed".format(sinTheta))

	else:
		thetaY = math.degrees( math.asin(rotMatrix[2][0]) )

	#Get the other two angles
	thetaX = math.degrees( math.atan2(rotMatrix[2][1], rotMatrix[2][2]) )
	thetaZ = math.degrees( math.atan2(rotMatrix[1][0], rotMatrix[0][0]) )

	return [thetaX, thetaY, thetaZ]

