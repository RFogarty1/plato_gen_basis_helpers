
import copy
import itertools as it
import math
import unittest

import gen_basis_helpers.shared.simple_vector_maths as tCode


class TestLenOneVector(unittest.TestCase):

	def testForVectorA(self):
		testVectA = [4.0,3.0,2.0,1.0]
		expLength = math.sqrt(30)
		actLength = tCode.getLenOneVector(testVectA)
		self.assertAlmostEqual(expLength, actLength)

class TestDistanceTwoVectors(unittest.TestCase):

	def testForVectorsAB(self):
		testVectA = [1.0, 2.0]
		testVectB = [3.0, 6.0]
		expDist = math.sqrt(20)
		actDist = tCode.getDistTwoVectors(testVectA, testVectB)
		self.assertAlmostEqual(expDist,actDist)


class TestGetAngleTwoVectors(unittest.TestCase):

	def testForVectorsAB(self):
		testVectA = [3,0,0]
		testVectB = [0,3,4]
		expAngle = 90
		actAngle = tCode.getAngleTwoVectors(testVectA,testVectB)
		self.assertAlmostEqual(expAngle,actAngle)


class TestGetUnitVectorFromInpVector(unittest.TestCase):

	def testForSimpleVectorA(self):
		vectA = [0,3,0]
		expVector = [0,1,0]
		actVector = tCode.getUnitVectorFromInpVector(vectA)
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expVector,actVector)]


class TestApplyRotationAroundAxisToInpCoords(unittest.TestCase):

	def setUp(self):
		self.coordsA = [ [1,1,2],
		                 [2,0,2] ]
		self.axis = [0,0,1]
		self.angle = 90
		self.inPlace= True

	def _runTestFunct(self):
		args = [self.coordsA, self.axis, self.angle]
		kwargs = {"inPlace":self.inPlace}
		return tCode.applyRotationAroundAxisToInpCoords(*args,**kwargs)

	def _checkExpAndActCoordsMatch(self, expCoords, actCoords):
		self.assertEqual( len(expCoords), len(actCoords) )
		for expCoord, actCoord in zip(expCoords,actCoords):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expCoord, actCoord)]

	def testExpectedSimpleInPlace(self):
		expCoords = [ [-1,1,2], [0,2,2] ]
		self._runTestFunct()
		actCoords = self.coordsA
		self._checkExpAndActCoordsMatch(expCoords, actCoords)

	def testExpectedNotInPlace(self):
		startCoords = copy.deepcopy(self.coordsA)
		self.inPlace = False
		expCoords = [ [-1,1,2], [0,2,2] ]
		actCoords = self._runTestFunct()
		self._checkExpAndActCoordsMatch(expCoords, actCoords)
		self._checkExpAndActCoordsMatch(startCoords,self.coordsA) #Check we havent modified the input coords


class TestApplyMultipleRotationsAroundAxes(unittest.TestCase):

	def setUp(self):
		self.coordsA = [ [1,1,2], [2,0,2] ]
		self.axisA = [0,0,1]
		self.axisB = [0,0,1]
		self.angleA = 45
		self.angleB = 45
		self.inPlace = True

	def _runTestFunct(self):
		axes = [self.axisA, self.axisB]
		angles = [self.angleA, self.angleB]
		args = [self.coordsA, axes, angles]
		kwargs = {"inPlace":self.inPlace}
		return tCode.applyMultipleRotationsAroundAxisToInpCoords(*args, **kwargs)

	def _checkExpAndActCoordsMatch(self, expCoords, actCoords):
		self.assertEqual( len(expCoords), len(actCoords) )
		for expCoord, actCoord in zip(expCoords,actCoords):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expCoord, actCoord)]

	def testExpectedSimpleInPlace(self):
		self._runTestFunct()
		expCoords = [ [-1,1,2], [0,2,2] ]
		actCoords = self.coordsA
		self._checkExpAndActCoordsMatch(expCoords, actCoords)

	def testExpectedNotInPlace(self):
		startCoords = copy.deepcopy(self.coordsA)
		self.inPlace = False
		expCoords = [ [-1,1,2], [0,2,2] ]
		actCoords = self._runTestFunct()
		self._checkExpAndActCoordsMatch(expCoords, actCoords)
		self._checkExpAndActCoordsMatch(startCoords,self.coordsA) #Check we havent modified the input coords


class TestGetStandardRotationAnglesFromRotationMatrix(unittest.TestCase):

	def setUp(self):
		self.thetaX = 50
		self.thetaY = 60
		self.thetaZ = -35
		self._createRotationMatrix()

	def _runTestFunct(self):
		return tCode.getStandardRotationAnglesFromRotationMatrix(self.rotMatrix)

	def _createRotationMatrix(self):
		axes = [ [1,0,0], [0,-1,0], [0,0,1] ]
		angles = [self.thetaX, self.thetaY, self.thetaZ] 
		self.rotMatrix = tCode.getMatrixForMultipleRotations(axes, angles).tolist()

	def testExpectedNonEdgeCase(self):
		expAngles = [self.thetaX, self.thetaY, self.thetaZ]
		actAngles = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expAngles,actAngles)]

	def testExpected_xAtMinus90Degrees(self):
		""" Likely not really important; Python seems to handle the div by zero error if its in tan anyway """
		self.thetaX = -90
		self._createRotationMatrix()
		self.rotMatrix[2][2] = 0
		self.rotMatrix[2][1] = -0.5
		expAngles = [self.thetaX, self.thetaY, self.thetaZ]
		actAngles = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expAngles,actAngles)]

	def testExpected_thetaY90Degrees_caseA(self):
		""" More specifically asin(theta_y)=1; slightly larger values could occur due to numerical imprecision """
		self.thetaY = 90
		self._createRotationMatrix()
		expAngles = [self.thetaX, self.thetaY, self.thetaZ]
		actAngles = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expAngles,actAngles)]

	def testExpectedThetaYMatrixElementGreaterThanOne(self):
		""" Could occur due to numerical error"""
		self.thetaY = 90
		self._createRotationMatrix()
		self.rotMatrix[2][0] = 1.0001
		expAngles = [self.thetaX, self.thetaY, self.thetaZ]
		actAngles = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expAngles,actAngles)]

	def testExpectedThetaYMatrixElementLessThanMinusOne(self):
		""" Could occur due to numerical error"""
		self.thetaY = -90
		self._createRotationMatrix()
		self.rotMatrix[2][0] = -1.0001
		expAngles = [self.thetaX, self.thetaY, self.thetaZ]
		actAngles = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expAngles,actAngles)]
	

	def testExpected_thetaY90Degrees_caseB(self):
		""" More specifically asin(theta_y)=1; slightly larger values could occur due to numerical imprecision """
		self.thetaY = 90
		self.thetaX += 20
		self.thetaZ += 15
		self._createRotationMatrix()
		expAngles = [self.thetaX, self.thetaY, self.thetaZ]
		actAngles = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expAngles,actAngles)]

	def testExpectedThetaXandZNear180(self):
		self.thetaX = 170
		self.thetaZ = 130
		self._createRotationMatrix()
		expAngles = [self.thetaX, self.thetaY, self.thetaZ]
		actAngles = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expAngles,actAngles)]

	def testExpectedThetaZRotatedBy180(self):
		self.thetaZ = 180
		self._createRotationMatrix()
		expAngles = [self.thetaX, self.thetaY, self.thetaZ]
		actAngles = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expAngles,actAngles)]


	@unittest.skip("Unclear what the result SHOULD be...")
	def testExpected_thetaYOutOfDomain(self):
		self.thetaY = 110
		self._createRotationMatrix()
		expAngles = [self.thetaX, 70, self.thetaZ]
		actAngles = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expAngles,actAngles)]
 

