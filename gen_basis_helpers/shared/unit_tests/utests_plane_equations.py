
import copy
import math
import itertools as it
import unittest
import unittest.mock as mock


import gen_basis_helpers.shared.plane_equations as tCode



class TestThreeDimPlaneEqn(unittest.TestCase):

	def setUp(self):
		a = 3
		b = 0
		c = 2
		d = 4
		self.planeParams = [a,b,c,d]
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.ThreeDimPlaneEquation(*self.planeParams)

	def testCalcDValueForXyz(self):
		inpXyz = [1,2,3]
		expDVal = 3 + 0 + 6
		actDVal = self.testObjA.calcDForInpXyz(inpXyz)
		self.assertEqual( expDVal, actDVal )

	def testCalcDistBetweenPlaneAndPoint_planeAlongXy(self):
		zDist = 4
		pointCoords = [0.0,0.0,zDist]
		self.planeParams = [0.0,0.0,1.0,0.0] #This is the xy plane which intercepts z at zero
		self.createTestObjs()
		expDistance = zDist
		actDistance = self.testObjA.getDistanceOfPointFromPlane(pointCoords)
		self.assertAlmostEqual(expDistance, actDistance)

	def testCalcDistanceBetweenPlaneAndPoint_011_plane(self):
		pointCoords = [0, 0, 0]
		self.planeParams = [0.0, 1.0, 1.0, 1.0] #Intercepts both y and z at 1
		self.createTestObjs()
		#z=0.5 and y=0.5 should be the closest point; hence the expected distance can be figured out easily
		expDist = math.sqrt( (0.5**2) + (0.5**2) )
		actDist = self.testObjA.getDistanceOfPointFromPlane(pointCoords)
		self.assertAlmostEqual(expDist, actDist)

	def testCalcFromTwoPositionVectors_xyPlane(self):
		expCoeffs = [0,0,1,0]
		inpVectors = [ [1,0,0], [0,1,0] ]
		outObj = tCode.ThreeDimPlaneEquation.fromTwoPositionVectors(*inpVectors, normaliseCoeffs=False)
		actCoeffs = outObj.coeffs

		for exp,act in it.zip_longest(expCoeffs,actCoeffs):
			self.assertAlmostEqual(exp,act)

	def testCalcFromTwoPositionVectors_xyPlaneNormalised(self):
		expCoeffs = [0,0,1,0]
		inpVectors = [ [2,0,0],[0,3,0] ]
		outObj = tCode.ThreeDimPlaneEquation.fromTwoPositionVectors(*inpVectors, normaliseCoeffs=True)
		actCoeffs = outObj.coeffs

		for exp,act in it.zip_longest(expCoeffs,actCoeffs):
			self.assertAlmostEqual(exp,act)

	def testCalcSignedDistanceBetweenPlaneAndPoint_011_plane_negative_dist(self):
		pointCoords = [0,0,0]
		self.planeParams = [0, 1.0, 1.0, 1.0]
		self.createTestObjs()
		#z=0.5 and y=0.5 should be the closest point; hence the expected distance can be figured out easily
		unSignedDist = math.sqrt( (0.5**2) + (0.5**2) )
		expDist = -1*unSignedDist #The normal vector points the opposite direction, hence signed distance is negative
		actDist = self.testObjA.getSignedDistanceOfPointFromPlane(pointCoords)
		self.assertAlmostEqual(expDist,actDist)

	def testCalcSignedDistanceBetweenPlaneAndPoint_011_plane_positive_dist(self):
		pointCoords = [0,0,0]
		self.planeParams = [0, -1.0, -1.0, -1.0]
		self.createTestObjs()
		expDist = math.sqrt( (0.5**2) + (0.5**2) )
		actDist = self.testObjA.getSignedDistanceOfPointFromPlane(pointCoords)
		self.assertAlmostEqual(expDist,actDist)

	def testTwoEqualPlanesCompareEqual(self):
		planeA = self.testObjA
		planeB = copy.deepcopy(planeA)
		self.assertEqual(planeA,planeB)

	def testTwoUnEqualPlanesCompareUnEqual(self):
		planeA = self.testObjA
		planeB = copy.deepcopy(planeA)
		planeB.d += 1
		self.assertNotEqual(planeA,planeB)


class TestGetOutOfPlaneDistTwoPoints(unittest.TestCase):

	def setUp(self):
		self.posA = [3, 5, 5]
		self.posB = [3, 7, 8]
		self.planeEqnA = tCode.ThreeDimPlaneEquation(0,0,1,6) #Put in the centre so that we need signed distances

	def _runTestFunct(self):
		return tCode.getOutOfPlaneDistTwoPoints(self.posA, self.posB, self.planeEqnA)

	def testCase_zNormalVector(self):
		expDist = 2
		actDist = self._runTestFunct()
		self.assertAlmostEqual(expDist, actDist)

	def testCase_yNormalVector(self):
		self.planeEqnA = tCode.ThreeDimPlaneEquation(0,1,0,0)
		expDist = 3
		actDist = self._runTestFunct()
		self.assertAlmostEqual(expDist, actDist)

	def testCase_xNormalVector(self):
		self.planeEqnA = tCode.ThreeDimPlaneEquation(1,0,0,0)
		expDist = math.sqrt( (2**2) + (3**2) )
		actDist = self._runTestFunct()
		self.assertAlmostEqual(expDist, actDist)


