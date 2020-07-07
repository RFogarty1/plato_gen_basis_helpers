
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


