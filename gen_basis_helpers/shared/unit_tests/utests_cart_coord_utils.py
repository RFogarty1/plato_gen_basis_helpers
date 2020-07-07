
import math
import itertools as it

import unittest
import unittest.mock as mock

import gen_basis_helpers.shared.cart_coord_utils as tCode

#Separate lower-level function to get the index
class TestGetNearestPointFunctions(unittest.TestCase):

	def setUp(self):
		self.inpPoint = [1,2,3]
		self.otherPoints = [ [0,0,0], [12,12,12], [1,3,4], [1,3,5] ]

	#Note the last point is closer than the first, which would lead to a wrong answer in original implementation if i forgot to update minDistance while looping over co-ordinates
	def testExpectedMathesActualForSimpleCaseA(self):
		expCoords = self.otherPoints[2]
		actCoords = tCode.getNearestCoordToInputPoint(self.inpPoint, self.otherPoints)
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expCoords, actCoords)]

	def testGetNearestDistanceToInpPoint(self):
		expDist = math.sqrt(2)
		actDist = tCode.getNearestDistanceToInputPoint(self.inpPoint,self.otherPoints)
		self.assertAlmostEqual(expDist,actDist)
