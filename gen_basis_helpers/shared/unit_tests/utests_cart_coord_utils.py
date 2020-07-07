
import math
import itertools as it

import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCell

import gen_basis_helpers.shared.plane_equations as planeEqnHelp
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


class TestFilterBasedOnWhetherInPlane(unittest.TestCase):

	def setUp(self):
		self.inpCoords = [ [0,1,2],
		                   [0,2,1],
		                   [0,0,0],
		                   [0,0,1] ]
		self.planeEqnA = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,1)
		self.inpPoint = [0,3,1]

	def testFiltersOutCoordsNotInZEqualOnePlane(self):
		expIndices = [1,3]
		actIndices = tCode.getFilteredIndicesForCoordsInInputPlane(self.inpCoords, self.planeEqnA)
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expIndices,actIndices)]

	def testFiltersOutCoordsInZPlane(self):
		expIndices = [0,2]
		actIndices = tCode.getFilteredIndicesForCoordsOutOfInputPlane(self.inpCoords, self.planeEqnA)
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expIndices,actIndices)]

	def testGetNearestInPlaneNeighbourIdx(self):
		expIndex = 1
		actIndex = tCode.getIdxOfNearestInPlanePointToInpPoint(self.inpPoint, self.inpCoords, self.planeEqnA)
		self.assertEqual(expIndex,actIndex)

	def testGetNearestOutOfPlaneNeighbourIdx(self):
		expIndex = 0
		actIndex = tCode.getIdxOfNearestOutOfPlanePointToInpPoint(self.inpPoint, self.inpCoords, self.planeEqnA)
		self.assertEqual(expIndex,actIndex)

	def testGetNearestInPlaneCoords(self):
		expCoords = self.inpCoords[1]
		actCoords = tCode.getCoordsOfNearestInPlanePointToInpPoint(self.inpPoint, self.inpCoords, self.planeEqnA)
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expCoords,actCoords)]

	def testGetNearestOutOfPlaneCoords(self):
		expCoords = self.inpCoords[0]
		actCoords = tCode.getCoordsOfNearestOutOfPlanePointToInpPoint(self.inpPoint, self.inpCoords, self.planeEqnA)
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expCoords, actCoords)]

	def testGetNearestInPlaneDistance(self):
		expDist = 1
		actDist = tCode.getDistanceToNearestInPlanePointToInpPoint(self.inpPoint, self.inpCoords, self.planeEqnA)
		self.assertAlmostEqual(expDist,actDist)

	def testGetNearestOutOfPlaneDistance(self):
		expDist = math.sqrt(5)
		actDist = tCode.getDistanceToNearestOutOfPlanePointToInpPoint(self.inpPoint, self.inpCoords, self.planeEqnA)
		self.assertAlmostEqual(expDist,actDist)


#Lots of this originally taken from the add_interstitials test code of the time
class TestUnitCellInterfaceFunctions(unittest.TestCase):

	def setUp(self):
		self.lattParams = [10,9,8]
		self.lattAngles = [90,90,90]
		self.fractPositions = [ [0.5,0.5,0.5], [0.51,0.5,0.5] ]
		self.atomList = ["X" for x in self.fractPositions]
		self.planeEqnA = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,1)
		self.testAtomIdx = 0
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"lattParams":self.lattParams, "lattAngles":self.lattAngles,
		             "fractCoords": self.fractPositions, "elementList":self.atomList}
		self.testCellA = uCell.UnitCell( **kwargDict )

	def testNearestNebInPlaneDistance_singleLayerWithImages(self):
		self.fractPositions[-1][-1] = self.fractPositions[0][-1]+0.1
		self.assertNotAlmostEqual(self.fractPositions[-1][-1],self.fractPositions[0][-1])
		expDist = self.lattParams[1] #Nearest neighbour has to be its own image
		actDist = tCode.getNearestInPlaneDistanceGivenInpCellAndAtomIdx(self.testCellA, self.testAtomIdx,self.planeEqnA)
		self.assertAlmostEqual(expDist,actDist)

	def testNearestNebInPlaneDistanceNoImages(self):
		self.fractPositions = [ [0.01,0.5,0.5], [0.99,0.5,0.5] ]
		self.createTestObjs()
		expDist = 0.98*self.lattParams[0]
		actDist = tCode.getNearestInPlaneDistanceGivenInpCellAndAtomIdx(self.testCellA, self.testAtomIdx,self.planeEqnA, includeImages=False)
		self.assertAlmostEqual(expDist, actDist)

	def testNearestNebInPlane_leftImage(self):
		""" Test this works when the nearest neighbour is at -ve x (this would have failed with an old, never commited, implementation """
		self.fractPositions = [ [0.01,0.5,0.5], [0.99,0.5,0.5] ]
		self.createTestObjs()
		expDist = 0.02*self.lattParams[0]
		actDist = tCode.getNearestInPlaneDistanceGivenInpCellAndAtomIdx(self.testCellA, self.testAtomIdx,self.planeEqnA)
		self.assertAlmostEqual(expDist,actDist)



