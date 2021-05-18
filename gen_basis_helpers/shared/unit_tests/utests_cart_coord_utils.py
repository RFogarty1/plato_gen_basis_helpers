
import copy
import math
import itertools as it
import types

import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCell

import gen_basis_helpers.shared.plane_equations as planeEqnHelp
import gen_basis_helpers.shared.cart_coord_utils as tCode


class TestGetDistanceOfAtomsFromPlaneEquation(unittest.TestCase):

	def setUp(self):
		self.lattParams = [10,10,10]
		self.lattAngles = [90,90,90]
		self.coordsA = [ [2,2,3,"X"],
		                 [3,3,3,"X"],
		                 [5,5,5,"Y"],
		                 [5,5,9,"Z"] ]
		self.planeEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,0)
		self.atomIndices = [idx for idx,unused in enumerate(self.coordsA)]
		self.signed = False
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCell.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coordsA

	def _runTestFunct(self):
		args = [self.cellA, self.planeEqn, self.atomIndices]
		kwargs = {"signed":self.signed}
		return tCode.getDistancesOfAtomsFromPlaneEquation_nearestImageAware(*args, **kwargs)

	def testUnsignedDistsA_pbcsImportant(self):
		expDists = [3,3,5,1]
		actDists = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expDists, actDists)]

	def testSignedDistsA_pbcsImportant(self):
		self.signed = True
		expDists = [3,3,5,-1]
		actDists = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expDists, actDists)]

	def testWithoutUsingAllIndices(self):
		self.atomIndices = [2,3]
		expDists = [5,1]
		actDists = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expDists, actDists)]


class TestGetAverageSurfacePlaneEqnForIndices(unittest.TestCase):

	def setUp(self):
		self.cellA =_loadCellA()
		self.atomIndices = [2,3]
		self.planeEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,0)

	def _runTestFunct(self):
		return tCode.getAveragePlaneEqnForAtomIndices(self.cellA, self.atomIndices, self.planeEqn)

	def testExpectedAverageEqnA(self):
		expD = 6
		expPlaneEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,expD)
		actPlaneEqn = self._runTestFunct()
		self.assertEqual(expPlaneEqn, actPlaneEqn)

	def testForSurfPlaneEqnInterface(self):
		expD = 6
		expPlaneEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,expD)
		actPlaneEqn = tCode.getAverageSurfacePlaneEquationForAtomIndices(self.cellA, self.atomIndices)
		self.assertEqual(expPlaneEqn, actPlaneEqn)

def _loadCellA():
	lattParams, lattAngles = [10,10,10], [90,90,90]
	cartCoords =  [ [2,2,3,"X"],
		            [3,3,3,"X"],
		            [5,5,5,"Y"],
		            [5,5,7,"Z"] ]
	outCell = uCell.UnitCell(lattParams=lattParams, lattAngles=lattAngles)
	outCell.cartCoords = cartCoords
	return outCell


class TestShiftCoordsToLeaveEqualVacAboveAndBelow(unittest.TestCase):

	def setUp(self):
		self.lattParams = [10,9,8]
		self.lattAngles = [90,90,90]
		self.cartCoordsA = [ [5,5,6,"X"],
		                     [5,5,7,"X"] ]
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = uCell.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.testObjA.cartCoords = self.cartCoordsA

	def _loadExpResultForSimpleOrthogonalCell(self):
		expCoords = [ [5, 5, 3.5, "X"],
		              [5, 5, 4.5, "X"] ]
		expCell = uCell.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		expCell.cartCoords = expCoords
		return expCell

	def testExpectedForOrthogonalCell(self):
		expCell = self._loadExpResultForSimpleOrthogonalCell()
		tCode.shiftCoordsToLeaveEqualVacAboveAndBelowSurface(self.testObjA)
		actCell = self.testObjA #works in place
		self.assertEqual(expCell,actCell)

	def testExpectedForOrthogonalCell_extraVacBelow(self):
		self.cartCoordsA = [ [5,5,1,"X"],
		                     [5,5,2,"X"] ]
		self.createTestObjs()
		expCell = self._loadExpResultForSimpleOrthogonalCell()
		tCode.shiftCoordsToLeaveEqualVacAboveAndBelowSurface(self.testObjA)
		actCell = self.testObjA
		self.assertEqual(expCell,actCell)

	def testExpectedForNoVacuumCase(self):
		self.cartCoordsA = [ [5,5,0,"X"],
		                     [5,5,8,"X"] ]
		self.createTestObjs()
		expCell = copy.deepcopy(self.testObjA)
		tCode.shiftCoordsToLeaveEqualVacAboveAndBelowSurface(self.testObjA)
		actCell = self.testObjA
		self.assertEqual(expCell,actCell)

	def testValueErrorIfAtomsOutsideCell(self):
		self.cartCoordsA = [ [5,5,8,"X"],
		                     [5,5,9,"X"] ]
		self.createTestObjs()
		with self.assertRaises(ValueError):
			tCode.shiftCoordsToLeaveEqualVacAboveAndBelowSurface(self.testObjA)

	@unittest.skip("Couldnt be bothered to figure out what the answer should be")
	def testExpectedValsForNonOrthogCell(self):
		self.lattParams = [10,10,10]
		self.lattAngles = [60,90,90]
		height = 10*math.cos(math.radians(60))
		self.cartCoordsA = [ [5,5,0.9*height,"X"],
		                     [5,5,0.8*height,"X"] ]
		self.createTestObjs()
		expCoords = [ [5, 5, 0.55*height, "X"],
		              [5, 5, 0.45*height, "X"] ]
		expCell = uCell.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		expCell.cartCoords = expCoords
		tCode.shiftCoordsToLeaveEqualVacAboveAndBelowSurface(self.testObjA)
		actCell = self.testObjA
		import pdb
		pdb.set_trace()
		self.assertEqual(expCell,actCell)



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

	def testGetSurfacePlanePointingSameDirAsCVector(self):
		""" Want the surface plane-equation with the normal vector pointing along c for the cubic cell (and the same direction for other cells)"""

		expPlaneEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,0)
		actPlaneEqn = tCode.getABPlaneEqnWithNormVectorSameDirAsC(self.testCellA.lattVects)
		self.assertEqual(expPlaneEqn,actPlaneEqn)

	def testGetSurfacePlanePointingSameDirAsCVector_negativeCVector(self):
		lattVects = self.testCellA.lattVects
		lattVects[-1] = [-1*x for x in lattVects[-1]]
		expPlaneEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,-1,0)
		actPlaneEqn = tCode.getABPlaneEqnWithNormVectorSameDirAsC(lattVects)
		self.assertEqual(expPlaneEqn,actPlaneEqn)

	def testGetSurfacePlanePointingSameDirAsCVectorUCellInterface_negativeCVector(self):
		lattVects = self.testCellA.lattVects
		lattVects[-1] = [-1*x for x in lattVects[-1]]
		stubInpCell = types.SimpleNamespace(lattVects=lattVects)
		expPlaneEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,-1,0)
		actPlaneEqn = tCode.getABPlaneEqnWithNormVectorSameDirAsC_uCellInterface(stubInpCell)
		self.assertEqual(expPlaneEqn, actPlaneEqn)


class TestGetDistMatrixBetweenTwoSetsOfSeparateCoords(unittest.TestCase):

	def setUp(self):
		self.coordsA = [ [0.0,0.0,0.0] ]

		self.coordsB = [ [1.0, 0.0, 0.0],
		                 [0.0, 2.0, 0.0],
		                 [0.0, 0.0, 3.0] ]

	def testActMatchesExpected_verySimple(self):
		expMatrix = [ [1,2,3] ]
		actMatrix = tCode._getDistMatrixBetweenTwoSetsOfSeparateCoords(self.coordsA, self.coordsB)
		for rIdx,unused in enumerate(expMatrix):
			for cIdx,unused in enumerate(expMatrix[rIdx]):
				exp,act = expMatrix[rIdx][cIdx], actMatrix[rIdx][cIdx]
				self.assertAlmostEqual(exp,act)

class TestGetDistMatrixSingleCoordSet(unittest.TestCase):

	def setUp(self):
		self.coordsA = [ [0.0,0.0,0.0],
		                 [0.0,1.0,0.0] ]

	def testActMatchesExpectedA(self):
		expMatrix = [ [0.0, 1.0],
		              [1.0, 0.0] ]
		actMatrix = tCode._getDistMatrixForSetOfCoords(self.coordsA)
		for rIdx,unused in enumerate(expMatrix):
			for cIdx,unused in enumerate(actMatrix):
				exp,act = expMatrix[rIdx][cIdx], actMatrix[rIdx][cIdx]
				self.assertAlmostEqual(exp,act)

class TestGetClosestDistanceBetweenTwoElements(unittest.TestCase):

	def setUp(self):
		self.lattParamsA = [2,2,2]
		self.lattAnglesA = [90,90,90]
		self.cartCoordsA = [ [1  , 1,  1, "A"],
		                     [0.1, 1,  1, "B"],
		                     [1.9, 1,  1, "C"] ]

		self.eleA, self.eleB = "A", "B"
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"lattParams":self.lattParamsA, "lattAngles":self.lattAnglesA}
		self.testCellA = uCell.UnitCell(**kwargDict)
		self.testCellA.cartCoords = self.cartCoordsA

	def testCentralCellOnly(self):
		expDist = 0.9
		actDist = tCode.getClosestDistanceBetweenTwoElementsForInpCell(self.testCellA, self.eleA, self.eleB, inclImages=True)
		self.assertAlmostEqual(expDist,actDist)

	def testExpectedFromCentralAndImages(self):
		expDist = 0.2
		actDist = tCode.getClosestDistanceBetweenTwoElementsForInpCell(self.testCellA, "B", "C")
		self.assertAlmostEqual(expDist, actDist)

	def testCentralForSameElementCentralOnly(self):
		self.cartCoordsA.append([0.8,1,1,"A"])
		self.createTestObjs()
		expDist = 0.2
		actDist = tCode.getClosestDistanceBetweenTwoElementsForInpCell(self.testCellA, "A", "A")
		self.assertAlmostEqual(expDist,actDist)

class TestGetHeightOfCell(unittest.TestCase):

	def setUp(self):
		self.lattParams = [10,9,8]
		self.lattAngles = [90,90,90]
		self.createTestObjs()

	def createTestObjs(self):
		self.testCell = uCell.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)

	def testForOrthogonalCell(self):
		expVal = self.lattParams[-1]
		actVal = tCode.getHeightOfCell_abSurface(self.testCell)
		self.assertAlmostEqual(expVal,actVal)

	def testForAlpha60(self):
		self.lattAngles = [60,90,90]
		self.createTestObjs()
		expVal = self.lattParams[-1]*math.sin(math.radians(self.lattAngles[0]))
		actVal = tCode.getHeightOfCell_abSurface(self.testCell)
		self.assertAlmostEqual(expVal, actVal)


