
import itertools as it
import math
import unittest
import unittest.mock as mock

import numpy as np

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.shared.plane_equations as planeEqnHelp
import gen_basis_helpers.analyse_md.calc_dists as tCode


class TestCalcNearestDistance(unittest.TestCase):

	def setUp(self):
		self.lattParams = [10,10,10]
		self.lattAngles = [90,90,90]
		self.coords = [ [5,5,5,"X"],
		                [5,5,6,"Y"],
		                [5,5,8,"Y"],
		                [5,5,1,"Z"] ]
		self.indicesA = None 
		self.indicesB = None
		self.minDist = 1e-2
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

	def _runTestFunct(self):
		args = [self.cellA]
		kwargs = {"indicesA":self.indicesA, "indicesB":self.indicesB, "minDist":self.minDist}
		return tCode.calcNearestDistanceBetweenTwoSetsOfIndices(*args, **kwargs)

	def testForSimpleCaseA(self):
		expVal = 1
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal,actVal)

	def testForTwoDiffSetsOfIndices(self):
		self.indicesA = [0,1]
		self.indicesB = [2,3]
		self.createTestObjs()

		expVal = 2
		actVal = self._runTestFunct()		
		self.assertAlmostEqual(expVal, actVal)

	def testCaseA_zeroExpected(self):
		self.minDist = -1e-2
		expVal = 0
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal)

class TestCalcDistMatrix(unittest.TestCase):

	def setUp(self):
		self.lattParams = [10,10,10]
		self.lattAngles = [90,90,90]
		self.coords = [ [5,5,5,"X"],
		                [5,5,6,"Y"] ]
		self.indicesA = None
		self.indicesB = None
		self.sparseMatrix = False
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

	def _runTestFunct(self):
		kwargs = {"indicesA":self.indicesA, "indicesB":self.indicesB, "sparseMatrix":self.sparseMatrix}
		return tCode.calcDistanceMatrixForCell_minImageConv(self.cellA, **kwargs)

	def testExpVal_pbcsIrrelevant(self):
		expArray = np.array( ([0, 1], [1,0]) )
		actArray = self._runTestFunct()
		self.assertTrue( np.allclose(expArray, actArray) )

	def testExpVal_pbcsMatter(self):
		self.coords[0][-2] = 1
		self.coords[1][-2] = 9
		self.createTestObjs()
		expArray = np.array( ([0,2],[2,0]) )
		actArray = self._runTestFunct()
		self.assertTrue( np.allclose(expArray,actArray) )

	def testWithOneIndicesListPassed(self):
		self.coords = [ [5,5,5,"X"],
		                [5,5,6,"Y"],
		                [5,5,7,"Z"] ]
		self.indicesA = [0,2]
		self.createTestObjs()
		expArray = np.array( ([0,2],[2,0]) )
		actArray = self._runTestFunct()
		self.assertTrue( np.allclose(expArray,actArray) )

	def testWithTwoIndicesListsPassed(self):
		self.coords = [ [5,5,5,"X"],
		                [5,5,6,"Y"],
		                [5,5,7,"Z"] ]
		self.indicesA, self.indicesB = [0,2], [0,1]
		self.createTestObjs()
		expArray = np.array( ([0,1], [2,1]) )
		actArray = self._runTestFunct()

		self.assertTrue( np.allclose(expArray,actArray) )

	def testWithTwoIndicesPassed_sparseMatrixTrueA(self):
		self.coords = [ [5,5,5,"X"],
		                [5,5,6,"Y"],
		                [5,5,7,"Z"] ]
		self.sparseMatrix = True
		self.indicesA, self.indicesB = [0,2], [0,1]
		self.createTestObjs()

		expMatrix = np.array( [ [0,1,2], [1,np.nan,1], [2,1,np.nan] ] )
		actMatrix = self._runTestFunct()

		self.assertTrue( np.allclose(expMatrix,actMatrix, equal_nan=True) )

	def testWithTwoIndicesPassed_sparseMatrixTrueB(self):
		""" This one has more np.nan than the other; with no overlapping indices """
		self.coords = [ [5,5,5,"X"],
		                [5,5,6,"Y"],
		                [5,5,7,"Z"] ]
		self.indicesA, self.indicesB = [0], [2]
		self.sparseMatrix = True
		self.createTestObjs()

		expMatrix = np.array( [ [np.nan, np.nan, 2],
		                        [np.nan, np.nan, np.nan],
		                        [2, np.nan, np.nan] ] )
		actMatrix = self._runTestFunct()

		self.assertTrue( np.allclose(expMatrix,actMatrix, equal_nan=True) )

	def testForEmptyIdxLists_sparseMatrixTrue(self):
		expMatrix = np.array( [ [np.nan, np.nan],
		                        [np.nan, np.nan] ] )
		self.indicesA, self.indicesB = list(), list()
		self.sparseMatrix = True
		self.createTestObjs()

		actMatrix = self._runTestFunct()

		self.assertTrue( np.allclose(expMatrix,actMatrix, equal_nan=True) )


class TestCalcSingleDistance(unittest.TestCase):

	def setUp(self):
		self.lattParams = [10,10,10]
		self.lattAngles = [90,90,90]
		self.coords = [ [5,5,5,"X"],
		                [5,5,6,"Y"] ]
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

	def _runFunct(self):
		return tCode.calcSingleDistBetweenCoords_minImageConv(self.cellA, self.coords[0], self.coords[1])

	def testExpVal_pbcsNoIssue(self):
		expDist = 1
		actDist = self._runFunct()
		self.assertAlmostEqual(expDist, actDist)

	def testExpVal_pbcsImportant(self):
		self.coords[0][-2] = 9
		self.coords[1][-2] = 1
		self.createTestObjs()
		expDist = 2
		actDist = self._runFunct()
		self.assertAlmostEqual(expDist, actDist, places=6)

	def testExpVal_pbcsImportant_oneAtomOutsideBox(self):
		self.coords[0][-2] = 13
		self.coords[1][-2] = 1
		self.createTestObjs()
		expDist = 2
		actDist = self._runFunct()
		self.assertAlmostEqual(expDist, actDist, places=6)


class TestCalcSingleAngle(unittest.TestCase):

	def setUp(self):
		self.lattParams = [10,10,10]
		self.lattAngles = [90,90,90]
		self.coordA = [5,4,5,"O"]
		self.coordB = [5,5,6,"H"]
		self.coordC = [5,5,4,"H"]

		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = [self.coordA, self.coordB, self.coordC]

	def _runTestFunct(self):
		return tCode.calcSingleAngleBetweenCoords_minImageConv(self.cellA, self.coordA, self.coordB, self.coordC)

	def testExpVal_pbcsNoIssue(self):
		expVal = 45
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal,actVal)

	def testExpVal_oneAtomOutsideBox(self):
		self.coordB[-2] += self.lattParams[-1]
		self.createTestObjs()
		expVal = 45
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal, places=5)


class TestGetNearestImageNebCoords(unittest.TestCase):

	def setUp(self):
		self.lattParams = [10,10,10]
		self.lattAngles = [90,90,90]
		self.coordA = [7,7,7]
		self.coordB = [7,7,8]
		self.cartToFractMatrix = None
		self.fractToCartMatrix = None
		self.createTestObjs()

	def createTestObjs(self):
		self.testCell = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		
	def _runTestFunct(self):
		kwargs = {"cartToFractMatrix":self.cartToFractMatrix, "fractToCartMatrix":self.fractToCartMatrix}
		return tCode.getNearestImageNebCoordsBasic(self.testCell, self.coordA, self.coordB) 

	def _checkExpAndActCoordsEqual(self, exp, act):
		self.assertEqual( len(exp), len(act) )
		for e,a in zip(exp,act):
			self.assertAlmostEqual(e,a)

	def testSimpleCaseWherePBCsIrrelevant(self):
		expCoord = self.coordB
		actCoord = self._runTestFunct()
		self._checkExpAndActCoordsEqual(expCoord,actCoord)

	def testWhenPBCsMatter(self):
		self.coordB[-1] += self.lattParams[-1]*-2
		self.createTestObjs()
		expCoord = [7,7,8]
		actCoord = self._runTestFunct()
		self._checkExpAndActCoordsEqual(expCoord,actCoord)

	def testWhenPBCsMatter_bothInCell(self):
		self.coordA = [7,7,9]
		self.coordB = [7,7,1]
		expCoord = [7,7,11]
		actCoord = self._runTestFunct()
		self._checkExpAndActCoordsEqual(expCoord, actCoord)

	#This failed because the tolerance for folding back was ridic low (1% of lattice paramter value)
	def testForBugFoundInProduction(self):
		lattVects = [ [19.25998959324, 0, 0],
		              [-9.62999479662, 16.6796402643698,0],
		              [0, 0, 43.41297727459] ]
		self.testCell = uCellHelp.UnitCell.fromLattVects(lattVects)

		self.coordA = [0.2462955954, -0.0845611422, 5.8293155133, "X"]
		self.coordB = [-1.65514675872001, -2.85419319106976, 27.8617855356, "Y"]

		expCoord = [-1.65514675872001, -2.85419319106976, -15.55119173899]
		actCoord = self._runTestFunct()
		self._checkExpAndActCoordsEqual(expCoord, actCoord)


class TestGetNearestImageNebCoordsMatrix(unittest.TestCase):

	def setUp(self):
		self.lattAngles, self.lattParams = [90,90,90], [10,10,10]
		self.coordA, self.coordB = [7,7,7,"X"], [7,7,18,"Y"]
		self.indicesA = [0,1]
		self.indicesB = [x for x in self.indicesA]
		self.createTestObjs()

	def createTestObjs(self):
		self.testCell = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.testCell.cartCoords = [ self.coordA, self.coordB ]

	def _runTestFunct(self):
		args = [self.testCell]
		kwargs = {"indicesA": self.indicesA, "indicesB":self.indicesB}
		return tCode.getNearestImageNebCoordsMatrixBasic(*args, **kwargs)

	def testFullMatrix_pbcsImportant(self):
		expMatrix = [ [[self.coordA , self.coordA ], [self.coordA, [7,7,8,"Y"] ]   ],
		              [[self.coordB, [7,7,17,"X"] ], [self.coordB, self.coordB]]   ]
		actMatrix = self._runTestFunct()

		for rIdx in range(len(expMatrix)):
			for cIdx in range(len(expMatrix[0])):
				[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expMatrix[rIdx][cIdx][0][:3],actMatrix[rIdx][cIdx][0][:3])]
				[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expMatrix[rIdx][cIdx][1][:3],actMatrix[rIdx][cIdx][1][:3])]
				self.assertEqual( expMatrix[rIdx][cIdx][0][-1], actMatrix[rIdx][cIdx][0][-1] )
				self.assertEqual( expMatrix[rIdx][cIdx][1][-1], actMatrix[rIdx][cIdx][1][-1] )


class TestCalcHozDistMatrix(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.cartCoordsA = [  [9,9,9, "X"],
		                      [1,1,9, "Y"],
		                      [1,2,1, "Z"] ] #Was 7, but setting to 1 should lead to same dists anyway
		self.indicesA = None
		self.indicesB = None
		self.minTotInterPlaneDist = 1e-7
		self.sparseMatrix = False
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoordsA

	def _runTestFunct(self):
		args = [self.cellA]
		kwargs = {"indicesA":self.indicesA, "indicesB":self.indicesB,
		          "minTotInterPlaneDist": self.minTotInterPlaneDist, "sparseMatrix":self.sparseMatrix}
		return tCode.calcHozDistMatrixForCell_minImageConv(*args,**kwargs)

	def _loadExpXyXzYzDists(self):
		xyDist = math.sqrt( 2**2 + 2**2 )
		xzDist = math.sqrt( 2**2 + 3**2 )
		yzDist = 1

		return xyDist, xzDist, yzDist

	def testExpValsA(self):
		xyDist = math.sqrt( 2**2 + 2**2 )
		xzDist = math.sqrt( 2**2 + 3**2 )
		yzDist = 1
		expMatrix = [ [0     , xyDist, xzDist],
		              [xyDist, 0     , yzDist],
		              [xzDist, yzDist, 0     ] ]

		actMatrix = self._runTestFunct()
		self.assertTrue( np.allclose( np.array(expMatrix), np.array(actMatrix) ) )
	
	def testExpVals_asymIndicesSpecified(self):
		self.indicesA = [0,2]
		self.indicesB = [0,1,2]

		xyDist, xzDist, yzDist = self._loadExpXyXzYzDists()
		expMatrix = [ [0     , xyDist, xzDist],
		              [xzDist, yzDist, 0     ] ]
		actMatrix = self._runTestFunct()
		self.assertTrue( np.allclose( np.array(expMatrix), np.array(actMatrix) ) )

	def testExpVals_reversedAsymIndicesSpecified(self):
		self.indicesA = [2,0]
		self.indicesB = [2,1,0]

		xyDist, xzDist, yzDist = self._loadExpXyXzYzDists()

		expMatrix = [ [0     , yzDist, xzDist],
		              [xzDist, xyDist, 0     ] ]

		actMatrix = self._runTestFunct()
		self.assertTrue( np.allclose( np.array(expMatrix), np.array(actMatrix) ) )

	def testExpVals_onlyIndicesASet(self):
		self.indicesA = [2,0]
		self.indicesB = None

		xyDist, xzDist, yzDist = self._loadExpXyXzYzDists()
		expMatrix = [ [0     , xzDist],
		              [xzDist, 0     ] ]
		actMatrix = self._runTestFunct()
		self.assertTrue( np.allclose( np.array(expMatrix), np.array(actMatrix) ) )

	def testExpVals_onlyIndicesA_sparseMatrix(self):
		self.indicesA, self.indicesB = [2,0], [0,1,2]
		self.sparseMatrix = True

		xyDist, xzDist, yzDist = self._loadExpXyXzYzDists()
		expMatrix = [ [0     , xyDist, xzDist],
		              [xyDist, np.nan, yzDist],
		              [xzDist, yzDist, 0     ] ]
		actMatrix = self._runTestFunct()
		self.assertTrue( np.allclose( np.array(expMatrix), np.array(actMatrix), equal_nan=True ) )

	@mock.patch("gen_basis_helpers.analyse_md.calc_dists.calcDistanceMatrixForCell_minImageConv")
	@mock.patch("gen_basis_helpers.analyse_md.calc_dists.getInterSurfPlaneSeparationTwoPositions")
	def testExpectedInterPlaneDistSlightlyLargerThanTotDist(self, mockCalcInterPlaneDist, mockGetDistMatrix):
		self.cartCoordsA = [ [5,5,5,"X"] ]
		self.createTestObjs()

		interPlaneDist, totalDist = 1e-9, 1e-10
		mockCalcInterPlaneDist.side_effect = lambda *args,**kwargs:  interPlaneDist
		mockGetDistMatrix.side_effect = lambda *args, **kwargs: [[totalDist]]


		expDistMatrix = [ [0] ]
		actDistMatrix = self._runTestFunct()
		self.assertEqual(expDistMatrix, actDistMatrix)


class TestGetInterSurfPlaneSeparationTwoPositions(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.posA = [5,5,9]
		self.posB = [3,3,7]

		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)

	def _runTestFunct(self):
		return tCode.getInterSurfPlaneSeparationTwoPositions(self.posA, self.posB, self.cellA)

	def testCase_pbcsDontMatter(self):
		expVal = 2
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal)

	def testCase_pbcsMatter(self):
		self.posB = [3,3,3]
		expVal = 4
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal)


class TestCalcDistanceFromSurfPlaneForCell(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.cartCoords = [ [2,2,2,"X"],
		                    [3,3,6,"Y"],
		                    [4,4,9,"Z"] ]

		self.indices = [0,1,2]
		self.planeEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,0)
		self.sparseMatrix = False
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

	def _runTestFunct(self):
		return tCode.calcDistancesFromSurfPlaneForCell(self.cellA, self.indices, self.planeEqn, sparseMatrix=self.sparseMatrix)

	def testExpectedDistsFromSurfacePlaneAtBottomOfCell(self):
		expDists = [2,4,1]
		actDists = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expDists,actDists)]

	def testExpectedWhenDefaultArgsUsed(self):
		self.planeEqn = None
		self.indices = None
		expDists = [2,4,1]
		actDists = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expDists,actDists)]

	def testExpectedDistsFromSurfPlaneAtCellMidpoint(self):
		self.planeEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,5)
		self.indices = [0,1,2]
		expDists = [3,1,4]
		actDists = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expDists, actDists)]

	def testExpectedVals_sparseMatrix(self):
		self.planeEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,5)
		self.indices = [0,2]
		self.sparseMatrix = True
		expDists = [3,np.nan,4]
		actDists = self._runTestFunct()
		self.assertTrue( np.allclose( np.array(expDists), np.array(actDists), equal_nan=True ) )


class TestCalcIterOfAnglesForInpIndices(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.cartCoords = [ [1,1,7,"X"],
		                    [1,1,9,"X"],
		                    [1,3,1,"Y"], #135 degrees once PBCs taken into acount
		                    [1,5,9,"Y"] ] #90 degrees with 0,1,this

		self.angleIndices = [ [0,1,2], [0,1,3] ]
		self.degrees = True
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

	def _runTestFunct(self):
		return tCode.getInterAtomicAnglesForInpGeom(self.cellA, self.angleIndices, degrees=self.degrees)

	def testExpected_orthogCell(self):
		expAngles = [180-45,90]
		actAngles = self._runTestFunct()
		[self.assertAlmostEqual(e,a, places=5) for e,a in it.zip_longest(expAngles, actAngles)]

	def testExpected_radians(self):
		self.degrees = False
		expAngles = [math.radians(x) for x in [180-45,90]]
		actAngles = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expAngles,actAngles)]

	def testExpected_noIndicesPassed(self):
		""" More convenient to just return no angles than to raise an error """
		self.angleIndices = list()
		expAngles = list()
		actAngles = self._runTestFunct()
		self.assertEqual(expAngles,actAngles)


#Same geom as "TestCalcNearestImageVectorMatrix"; meaning easy to steal the expected reults from there
class TestCalcNearestImageVectorsForIdxPairs(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.cartCoords = [ [1,1,7,"A"],
		                    [1,1,9,"B"],
		                    [1,3,1,"C"] ]

		self.idxPairs = [  [1,2], [0,1], [0,2], [2,1], [1,1] ]
		self.sparseMatrix = False

		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

	def _runTestFunct(self):
		return tCode.getNearestImageVectorsForIdxPairs(self.cellA, self.idxPairs, sparseMatrix=self.sparseMatrix)

	def testExpectedCaseA_listReturned(self):
		aToB, aToC, bToC = [0,0,2] , [0,2,4]  , [0,2,2]
		bToA, cToA, cToB = [0,0,-2], [0,-2,-4], [0,-2,-2]

		expVals = [ bToC, aToB, aToC, cToB, [0,0,0] ] 
		actVals = self._runTestFunct()

		for exp,act in it.zip_longest(expVals,actVals):
			self.assertTrue( np.allclose( np.array(exp), np.array(act) ) )

	def testExpectedCaseA_sparseMatrixReturned(self):
		self.sparseMatrix=True

		expMatrix = np.empty( (3,3,3) )
		expMatrix[:] = np.nan
		expMatrix[tuple([1,2])], expMatrix[tuple([2,1])] = [0,2,2], [0,-2,-2]
		expMatrix[tuple([0,1])], expMatrix[tuple([1,0])] = [0,0,2], [0,0,-2]
		expMatrix[tuple([0,2])], expMatrix[tuple([2,0])] = [0,2,4], [0,-2,-4]
		expMatrix[tuple([1,1])] = [0,0,0]

		actMatrix = self._runTestFunct()

		self.assertTrue( np.allclose(expMatrix,actMatrix,equal_nan=True) )


class TestCalcNearestImageVectorMatrix(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.cartCoords = [ [1,1,7,"A"],
		                    [1,1,9,"B"],
		                    [1,3,1,"C"] ]

		self.indicesA = None
		self.indicesB = None
		self.sparseMatrix = False

		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

	def _runTestFunct(self):
		kwargs = {"indicesA":self.indicesA, "indicesB":self.indicesB, "sparseMatrix":self.sparseMatrix}
		return tCode.getNearestImageVectorMatrixBasic(self.cellA, **kwargs)

	def _loadStandardExpectedAllIndices(self):
		zeroVector = [0,0,0]
		aToB, aToC, bToC = [0,0,2] , [0,2,4]  , [0,2,2]
		bToA, cToA, cToB = [0,0,-2], [0,-2,-4], [0,-2,-2]
		outMatrix = [ [zeroVector, aToB      , aToC],
		              [bToA      , zeroVector, bToC],
		              [cToA      , cToB      , zeroVector] ]
		return outMatrix

	def testExpectedCaseA_allIndices(self):
		expMatrix = self._loadStandardExpectedAllIndices()
		actMatrix = self._runTestFunct()
		self.assertTrue( np.allclose( np.array(expMatrix), np.array(actMatrix) ) )

	def testExpectedIndicesABPassed(self):
		self.indicesA = [0,1]
		self.indicesB = [2]
		fullMatrix = self._loadStandardExpectedAllIndices()
		expMatrix = [ [fullMatrix[0][2]],
		              [fullMatrix[1][2]] ] 
		actMatrix = self._runTestFunct()
		self.assertTrue( np.allclose(np.array(expMatrix), np.array(actMatrix)) )

	def testExpected_indicesAB_sparseMatrix(self):
		#Setup options
		self.indicesA = [0,1]
		self.indicesB = [2]
		self.sparseMatrix = True

		#Load expected
		fullMatrix = self._loadStandardExpectedAllIndices()
		expMatrix = np.empty( (3,3,3) )
		expMatrix[:] = np.nan

		expMatrix[ tuple([0,2]) ] = fullMatrix[0][2]
		expMatrix[ tuple([1,2]) ] = fullMatrix[1][2]
		expMatrix[ tuple([2,0]) ] = fullMatrix[2][0]
		expMatrix[ tuple([2,1]) ] = fullMatrix[2][1]

		#check exp vs actual
		actMatrix = self._runTestFunct()

		self.assertTrue( np.allclose(np.array(expMatrix), np.array(actMatrix), equal_nan=True) )

class TestCalcAnglesWithArbitraryVector(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.diatomA = [ [0,0,4,"X"], [0,0,5,"Y"] ]
		self.diatomB = [ [0,0,9,"X"], [0,2,1,"Y"] ] #Second should go to [0,11,2]

		#
		self.cartCoords = [ [0,0,0,"A"] ] + self.diatomA + [ [5,5,5,"B"] ] + self.diatomB  
		self.inpIndices = [ [1,2], [4,5] ]
		self.inpVector = [0,0,1]
		self.bothDirs = True #For every [1,2] also calculate the [2,1] angle

		self.createTestObjs()

	def createTestObjs(self):
		self.geomA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.geomA.cartCoords = self.cartCoords

	def _runTestFunct(self):
		currArgs = [self.geomA, self.inpIndices, self.inpVector]
		return tCode.calcSparseDiatomAngleWithArbVectorMatrix(*currArgs, bothDirs=self.bothDirs)

	def testExpOutputA(self):
		expMatrix = np.empty( (6,6) )
		expMatrix[:] = np.nan
		expMatrix[1][2], expMatrix[2][1] = 0,180 
		expMatrix[4][5], expMatrix[5][4] = 45, 135
		actMatrix = self._runTestFunct()
		self.assertTrue( np.allclose(np.array(expMatrix), np.array(actMatrix), equal_nan=True) )

	def testExpOutputB(self):
		self.bothDirs = False
		expMatrix = np.empty( (6,6) )
		expMatrix[:] = np.nan
		expMatrix[1][2], expMatrix[4][5] = 0, 45
		actMatrix = self._runTestFunct()
		self.assertTrue( np.allclose(np.array(expMatrix), np.array(actMatrix), equal_nan=True) )



