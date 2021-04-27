

import unittest
import unittest.mock

import numpy as np

import plato_pylib.shared.ucell_class as uCellHelp

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
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

	def _runTestFunct(self):
		kwargs = {"indicesA":self.indicesA, "indicesB":self.indicesB}
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

