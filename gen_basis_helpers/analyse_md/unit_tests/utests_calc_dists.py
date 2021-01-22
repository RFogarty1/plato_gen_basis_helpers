

import unittest
import unittest.mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.calc_dists as tCode

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


