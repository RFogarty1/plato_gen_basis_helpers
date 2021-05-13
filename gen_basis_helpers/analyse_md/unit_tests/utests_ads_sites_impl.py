

import itertools as it
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.ads_sites_impl as tCode


class TestTopSite(unittest.TestCase):

	def setUp(self):
		self.inpIdx = 1
		self.runFunctInpIdx = None
		self.inpCartCoords = None
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = _createTestCellA()
		self.testObjA = tCode.TopStandard( self.inpIdx )

	def _runTestFunct(self):
		kwargs = {"inpCartCoords":self.inpCartCoords, "inpIdx":self.runFunctInpIdx}
		return self.testObjA.positionFromGeom(self.cellA, **kwargs)

	def testExpectedValsA(self):
		expPos = [4,4,6]
		actPos = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expPos, actPos)]

class TestAdsSitesBridgeImpl(unittest.TestCase):

	def setUp(self):
		self.atomIndices = [1,2]
		self.inpCartCoords = None
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = _createTestCellA()
		self.testObjA = tCode.BridgeStandard(self.atomIndices)

	def _runTestFunct(self):
		kwargs = {"inpCartCoords":self.inpCartCoords}
		return self.testObjA.positionFromGeom(self.cellA, **kwargs)

	def testExpectedVals(self):
		expPos = [4.5,4.5,6.5]
		actPos = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expPos, actPos)]

	def testExpectedValsWhenPbcsImportant(self):
		self.atomIndices = [0,1]
		self.createTestObjs()
		self.cellA = createTestCellB_pbcsMatter()
		expPos = [10,10,10]
		actPos = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expPos, actPos)]

class TestAdsSitesHollowImpl(unittest.TestCase):

	def setUp(self):
		self.atomIndices = [0,1,2]
		self.inpCartCoords = None
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = _createTestCellA()
		self.testObjA = tCode.HollowStandard(self.atomIndices)

	def _runTestFunct(self):
		return self.testObjA.positionFromGeom(self.cellA, inpCartCoords=self.inpCartCoords)

	def testExpectedVals(self):
		expPos = [4,4,6]
		actPos = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expPos, actPos)]

	@unittest.skip("")
	def testRaisesIfLenAtomIndicesWrong(self):
		self.assertTrue(False)

	def testExpectedValsWhenPbcsImportant(self):
		self.atomIndices = [0,1,2]
		self.createTestObjs()
		self.cellA = createTestCellB_pbcsMatter()
		expPos = [31/3, 29/3, 29/3]
		actPos = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expPos,actPos)]


def _createTestCellA():
	lattParams, lattAngles = [10,10,10], [90,90,90]
	cartCoords = [ [3,3,5,"X"],
	               [4,4,6,"Y"],
	               [5,5,7,"Z"] ]
	outCell = uCellHelp.UnitCell(lattParams=lattParams, lattAngles=lattAngles)
	outCell.cartCoords = cartCoords
	return outCell


def createTestCellB_pbcsMatter():
	lattParams, lattAngles = [10,10,10], [90,90,90]
	cartCoords = [ [9,9,9,"X"],
	               [1,1,1,"Y"],
	               [1,9,9,"Z"] ]
	outCell = uCellHelp.UnitCell(lattParams=lattParams, lattAngles=lattAngles)
	outCell.cartCoords = cartCoords
	return outCell

