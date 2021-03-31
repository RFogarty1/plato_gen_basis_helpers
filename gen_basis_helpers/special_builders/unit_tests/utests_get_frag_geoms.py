
import copy
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.special_builders.get_frag_geoms as tCode


class TestGetPairFragNames(unittest.TestCase):

	def setUp(self):
		self.fragNames = ["fragA", "fragB", "fragC"]
		self.maxSep = None
		
	def _runTestFunct(self):
		return tCode.getPairFragNamesFromIndividualFragNames(self.fragNames, maxSep=self.maxSep)

	def testExpectedValsNoSep(self):
		expVals = ["fragA-fragB", "fragA-fragC", "fragB-fragC"]
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)

	def testExpectedValsSepOne(self):
		self.maxSep = 1
		expVals = ["fragA-fragB", "fragB-fragC"]
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)


class TestGetIsolatedFragGeomsFromFragIndices(unittest.TestCase):


	def setUp(self):
		self.fragIndices = [ [0,2], [1,3] ]
		self.createTestObjs()

	def createTestObjs(self):
		self.inpCell = _loadStandardCellA()

	def _runTestFunct(self):
		return tCode.getIsolatedFragGeomsFromFragmentIndices(self.inpCell, self.fragIndices)

	def testExpectedGeomsA(self):
		expCells = [copy.deepcopy(self.inpCell) for x in self.fragIndices]
		expCells[0].cartCoords = [  [3,3,3,"X"], [5,5,5,"Z"] ]
		expCells[1].cartCoords = [  [4,4,4,"Y"], [6,6,6,"A"] ]

		actCells = self._runTestFunct()
		self.assertEqual(expCells,actCells)


class TestGetAllPairsOfFragGeomsFromFragIndices(unittest.TestCase):

	def setUp(self):
		self.fragIndices = [ [0], [1,3], [2] ]
		self.maxSep = None
		self.createTestObjs()

	def createTestObjs(self):
		self.inpCell = _loadStandardCellA()

	def _runTestFunct(self):
		return tCode.getFragPairGeomsFromIndices(self.inpCell, self.fragIndices, maxSep=self.maxSep)

	def testExpectedGeomsA_noSep(self):
		expNumbGeoms = 3
		expCells = [copy.deepcopy(self.inpCell) for x in range(expNumbGeoms)]
		expCells[0].cartCoords = [ [3,3,3,"X"], [4,4,4,"Y"], [6,6,6,"A"] ]
		expCells[1].cartCoords = [ [3,3,3,"X"], [5,5,5,"Z"] ]
		expCells[2].cartCoords = [ [4,4,4,"Y"], [5,5,5,"Z"], [6,6,6,"A"] ]
		actCells = self._runTestFunct()
		self.assertEqual(expCells, actCells)

	def testExpectedGeomsA_maxSepOne(self):
		self.maxSep = 1
		expNumbGeoms = 2
		expCells = [copy.deepcopy(self.inpCell) for x in range(expNumbGeoms)]
		expCells[0].cartCoords = [ [3,3,3,"X"], [4,4,4,"Y"], [6,6,6,"A"] ]
		expCells[1].cartCoords = [ [4,4,4,"Y"], [5,5,5,"Z"], [6,6,6,"A"] ]
		actCells = self._runTestFunct()
		self.assertEqual(expCells, actCells)


def _loadStandardCellA():
	lattParams, lattAngles = [10,10,10], [90,90,90]
	cartCoordsA = [ [3,3,3,"X"], [4,4,4,"Y"], [5,5,5,"Z"], [6,6,6,"A"] ]
	currKwargs = {"lattParams":lattParams, "lattAngles":lattAngles}
	inpCell = uCellHelp.UnitCell(**currKwargs)
	inpCell.cartCoords = cartCoordsA
	return inpCell

