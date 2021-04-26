
import copy
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.traj_core as trajCoreHelp
import gen_basis_helpers.analyse_md.surf_norm_z.shift_trajs as tCode


class TestShiftToCentre_trajInterface(unittest.TestCase):

	def setUp(self):
		self.steps = [1,2]
		self.uCellObjs = [3,4]
		self.inpIndices = [5,6]
		self.targZ = 7
		self.foldAfter = False
		self.createTestObjs()

	def createTestObjs(self):
		self.stepA = trajCoreHelp.TrajStepBase(unitCell=self.uCellObjs[0], step=self.steps[0])
		self.stepB = trajCoreHelp.TrajStepBase(unitCell=self.uCellObjs[1], step=self.steps[1])
		self.trajA = trajCoreHelp.TrajectoryInMemory([self.stepA, self.stepB])

	def _runTestFunct(self):
		args = [self.trajA, self.inpIndices, self.targZ]
		kwargs = {"foldAfter":self.foldAfter}
		tCode.shiftUnitCellToCentreAverageOfZIndices_trajInterface(*args, **kwargs)

	@mock.patch("gen_basis_helpers.analyse_md.surf_norm_z.shift_trajs._shiftUnitCellToCentreAverageZOfIndices")
	def testExpCallsA(self, mockShiftFunct):
		self._runTestFunct()

		mockShiftFunct.assert_any_call( self.uCellObjs[0], self.inpIndices, self.targZ, foldAfter=self.foldAfter)
		mockShiftFunct.assert_called_with( self.uCellObjs[1], self.inpIndices, self.targZ, foldAfter=self.foldAfter)


class TestShiftToCentreAverageZOfAtomIndices_functInterface(unittest.TestCase):

	def setUp(self):
		self.inpIndices = [1]
		self.foldAfter = False
		self.targZ = 5
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = _loadStandardUnitCellA()

	def _runTestFunct(self):
		args = [self.cellA, self.inpIndices, self.targZ]
		return tCode._shiftUnitCellToCentreAverageZOfIndices(*args, foldAfter=self.foldAfter)

	def testSingleIndiceA(self):
		expCartCoords = [ [2,2,4 , "X"],
	                      [3,3,5 , "Y"],
	                      [4,4,6 , "X"],
	                      [8,2,10, "X"] ]
		self._runTestFunct()
		actCartCoords = self.cellA.cartCoords
		self._checkExpAndActCartCoordsMatch(expCartCoords, actCartCoords)

	def testTwoIndices(self):
		expCartCoords = [ [2,2,3, "X"],
		                  [3,3,4, "Y"],
		                  [4,4,5, "X"],
		                  [8,2,9, "X"] ]
		self.inpIndices = [0,3] #Average is 5
		self.targZ = 6
		self._runTestFunct()
		actCartCoords = self.cellA.cartCoords
		self._checkExpAndActCartCoordsMatch(expCartCoords, actCartCoords)

	def _checkExpAndActCartCoordsMatch(self, expCoords, actCoords):
		expCell, actCell = copy.deepcopy(self.cellA), copy.deepcopy(self.cellA)
		expCell.cartCoords, actCell.cartCoords = expCoords, actCoords
		self.assertEqual(expCell, actCell)

def _loadStandardUnitCellA():
	outCell = uCellHelp.UnitCell(lattParams=[10,10,10], lattAngles=[90,90,90])

	cartCoords = [ [2,2,2, "X"],
	               [3,3,3, "Y"],
	               [4,4,4, "X"],
	               [8,2,8, "X"] ]

	outCell.cartCoords = cartCoords
	return outCell

