
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.special_builders.water_box as tCode


class TestGetLatticeParameterForTargetDensity(unittest.TestCase):

	def setUp(self):
		self.lattParamsA = [1,1,2]
		self.lattAnglesA = [90,90,90]
		self.targDensity = 10
		self.nWater = 1
		self.varyLattParam = "c"
		self.massDictA = {"O":1,"H":0}
		self.createTestObjs()

	def createTestObjs(self):
		self.testCellA = uCellHelp.UnitCell(lattParams=self.lattParamsA, lattAngles=self.lattAnglesA)

	def runTestFunct(self):
		args = [self.testCellA, self.nWater, self.targDensity]
		kwargs = {"lattParam":self.varyLattParam, "massDict":self.massDictA}
		return tCode.findLatticeParameterToGetTargetDensityForNWater(*args, **kwargs)

	def testExpValForSimpleCell(self):
		expVal = 0.1
		actVal = self.runTestFunct()
		self.assertAlmostEqual(expVal, actVal)



