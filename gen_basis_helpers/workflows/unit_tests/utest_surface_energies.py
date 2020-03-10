
import types
import unittest
import unittest.mock as mock

import gen_basis_helpers.workflows.surface_energies as tCode


class TestSurfaceEnergyWorkflow(unittest.TestCase):


	def setUp(self):
		self.bulkEnergy = 4
		self.bulkNumbAtoms = 1
		self.surfEnergy = 6
		self.surfNumbAtoms = 2
		self.surfaceArea = 3
		self.surfCalcObj = mock.Mock()
		self.bulkCalcObj = mock.Mock()
		self.surfaceAreaFromUnitCell = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		energiesObjBulk = types.SimpleNamespace( electronicTotalE=self.bulkEnergy )
		energiesObjSurface = types.SimpleNamespace( electronicTotalE=self.surfEnergy )
		self.bulkCalcObj.parsedFile.energies = energiesObjBulk
		self.bulkCalcObj.parsedFile.numbAtoms = self.bulkNumbAtoms

		self.surfCalcObj.parsedFile.energies = energiesObjSurface
		self.surfCalcObj.parsedFile.numbAtoms = self.surfNumbAtoms

		self.surfaceAreaFromUnitCell.side_effect = lambda x: self.surfaceArea

		self.testObjA = tCode.SurfaceEnergyWorkflow(self.surfCalcObj, self.bulkCalcObj,
		                                            self.surfaceAreaFromUnitCell)

	def testEnergyPerAtomBulk(self):
		expValue = self.bulkEnergy / self.bulkNumbAtoms
		actValue = self.testObjA._energyPerAtomBulk
		self.assertEqual(expValue,actValue) 


	def testEnergyPerAtomSurface(self):
		expValue = self.surfEnergy / self.surfNumbAtoms
		actValue = self.testObjA._energyPerAtomSurface
		self.assertEqual(expValue,actValue)

	def testRunGivesExpectedVal(self):
		self.testObjA.run()
		expVal = -(1/3)
		actVal = self.testObjA.output[0].surfaceEnergy
		self.assertAlmostEqual(expVal,actVal)

	def testRunGivesExpectedEPerAtomVals(self):
		self.testObjA.run()
		expSurfEPerAtom, expBulkEPerAtom = self.surfEnergy/self.surfNumbAtoms, self.bulkEnergy/self.bulkNumbAtoms
		actSurfEPerAtom, actBulkEPerAtom = self.testObjA.output[0].surfEPerAtom, self.testObjA.output[0].bulkEPerAtom
		self.assertAlmostEqual(expSurfEPerAtom, actSurfEPerAtom)
		self.assertAlmostEqual(expBulkEPerAtom, actBulkEPerAtom)

