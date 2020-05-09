
import types
import unittest
import unittest.mock as mock

import gen_basis_helpers.workflows.geom_opt_workflow as tCode


#Quite a lot of duplication with the total energy workflow tests
class TestGeomOptWorkflow(unittest.TestCase):

	def setUp(self):
		self.outGeom = mock.Mock()
		self.outEnergy = 30
		self.createTestObjs()

	def createTestObjs(self):
		self.calcObjA = mock.Mock()
		energiesObj = types.SimpleNamespace(electronicTotalE=self.outEnergy)
		self.parsedFileA = types.SimpleNamespace(energies=energiesObj, unitCell=self.outGeom)
		self.calcObjA.parsedFile = self.parsedFileA
		self.testObjA = tCode.GeomOptWorkflow(self.calcObjA)

	def testWritesInpFilesUponInitialisation(self):
		self.calcObjA.writeFile.assert_called_once_with()

	def testExpectedGeomGiven(self):
		self.testObjA.run()
		expGeom = self.outGeom
		actGeom = self.testObjA.output[0].geom
		self.assertEqual(expGeom,actGeom)

	def testExpectedEnergyGiven(self):
		self.testObjA.run()
		expEnergy = self.outEnergy
		actEnergy = self.testObjA.output[0].energy
		self.assertEqual(expEnergy, actEnergy)

