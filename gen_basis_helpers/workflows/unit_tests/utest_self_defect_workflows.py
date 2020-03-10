
import types

import unittest
import unittest.mock as mock

import gen_basis_helpers.workflows.self_point_defects as tCode

class TestSelfPointDefectWorkflow(unittest.TestCase):

	def setUp(self):
		self.eDefect = 20
		self.eBulk = 30
		self.nAtomsDefect = 10
		self.nAtomsBulk = 2

		self.eBulkPerAtom = self.eBulk / self.nAtomsBulk
		self.eDefectPerAtom = self.eDefect/ self.nAtomsDefect
		self.createTestObjs()

	def createTestObjs(self):
		#Get expected parsed files
		energyDefect = types.SimpleNamespace(electronicTotalE=self.eDefect)
		energyBulk = types.SimpleNamespace(electronicTotalE=self.eBulk)
		parsedFileDefect = types.SimpleNamespace(energies=energyDefect, numbAtoms=self.nAtomsDefect)
		parsedFileBulk = types.SimpleNamespace(energies=energyBulk, numbAtoms=self.nAtomsBulk)
		#Create mock calc objects and attach relevant parsed files
		self.bulkObjA = mock.Mock()
		self.defectObjA = mock.Mock()
		self.bulkObjA.parsedFile = parsedFileBulk
		self.defectObjA.parsedFile = parsedFileDefect

		self.testObjA = tCode.SelfPointDefectWorkflow(self.defectObjA,self.bulkObjA)

	def testWriteFilesCalledUponInitiation(self):
		self.bulkObjA.writeFile.assert_called_once_with()
		self.defectObjA.writeFile.assert_called_once_with()

	def testEPerAtomAsExpected(self):
		actBulkEnergyPerAtom = self.testObjA._bulkEPerAtom
		actDefectEnergyPerAtom = self.testObjA._defectEPerAtom
		self.assertAlmostEqual(self.eBulkPerAtom,actBulkEnergyPerAtom)
		self.assertAlmostEqual(self.eDefectPerAtom,actDefectEnergyPerAtom)

	def testRunGivesExpectedValues(self):
		expKwargDict = {"bulkEPerAtom":self.eBulkPerAtom, "defectEPerAtom":self.eDefectPerAtom, "defectE":-130}
		expOutput = types.SimpleNamespace(**expKwargDict)
		self.testObjA.run()
		actOutputAll = self.testObjA.output
		self.assertTrue( len(actOutputAll)==1 )
		self.assertEqual(expOutput, actOutputAll[0])

