

import unittest
import unittest.mock as mock

import gen_basis_helpers.job_helpers.cp2k.solid_adsorption_energy_help as tCode

class TestSolidAdsorptionEnergyInputCreator(unittest.TestCase):

	def setUp(self):
		self.absGridCutoff = 2000
		self.addedMOsBulk = 10
		self.basisObjs = [mock.Mock(), mock.Mock()]
		self.cp2kMethodStr = "fake_method_str"
		self.relGridCutoff = 4000
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"absGridCutoff":self.absGridCutoff, "addedMOsBulk":self.addedMOsBulk,
		             "basisObjs":self.basisObjs, "cp2kMethodStr":self.cp2kMethodStr,
		             "relGridCutoff":self.relGridCutoff}
		self.testObjA = tCode.SolidAdsorptionEnergyStandardInputCreator(**kwargDict)

	@mock.patch("gen_basis_helpers.job_helpers.cp2k.solid_adsorption_energy_help.solidAdsorb.CodeSpecificStandardInputCreatorTemplate._modifyBulkCreatorWithSharedOptions")
	def testExpectedBulkModifications(self, mockedSuperMethod):
		testCreator = mock.Mock()
		self.testObjA._modifyBulkCreatorWithSharedOptions(testCreator)

		expAttrDict = {"addedMOs":self.addedMOsBulk}

		for key in expAttrDict.keys():
			self.assertEqual( expAttrDict[key], getattr(testCreator,key) )


	def testExpectedSharedOptionsPassedToBulkCreator(self):
		expAttrDict = {"methodStr":self.cp2kMethodStr, "absGridCutoff": self.absGridCutoff,
		               "basisObjs":self.basisObjs, "relGridCutoff":self.relGridCutoff}

		outFactory = self.testObjA._getBaseCreator()
		for key in expAttrDict.keys():
			self.assertEqual( expAttrDict[key], getattr(outFactory, key) )


