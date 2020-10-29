

import unittest
import unittest.mock as mock

import gen_basis_helpers.cp2k.cp2k_creator as creatorHelp
import gen_basis_helpers.cp2k.cp2k_creator_to_dict as tCode

class TestGetDictFromCreatorFuncts(unittest.TestCase):

	def setUp(self):
		self.addedMOs = 4
		self.absGridCutoff = 600
		self.relGridCutoff = 400
		self.charge = 0
		self.kPts = [12,12,1]
		self.createTestObjs()

	def createTestObjs(self):
		self.kwargDictA = {"absGridCutoff":self.absGridCutoff, "addedMOs":self.addedMOs,
		                   "charge": self.charge, "kPts":self.kPts} 
		self.creatorObjA = creatorHelp.CP2KCalcObjFactoryStandard(**self.kwargDictA)

	def testGetBasicSettingsFunct(self):
		expDict = self.kwargDictA
		actDict = tCode.getSelectedBasicInfoDictFromCreatorObj(self.creatorObjA)
		self.assertEqual(expDict,actDict)


