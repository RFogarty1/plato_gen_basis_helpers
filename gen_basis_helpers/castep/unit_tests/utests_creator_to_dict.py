
import unittest
import unittest.mock as mock

import gen_basis_helpers.castep.castep_creator as creatorHelp
import gen_basis_helpers.castep.castep_creator_to_dict as tCode


class TestGetDictFromCreatorFuncts(unittest.TestCase):

	def setUp(self):
		self.cutoffEnergy = 500
		self.kPts = [3,3,3]
		self.charge = 2
		self.createTestObjs()

	def createTestObjs(self):
		self.kwargDictA = {"cutoffEnergy":self.cutoffEnergy, "kPts":self.kPts, "charge":self.charge}
		self.creatorObjA = creatorHelp.CastepCalcObjFactoryStandard(**self.kwargDictA)

	def testGetBasicSettingsFunctA(self):
		expDict = self.kwargDictA
		actDict = tCode.getSelectedBasicInfoDictFromCreatorObj(self.creatorObjA)
		self.assertEqual(expDict,actDict)

