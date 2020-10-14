

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

class TestGetDictFromObjCreatorClass(unittest.TestCase):

	def setUp(self):
		self.retDictA = {"key_a":"val_a" , "key_b":"val_b"}
		self.retDictB = {"key_a":"val_a2", "key_c":"val_c"}
		self.mockCreatorA = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		self.getterFunctA = mock.Mock()
		self.getterFunctB = mock.Mock()
		self.getterFunctA.side_effect = lambda x: self.retDictA 
		self.getterFunctB.side_effect = lambda x: self.retDictB
		self.delFunctA = lambda x: x.pop("key_a")

	def testWithSingleGetterFunct(self):
		testObj = tCode.GetDictFromCreatorObj([self.getterFunctA])
		expOutput = self.retDictA
		actOutput = testObj( self.mockCreatorA )
		self.assertEqual(expOutput, actOutput)

	def testWithTwoGetterFuncts_secondOverwriting(self):
		testObj = tCode.GetDictFromCreatorObj([self.getterFunctA, self.getterFunctB])
		expOutput = {"key_a":"val_a2", "key_b":"val_b", "key_c":"val_c"}
		actOutput = testObj( self.mockCreatorA )
		self.assertEqual(expOutput, actOutput)

	def testWithFinalFunctDeletingKeys(self):
		testObj = tCode.GetDictFromCreatorObj([self.getterFunctA], modFinalDictFuncts=[self.delFunctA])
		expOutput = {"key_b":"val_b"}
		actOutput = testObj( self.mockCreatorA )
		self.assertEqual(expOutput, actOutput)

