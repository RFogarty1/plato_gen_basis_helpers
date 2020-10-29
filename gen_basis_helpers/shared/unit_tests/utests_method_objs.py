

import unittest
import unittest.mock as mock

import gen_basis_helpers.shared.method_objs as tCode



class TestStandardParsedFileObj(unittest.TestCase):

	def setUp(self):
		self.energies = mock.Mock()
		self.numbAtoms = mock.Mock()
		self.unitCell = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		self.kwargDict = {"energies":self.energies, "numbAtoms":self.numbAtoms, "unitCell":self.unitCell}

	def testInitFromKwargDict(self):
		testObj = tCode.StandardParsedOutputFile(**self.kwargDict)
		for attr in self.kwargDict.keys():
			objAttr = getattr(testObj, attr)
			self.assertEqual( self.kwargDict[attr], objAttr )


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

