

import unittest
import unittest.mock as mock

import gen_basis_helpers.cp2k.cp2k_creator as creatorHelp
import gen_basis_helpers.cp2k.cp2k_basis_obj as basisObjHelp

import gen_basis_helpers.cp2k.cp2k_creator_to_dict as tCode

class TestGetDictFromCreatorFuncts(unittest.TestCase):

	def setUp(self):
		self.addedMOs = 4
		self.absGridCutoff = 600
		self.relGridCutoff = 400
		self.charge = 0
		self.kPts = [12,12,1]
		self.basisObjs = [mock.Mock()]
		self.createTestObjs()

	def createTestObjs(self):
		self.kwargDictA = {"absGridCutoff":self.absGridCutoff, "addedMOs":self.addedMOs,
		                   "basisObjs":self.basisObjs,
		                   "charge": self.charge, "kPts":self.kPts} 
		self.creatorObjA = creatorHelp.CP2KCalcObjFactoryStandard(**self.kwargDictA)

	def testGetBasicSettingsFunct(self):
		expDict = {"absGridCutoff":self.absGridCutoff, "addedMOs":self.addedMOs,
		                   "charge": self.charge, "kPts":self.kPts} 
		actDict = tCode.getSelectedBasicInfoDictFromCreatorObj(self.creatorObjA)
		self.assertEqual(expDict,actDict)

	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator_to_dict.getDictFromBasisObjList")
	def testGetBasisObjs(self, mockedDictFromBasisObjs):
		expObj = mock.Mock()
		expDict = {"basisObjs":expObj}
		mockedDictFromBasisObjs.side_effect = lambda *args: expObj
		actDict = tCode._getBasisObjsInfoFromCreatorObj(self.creatorObjA)
		mockedDictFromBasisObjs.assert_called_with(self.basisObjs)
		self.assertEqual(expDict,actDict)


class TestGetBasisDictFromBasisObjsList(unittest.TestCase):

	def setUp(self):
		self.eleA, self.eleB = "Mg", "O"
		self.basisA, self.basisB = "basis_set_a", "basis_set_b"
		self.potA, self.potB = "pot_a", "pot_b"
		self.basisFileA, self.basisFileB = "basis_file_a", "basis_file_b"
		self.ghostA, self.ghostB = False, True
		self.kindA, self.kindB = "Mg", "O_other" 
		self.createTestObjs()

	def createTestObjs(self):
		currKwargs = {"element":self.eleA, "basis":self.basisA, "potential":self.potA,
		              "basisFile":self.basisFileA, "potFile":self.potA, "ghost":self.ghostA,
		              "kind":self.kindA}
		self.testObjA = basisObjHelp.CP2KBasisObjStandard(**currKwargs)
		currKwargs = {"element":self.eleB, "basis":self.basisB, "potential":self.potB,
		              "basisFile":self.basisFileB, "potFile":self.potB, "ghost":self.ghostB,
		              "kind":self.kindB}
		self.testObjB = basisObjHelp.CP2KBasisObjStandard(**currKwargs)
		self.basisObjs = [self.testObjA, self.testObjB]

	def testToAndFromDictConsistentForCaseA(self):
		basisObjDict = tCode.getDictFromBasisObjList(self.basisObjs)
		expBasisObjsList = self.basisObjs
		actBasisObjsList = tCode.getBasisObjListFromDict(basisObjDict)

		#compare actual/expected
		sortedExpList = sorted(expBasisObjsList, key=lambda x:x.basis)
		sortedActList = sorted(actBasisObjsList, key=lambda x:x.basis)
		self.assertEqual(sortedExpList, sortedActList)


