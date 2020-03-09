

import itertools as it
import unittest
import unittest.mock as mock


import gen_basis_helpers.cp2k.basis_register as tCode
import gen_basis_helpers.cp2k.private.get_default_basis_strs as defBasisStrs


class TestRegisteringBasisStrs(unittest.TestCase):

	def setUp(self):
		pass

	def testDefaultStrsPresentOnImport(self):
		defStrs = set(defBasisStrs.defaultBasisStrsToCreators.keys())
		for x in defStrs:
			self.assertTrue( x in tCode.getRegisteredCP2KObjCreatorStrings() )

	def testCreatorLedToExpectedObj(self):
		expObj = defBasisStrs._createStandardBasisObjForUnitTests()
		expObj.element = "Mg"
		actObj = tCode.createCP2KBasisObjFromEleAndBasisStr("Mg", "utests-basis")
		self.assertEqual(expObj,actObj)

	@mock.patch("gen_basis_helpers.cp2k.basis_register.createCP2KBasisObjFromEleAndBasisStr")
	def testCreateFromStrDictMakesExpectedCalls(self, mockedBasisObjGetter):
		testEleKeys =   ["Mg",      "H"]
		testBasisKeys = ["basis_a","basis_b"]
		testDict = {k:v for k,v in it.zip_longest(testEleKeys, testBasisKeys)}
		tCode.createCP2KBasisObjsFromStrDict(testDict)
		for k,v in testDict.items():
			mockedBasisObjGetter.assert_any_call(k,v)



