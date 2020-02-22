
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


