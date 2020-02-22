

import unittest
import unittest.mock as mock

import gen_basis_helpers.cp2k.method_register as tCode
import gen_basis_helpers.cp2k.private.get_default_method_strs as defMethStrs

class TestRegisteringMethodStrs(unittest.TestCase):

	def setUp(self):
		pass

	def testDecoratorAddsStr(self):
		testStr, testMethod = "test", mock.Mock()
		self.assertFalse( testStr in tCode.getRegisteredCP2KObjCreatorStrings() )
		tCode.registerCP2KObjCreatorToMethodStr(testMethod, testStr)
		self.assertTrue( testStr in tCode.getRegisteredCP2KObjCreatorStrings() )

	def testDefaultStrsPresentOnImport(self):
		""" Test that certain default strings are present in our CP2K metod object "database" after import """
		defStrs = set(defMethStrs.defaultMethodStrsToObjCreators.keys())
		for x in defStrs:
			self.assertTrue( x in tCode.getRegisteredCP2KObjCreatorStrings() )

	def testDecoratorFunctWorked(self):
		self.assertTrue( "FakeDecoStr".lower() in tCode.getRegisteredCP2KObjCreatorStrings() )

	def testStrCantBeRegisteredTwiceByDefault(self):
		testStr, testMethod = "testStr", mock.Mock()
		self.assertFalse( testStr.lower() in tCode.getRegisteredCP2KObjCreatorStrings() )
		tCode.registerCP2KObjCreatorToMethodStr(testMethod, testStr)
		with self.assertRaises(AssertionError):
			tCode.registerCP2KObjCreatorToMethodStr(testMethod, testStr)

	def testOverwriteMethodWorks(self):
		testStr, testMethod = "testStr123", mock.Mock()
		tCode.registerCP2KObjCreatorToMethodStr(testMethod, testStr, overwrite=True)
		try:
			tCode.registerCP2KObjCreatorToMethodStr(testMethod, testStr, overwrite=True)
		except AssertionError:
			self.fail("Assertion error unexpectedly raised")


@tCode.decoRegisterCP2KObjCreatorToMethodStr("FakeDecoStr")
def fakeFunct():
	pass




