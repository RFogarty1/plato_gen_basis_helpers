

import unittest
import unittest.mock as mock

import gen_basis_helpers.shared.register_key_decorator as tCode


class TestKeyRegisterDecorator(unittest.TestCase):

	def setUp(self):
		self.forceKeysToCase = "lower"
		self.testRegister = dict()
		self.createTestObjs()

	def createTestObjs(self):
		decoKwargs = {"forceKeysToCase":self.forceKeysToCase}
		self.testDecoFunct = tCode.RegisterKeyValDecorator(self.testRegister,**decoKwargs)

	def testFirstValueRegistersWithoutProblems(self):
		testKey, testVal = "test_key", "test_val"
		self.testDecoFunct(testKey, testVal)
		self.assertEqual( testVal, self.testRegister[testKey] )

	def testLowerCaseKeysEnforcedIfRequested(self):
		testKey, testVal = "test_KEY".upper(), "test_val"
		self.testDecoFunct(testKey,testVal)
		self.assertEqual( testVal, self.testRegister[testKey.lower()] )

	#This test shows that ANY attribute can be temporarily set at runtime due to the way
	#this features implemented
	def testKeyEnforcementCanBeOverwrittenAtDecoTime(self):
		testKey, testVal = "test_KEY", "test_val"
		testForceKeysToCase = "upper"
		expKey = testKey.upper()
		self.assertNotEqual(testForceKeysToCase, self.testDecoFunct.forceKeysToCase) 
		self.testDecoFunct(testKey, testVal, forceKeysToCase=testForceKeysToCase)
		self.assertEqual(testVal, self.testRegister[expKey])

	def testRaisesWhenAddingSameKeyTwice(self):
		testKeyA, testKeyB = "fake_key".lower(), "fake_key".upper()
		testVal = "testVal"
		self.forceKeyToCase = "upper" #Means both keys will be effectively the same
		self.createTestObjs()
		self.testDecoFunct(testKeyA,testVal)
		with self.assertRaises(AssertionError):
			self.testDecoFunct(testKeyB, testVal)




	

