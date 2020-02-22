
import copy
import unittest

import gen_basis_helpers.cp2k.cp2k_basis_obj as tCode

class TestStandardBasisObj(unittest.TestCase):

	def setUp(self):
		self.testAttrs = ["basis", "element", "potential", "potFile", "basisFile"]
		for attr in self.testAttrs:
			setattr(self, attr, "test_" + attr)

		self.createTestObjs()

	def createTestObjs(self):
		self.kwargDict = {k:getattr(self,k) for k in self.testAttrs}
		self.testObjA = tCode.CP2KBasisObjStandard(**self.kwargDict)

	def testAllPropsSet(self):
		for attr in self.testAttrs:
			self.assertEqual( getattr(self,attr), getattr(self.testObjA,attr) )

	#Attribute error seems to make sense
	def testArgsRequired(self):
		for attr in self.testAttrs:
			currKwargDict = dict(self.kwargDict)
			currKwargDict.pop(attr)
			with self.assertRaises(AttributeError):
				tCode.CP2KBasisObjStandard(**currKwargDict)


	def testEqualityTrueForEqualStructs(self):
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = self.testObjA
		self.assertEqual(objA, objB)

	def testEqualityFalseForUnequalStructs(self):
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = self.testObjA
		objB.element = "random-ele_str"
		self.assertNotEqual( objA.element, objB.element )
		self.assertNotEqual(objA,objB)

