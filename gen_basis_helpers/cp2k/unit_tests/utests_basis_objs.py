
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


class TestGetGhostVersionsOfInputBasisObjs(unittest.TestCase):

	def setUp(self):
		self.eleA, self.basisA = "eleA", "basisA"
		self.potA, self.basFileA, self.potFileA = "potA", "basFileA", "potFileA"
		self.kindA = None
		self.ghostA = False
		self.createTestObjs()

	def createTestObjs(self):
		kwargDictA = {"element":self.eleA, "basis":self.basisA, "potential":self.potA,
		              "basisFile": self.basFileA, "potFile":self.potFileA, "ghost":self.ghostA, "kind":self.kindA}
		self.basisObjA = tCode.CP2KBasisObjStandard(**kwargDictA)

	def testGhostVersionFromNonGhostA(self):
		startObj = copy.deepcopy(self.basisObjA)
		self.eleA += "_ghost"
		self.kindA = str(self.eleA)
		self.ghostA = True
		self.createTestObjs()
		expObj = self.basisObjA
		actObj = tCode.getStandardGhostVersionOfBasisObj(startObj)
		self.assertEqual(expObj,actObj)

	def testNoChangesWhenGhostObjGiven(self):
		self.ghostA = True
		self.createTestObjs()
		expObj = self.basisObjA
		actObj = tCode.getStandardGhostVersionOfBasisObj(self.basisObjA)
		self.assertEqual(expObj, actObj)
	
	def testInputListOfGhostBasisObjsReturnedUnchanged(self):
		self.ghostA = True
		self.createTestObjs()
		basisIterA = [self.basisObjA]
		expOutput = basisIterA #Since its already a "ghost" type, we expect it unchanged
		actOutput = tCode.getBasisObjsWithGhostVersionsIncluded(basisIterA)
		self.assertEqual(expOutput, actOutput)

	def testInputListOfBasisObjsGivesExpectedOutputA(self):
		basisIterA = [self.basisObjA]
		expExtra = copy.deepcopy(self.basisObjA)
		expExtra.element = expExtra.element + "_ghost"
		expExtra.kind = expExtra.kind + "_ghost"
		expExtra.ghost = True
		expOutput = [self.basisObjA, expExtra]
		actOutput = tCode.getBasisObjsWithGhostVersionsIncluded(basisIterA)
		self.assertEqual(expOutput, actOutput)




