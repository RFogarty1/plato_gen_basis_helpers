

import copy
import types
import unittest
import unittest.mock as mock

import gen_basis_helpers.shared.creator_resetable_kwargs as tCode

class DudConcreteCreator(tCode.CreatorWithResetableKwargsTemplate):

	registeredKwargs = set( tCode.CreatorWithResetableKwargsTemplate.registeredKwargs ) #Note we NEED to copy tihs. Though we could make a blank set in truth....
	registeredKwargs.add("attrA")
	registeredKwargs.add("attrB")

	def _createFromSelf(self):
		return types.SimpleNamespace(attrA=self.attrA, attrB=self.attrB)



class TestTemplateCreator(unittest.TestCase):

	def setUp(self):
		self.kwargDict = {"attrA":"attA", "attrB":None}
		self.createTestObj()

	def createTestObj(self):
		self.testObjA = DudConcreteCreator( **self.kwargDict )

	def testProperKwargSetByInit(self):
		for key in self.kwargDict:
			self.assertEqual( self.kwargDict[key], getattr(self.testObjA,key) )

	def testErrorOnInitFromUnregisteredKwarg(self):
		self.kwargDict["not_really_a_kwarg"] = None
		with self.assertRaises(KeyError):
			self.createTestObj()

	def testErrorOnCreateFromUnregisteredKwarg(self):
		with self.assertRaises(KeyError):
			self.testObjA.create(this_is_not="a_registeredKwarg")

	def testExpectedObjectCreatedWhenNoArgsPassed(self):
		expObj = types.SimpleNamespace(**self.kwargDict)
		actObj = self.testObjA.create()
		self.assertEqual(expObj,actObj)

	def testExpectedObjCreatedWhenPassingArgs(self):
		modKwargDict = {"attrB":"new_attr_b"}
		expDict = dict(self.kwargDict)
		expDict.update(modKwargDict)
		expObj = types.SimpleNamespace(**expDict)
		actObj = self.testObjA.create(**expDict)
		self.assertEqual(expObj, actObj)

	def testStateUnmodifiedAfterPassingAKwargToCreate(self):
		startAttrA = copy.deepcopy(self.testObjA.attrA)
		testAttrA = "hopefully_different"
		self.assertNotEqual(startAttrA,testAttrA)
		self.testObjA.create(attrA=testAttrA)
		self.assertEqual(startAttrA, self.testObjA.attrA)		


