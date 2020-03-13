#!/usr/bin/python3

import copy
import types
import unittest

import gen_basis_helpers.shared.label_objs as labelObjs
import gen_basis_helpers.shared.misc_utils as tCode

#Note: I'm pretty sure this is tested indirectly elsewhere 
# (on specific classes which use it, as it was factored out of
# a single-class implementation)
class TestGetObjsWithComponentsWrapper(unittest.TestCase):

	def setUp(self):
		self.labelsA = types.SimpleNamespace( components=["testStringA"], labelNames=["testA"] )
		self.labelsB = types.SimpleNamespace( components=["testThingA"], labelNames=["testB"] )
		self.createTestObjLeaf()

	def createTestObjLeaf(self):
		self.testObjLeaf = types.SimpleNamespace( label = [self.labelsA] )	
#		decorator = tCode.getObjectsWithComponentsInstanceWrapper( isComposite=False )
		tCode.wrapInstanceWithGetObjectsWithComponentsSearch( self.testObjLeaf, isComposite=False )


	def testSearchWithPartialMatches_leafObj(self):
		compSearch = "test" #This should only partially match a component in labelsA
		expNumbOutObjs = 1
		expOutObjs = [ types.SimpleNamespace(label=[self.labelsA]) ]
		actOutObjs = self.testObjLeaf.getObjectsWithComponents( [compSearch], partialMatch=True )
		self.assertEqual( expNumbOutObjs, len(actOutObjs) )
		self.assertEqual(expOutObjs[0].label, actOutObjs[0].label)

	def testSearchWithPartialMatchesOffFails_leafObj(self):
		compSearch = "test"
		expNumbOutObjs = 0
		actOutObjs = self.testObjLeaf.getObjectsWithComponents( [compSearch], partialMatch=False )
		self.assertEqual ( expNumbOutObjs, len(actOutObjs) )



class DudClass():
	def __init__(self, **kwargs):
		for x in kwargs.keys():
			setattr(self, x, kwargs[x])

class TestStandardComponentAttr(unittest.TestCase):

	def setUp(self):
		self.leafA = types.SimpleNamespace( testComp = ["testA"] )
		self.leafB = types.SimpleNamespace( testComp = ["testB"] )
		self.leafC = types.SimpleNamespace( testComp = ["testC"] )
		self.createTestCompObjWithLeafsOnly()
		self.createTestCompObjLeafsAndComp()

	def createTestCompObjWithLeafsOnly(self):
		setattr(DudClass, "testComp", tCode.StandardComponentDescriptor("testComp"))
		self.testCompObjLeafsOnly = DudClass( objs=[self.leafA,self.leafB] )

	def createTestCompObjLeafsAndComp(self):
		self.createTestCompObjWithLeafsOnly()
		self.testCompObjLeafAndComp = DudClass( objs = [copy.deepcopy(self.testCompObjLeafsOnly), self.leafC] )

	def testAttrObjWithLeafsOnly(self):
		expResult = ["testA","testB"]
		actResult = getattr(self.testCompObjLeafsOnly,"testComp")
		self.assertEqual(expResult, actResult)

	def testAttrObjWithLeafAndComp(self):
		expResult = ["testA","testB","testC"]
		actResult = getattr(self.testCompObjLeafAndComp, "testComp")
		self.assertEqual(expResult, actResult)

class DudClassForTestingUniqueLabelEnforcement():

	@tCode.getAssertAllLabelsUniqueUponCreationClassInitializerWrapper()
	def __init__(self, labelList):
		self._labels = list(labelList)

	@property
	def label(self):
		return self._labels

class TestAssertUniqueLabelsForCompositeCreation(unittest.TestCase):

	def setUp(self):
		self.labelA = labelObjs.StandardLabel(eleKey="eleA",structKey="structA",methodKey="methB")
#		self.labelB = labelObjs.StandardLabel(eleKey="eleB",structKey="structB",methodKey="methB")

	def testCreatingCompositeWithEqualLabels(self):
		copiedLabelA = copy.deepcopy(self.labelA)
		with self.assertRaises(AssertionError):
			DudClassForTestingUniqueLabelEnforcement([self.labelA, copiedLabelA])



class TestTemporarilySetInstanceAttrs(unittest.TestCase):

	def setUp(self):
		self.testAttrA = "random_val_a"
		self.testAttrB = "random_val_b"
		self.createTestObjs()

	def createTestObjs(self):
		self.testInstanceA = types.SimpleNamespace(testAttrA=self.testAttrA, testAttrB=self.testAttrB)

	def testArgsAreSetWithinManagerThenReset(self):
		newAttrAVal = "new_attr_a_val"
		self.assertNotEqual(newAttrAVal, self.testInstanceA.testAttrA)
		self.assertEqual(self.testAttrA,self.testInstanceA.testAttrA)
		kwargDict = {"testAttrA":newAttrAVal}
		with tCode.temporarilySetInstanceAttrs(self.testInstanceA, kwargDict):
			self.assertEqual(newAttrAVal, self.testInstanceA.testAttrA)
		self.assertEqual(self.testAttrA, self.testInstanceA.testAttrA)

	#Probably NEVER end up changin behavior even as kwarg-option. It would be very tough to figure out
	#what to set the attr to afterwards
	def testDynamicAttrSettingRaisesAttribErrorByDefault(self):
		fakeKey = "fakeKey"
		with self.assertRaises(AttributeError):
			with tCode.temporarilySetInstanceAttrs(self.testInstanceA,{fakeKey:None}):
				pass



if __name__ == '__main__':
	unittest.main()

