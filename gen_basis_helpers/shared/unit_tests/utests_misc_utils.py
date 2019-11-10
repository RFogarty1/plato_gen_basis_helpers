#!/usr/bin/python3

import types
import unittest

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

if __name__ == '__main__':
	unittest.main()

