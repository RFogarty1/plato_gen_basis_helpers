
import copy

import numpy as np

import unittest
import unittest.mock as mock

import gen_basis_helpers.workflows.elastic_workflows as tCode



class TestStrainObject(unittest.TestCase):

	def setUp(self):
		self.strainValsA = [1,2,3,4,5,6]
		self.strainValsB = [2,2,2,3,3,3]
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.CrystalStrain(self.strainValsA)
		self.testObjB = tCode.CrystalStrain(self.strainValsB)

	@mock.patch("gen_basis_helpers.workflows.elastic_workflows._getUnitStrainMatrix")
	def testExpectedStrainMatrixWithMocking(self, mockedGetter):
		mockedGetter.side_effect = lambda x: x if (x>0 and x<7) else None
		expStrainMatrix = sum([x*x for x in self.strainValsA]) #In reality its the sum of some matrices
		actStrainMatrix = self.testObjA.strainMatrix
		self.assertEqual(expStrainMatrix, actStrainMatrix)

	def testExpectedStrainMatrixWithoutMocking(self):
		expMatrix = np.array( [ [1,3,2.5], [3,2,2], [2.5,2,3] ] )
		actMatrix = self.testObjA.strainMatrix
		self.assertTrue(np.allclose(expMatrix,actMatrix))

	def testEquality_equalObjsCompareEqual(self):
		self.strainValsA = [x*0.1 for x in self.strainValsA] #Want non-integers
		self.strainValsB = list(self.strainValsB)
		self.createTestObjs()
		self.assertEqual( self.testObjA, self.testObjB )

	def testEquality_unequalObjsCompareUnequal(self):
		self.assertNotEqual( self.testObjA, self.testObjB )

	def testToStrAsExpectedWhenAllSet(self):
		expStr = "+".join(["{}eps{}".format(x,idx) for idx,x in enumerate(self.strainValsB,1)])
		actStr = self.testObjB.toStr()
		self.assertEqual(expStr,actStr)

	def testToStrAsExpectedWhenTwoSet(self):
		self.strainValsA = [0,1,0,0,2,0]
		self.createTestObjs()
		expStr = "1eps2+2eps5"
		actStr = self.testObjA.toStr()
		self.assertEqual(expStr,actStr)





