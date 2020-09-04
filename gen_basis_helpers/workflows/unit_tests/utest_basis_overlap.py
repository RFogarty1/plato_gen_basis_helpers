
import itertools as it
import math
import types
import unittest
import unittest.mock as mock

import plato_pylib.plato.parse_gau_files as parseGau
import gen_basis_helpers.workflows.basis_overlap_workflow as tCode

class TestBasisFunctSelfOverlapWorkflow(unittest.TestCase):

	def setUp(self):
		self.distA = 0
		self.angMom = 0
		self.exponentsA = [1,2] 
		self.coeffsA = [ [1,2] ]
		self.createTestObjs()

	def createTestObjs(self):
		self.basisObjA = parseGau.GauPolyBasis(self.exponentsA,self.coeffsA)
		self.testObjA = tCode.BasisFunctSelfOverlapAtDistWorkflow(self.basisObjA, self.angMom, self.distA)

	def testRaisesForTooHighAngMom(self):
		self.angMom = 1 #To start i'm only doing the l=0 case
		self.createTestObjs()
		with self.assertRaises(ValueError):
			self.testObjA.run()

	def testRaisesForMultiplePoly(self):
		self.coeffsA = [ [1,2], [3,4] ]
		self.createTestObjs()
		with self.assertRaises(AssertionError):
			self.testObjA.run()

	@mock.patch("gen_basis_helpers.workflows.basis_overlap_workflow.sIntHelp.getSelfOverlapMcWedaWeightFromGauPolyBasis")
	def testRunCallsExpectedFuncts(self, mockedGetOverlapInt):
		expOverlap = 3
		mockedGetOverlapInt.side_effect = lambda *args: expOverlap
		self.testObjA.run()
		actOverlap = self.testObjA.output[0].overlap
		mockedGetOverlapInt.assert_called_with(self.basisObjA, self.distA)
		self.assertEqual(expOverlap,actOverlap)

	def testRunGivesExpectedForDistZero(self):
		self.exponentsA = [1]
		normCoeff = math.sqrt( 1 / (math.sqrt( math.pi/(2) )**3) )
		self.coeffsA = [ [normCoeff] ]
		self.createTestObjs()
		self.testObjA.run()
		expVal = 1
		actVal = self.testObjA.output[0].overlap
		self.assertAlmostEqual(expVal,actVal)


class TestSelfOverlapCoeffUpdater(unittest.TestCase):

	def setUp(self):
		self.workflowA = mock.Mock()
		self.coeffToPolyMapper = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.BasisFunctSelfOverlapCoeffUpdater(self.workflowA, self.coeffToPolyMapper)

	def testUpdateCoeffsCallsExpectedFuncts(self):
		testCoeffs = [1,2,3]
		expBasisObj = mock.Mock()
		self.workflowA.basisObj = mock.Mock()
		self.coeffToPolyMapper.side_effect = lambda coeffs: expBasisObj
		self.assertNotEqual(expBasisObj,self.workflowA.basisObj)

		self.testObjA.updateCoeffs(testCoeffs)
		self.coeffToPolyMapper.assert_called_with(testCoeffs)
		self.assertEqual(expBasisObj, self.workflowA.basisObj)



class TestSelfOverlapWorkflowOutputToObjFunct(unittest.TestCase):

	def setUp(self):
		self.targOverlapsA = [4]
		self.targOverlapsB = [5,7]
		self.overlapsA = [2]
		self.overlapsB = [2,3]
		self.targActMapFunct = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		self.outputA = [types.SimpleNamespace(overlap=x) for x in self.overlapsA]
		self.outputB = [types.SimpleNamespace(overlap=x) for x in self.overlapsB]
		self.testObjA = tCode.SelfOverlapWorkflowToObjFunctValStandard(self.targOverlapsA, self.targActMapFunct)
		self.testObjB = tCode.SelfOverlapWorkflowToObjFunctValStandard(self.targOverlapsB, self.targActMapFunct)

	def testExpectedValuesFromOutputSingleVal(self):
		expOverlaps = self.overlapsA
		actOverlaps = self.testObjA._getValuesFromOutput(self.outputA)
		self.assertEqual(expOverlaps, actOverlaps)

	def testExpectedValuesFromOutputMultiVals(self):
		expOverlaps = self.overlapsB
		actOverlaps = self.testObjB._getValuesFromOutput(self.outputB)
		self.assertEqual(expOverlaps, actOverlaps)

	def testExpectedOutputB(self):
		self.targActMapFunct.side_effect = lambda targVals, actVals: sum([abs(x-y) for x,y in it.zip_longest(targVals,actVals)])
		expVal = 7
		actVal = self.testObjB(self.outputB)
		self.assertEqual(expVal,actVal)




