
import math
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

