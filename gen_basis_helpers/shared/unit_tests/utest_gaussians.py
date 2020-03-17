
import copy
import itertools as it
import unittest
import unittest.mock as mock

import gen_basis_helpers.shared.gaussians as tCode



class TestGauPrimComposite(unittest.TestCase):

	def setUp(self):
		self.leafA = mock.Mock()
		self.leafB = mock.Mock()
		self.leafC = mock.Mock()
		
		#Put default, simple functions on the mocks
		self.distFunctFactors = [1,2,3]

		self.mockedIntegralsOverSpace = [1,2,3]

		self.createTestObjs()

	def createTestObjs(self):

		#Mock the various functions as needed
		self.leafA.evalFunctAtDists.side_effect = lambda x: [self.distFunctFactors[0]*a for a in x]
		self.leafB.evalFunctAtDists.side_effect = lambda x: [self.distFunctFactors[1]*a for a in x]
		self.leafC.evalFunctAtDists.side_effect = lambda x: [self.distFunctFactors[2]*a for a in x]

		self.leafA.getIntegralAllSpace.side_effect = lambda : self.mockedIntegralsOverSpace[0]
		self.leafB.getIntegralAllSpace.side_effect = lambda : self.mockedIntegralsOverSpace[1]
		self.leafC.getIntegralAllSpace.side_effect = lambda : self.mocekdIntegralsOverSpace[2]

		self.leafA.leaves = [self.leafA]
		self.leafB.leaves = [self.leafB]
		self.leafC.leaves = [self.leafC]

		self.allLeaves = [self.leafA, self.leafB, self.leafC]
		self.compObjA = tCode.GauPrimComposite([self.leafA, self.leafB])
		self.compObjB = tCode.GauPrimComposite([self.compObjA,self.leafC])

	def testEvalAtDists_compWithOnlyLeaves(self):
		testDistances = [1,2,3]
		totalFactor = sum(self.distFunctFactors[:2])
		expVals = [totalFactor*x for x in testDistances]
		actVals = self.compObjA.evalFunctAtDists(testDistances)
		for exp,act in it.zip_longest(expVals,actVals):
			self.assertAlmostEqual(exp,act)

	def testEvalAtDists_compWithCompAndLeaves(self):
		testDistances = [1,2,3]
		totalFactor = sum(self.distFunctFactors)
		expVals = [totalFactor*x for x in testDistances]
		actVals = self.compObjB.evalFunctAtDists(testDistances)
		for exp,act in it.zip_longest(expVals,actVals):
			self.assertAlmostEqual(exp,act)

	def testGetIntegralAllSpace_compWithOnlyLeaves(self):
		expVal = sum(self.mockedIntegralsOverSpace[:2])
		actVal = self.compObjA.getIntegralAllSpace()
		self.assertEqual(expVal, actVal)

	#At this point i simply want this to be impossible; later i may add a way to set them as a list
	def testErrorWhenSettingCoeffs(self):
		with self.assertRaises(NotImplementedError):
			self.compObjA.c = 4

	def testErrorWhenSettingExponents(self):
		with self.assertRaises(NotImplementedError):
			self.compObjA.a = 4

	def testLeavesProperty(self):
		expLeaves = self.allLeaves
		actLeaves = self.compObjB.leaves
		self.assertEqual(expLeaves, actLeaves)

class TestGauCompositeEquality(unittest.TestCase):

	def setUp(self):
		self.expA = 0.1
		self.coeffA = 3
		self.posA = [0.0,1.0,0.0]
		self.createTestObjs()

	def createTestObjs(self):
		self.leafA = tCode.GauPrim( self.expA, self.coeffA, self.posA )
		self.compositeA = tCode.GauPrimComposite( [self.leafA] ) 

	def testLeafComparesEqualWithEquivComposite(self):
		self.assertEqual(self.leafA,self.compositeA)

	def testCompositeComparesEqualWithEquivLeaf(self):
		self.assertEqual(self.compositeA,self.leafA)

	def testLengthOneCompositeComparesUnequalWithNonEquivLeaf(self):
		objA = copy.deepcopy(self.compositeA)
		testCoeff = 5
		self.assertNotAlmostEqual(testCoeff,self.coeffA)
		self.coeffA = testCoeff
		self.createTestObjs()
		objB = self.leafA
		self.assertNotEqual(objA,objB)

	def testLeafNonEqualWithLenTwoComposite(self):
		objA = tCode.GauPrimComposite([self.leafA,self.leafA])
		objB = self.leafA
		self.assertNotEqual(objB, objA) #Order is important; we want to test the __eq__ on the leaf

	def testEquivCompositesCompareEqual_copyObj(self):
		objA = copy.deepcopy(self.compositeA)
		self.createTestObjs()
		objB = self.compositeA
		self.assertTrue(objA is not objB)
		self.assertEqual(objA,objB)

	def testDiffCompositesCompareDifferent_diffLengthLeaves(self):
		objA = copy.deepcopy(self.compositeA)
		objB = tCode.GauPrimComposite( [self.compositeA, self.leafA] )
		self.assertNotEqual(objA,objB)

	def testDiffCompositesCompareUnequal_diffValLeaf(self):
		objA = copy.deepcopy(self.compositeA)
		testExp = 0.5
		self.assertNotAlmostEqual(self.expA,testExp)
		self.expA = testExp
		self.createTestObjs()
		objB = self.compositeA
		self.assertNotEqual(objA,objB)

class TestGauPrimClass(unittest.TestCase):

	def setUp(self):
		self.expA = 0.1
		self.coeffA = 3
		self.posA = [0.0,1.0,0.0]

		self.expB = 0.2
		self.coeffB = 4
		self.posB = [0.0,1.0,2.0]
		self.createTestObjs()

	def createTestObjs(self):
		self.gauPrimA = tCode.GauPrim( self.expA, self.coeffA, self.posA )

	def testEvalFunctAtPositions(self):
		testPositions = [ [1.0,2.0,3.0], [2.0,0.0,2.0] ]
		expVals = [0.998613251094239, 1.2197089792218]
		actVals = [self.gauPrimA.evalFunctAtPos(x) for x in testPositions]

		for exp,act in it.zip_longest(expVals,actVals):
			self.assertAlmostEqual(exp,act)

	def testEvalFunctionAtDists(self):
		testDists = [1,2,3]
		expVals = [2.71451225410788, 2.01096013810692, 1.2197089792218]
		actVals = self.gauPrimA.evalFunctAtDists(testDists)
		for exp,act in it.zip_longest(expVals,actVals):
			self.assertAlmostEqual(exp,act)

	def testEqualObjsCompareEqual_copiedObj(self):
		objA = copy.deepcopy(self.gauPrimA)
		self.createTestObjs()
		objB = self.gauPrimA
		self.assertTrue(objA is not objB)
		self.assertEqual(objA,objB)		

	
	def testUnEqualObjsCompareUnEqual_diffExp(self):
		objA = copy.deepcopy(self.gauPrimA)
		testExpA = 0.5
		self.assertNotAlmostEqual(testExpA,self.expA)
		self.expA = testExpA
		self.createTestObjs()
		objB = self.gauPrimA
		self.assertNotEqual(objA,objB)

	def testUnequalCompareUnequal_diffPos(self):
		objA = copy.deepcopy(self.gauPrimA)
		self.posA[1] += 0.5
		self.createTestObjs()
		objB = self.gauPrimA
		self.assertNotEqual(objA,objB)


