
import itertools as it
import math
import unittest
import unittest.mock as mock

import plato_pylib.plato.parse_gau_files as parseGau

import gen_basis_helpers.fit_cp2k_basis.basis_coeff_mappers as tCode


class TestCoeffsToExponentsAndCoeffsStandard(unittest.TestCase):

	def setUp(self):
		self.fixedCoeffs = list()
		self.fixedExponents = [5]
		self.testCoeffsA = [1,2,3]
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.FitCoeffsToBasisFunctionExponentsAndCoeffsMixedOptStandard(self.fixedCoeffs, self.fixedExponents)

	def _runFunct(self):
		return self.testObjA(self.testCoeffsA)

	def testExpectedForNoFixedCoeffs(self):
		expExponents = self.fixedExponents + [self.testCoeffsA[0]]
		expCoeffs = self.testCoeffsA[1:]
		actExponents, actCoeffs  = self._runFunct()
		self.assertEqual(expCoeffs, actCoeffs)
		self.assertEqual(expExponents, actExponents)

	def testExpectedForAllFree(self):
		self.fixedCoeffs, self.fixedExponents = list(), list()
		expExponents, expCoeffs = [1,2], [3,4]
		self.testCoeffsA = expExponents + expCoeffs
		self.createTestObjs()
		actExponents, actCoeffs = self._runFunct()
		self.assertEqual(expExponents, actExponents)
		self.assertEqual(expCoeffs, actCoeffs)

	def testExpecetedForMixOfFixedCoeffsAndExponents(self):
		self.fixedExponents, self.fixedCoeffs = [1],[3]
		self.testCoeffsA = [2,4]
		self.createTestObjs()
		expExponents, expCoeffs = [1,2], [3,4]
		actExponents, actCoeffs = self._runFunct()
		self.assertEqual(expExponents, actExponents)
		self.assertEqual(expCoeffs, actCoeffs)

	def testExpectedForNoFixedExponents(self):
		self.fixedExponents, self.fixedCoeffs = list(), [4]
		self.testCoeffsA = [1,2,3]
		expExponents, expCoeffs = [1,2], [4,3]
		self.createTestObjs()
		actExponents, actCoeffs = self._runFunct()
		self.assertEqual(expExponents, actExponents)
		self.assertEqual(expCoeffs, actCoeffs)

	def testRaisesForOddNumberOfCoefficients(self):
		""" Check we throw an error if we get an odd number of total coefficients """
		self.testCoeffsA = [1,2]
		with self.assertRaises(AssertionError):
			self.assertTrue(False)


class TestCoeffsToBasisFunctionMixedOpt(unittest.TestCase):

	def setUp(self):
		self.fixedExponents = [3]
		self.fixedCoeffs = [1,2]
		self.testCoeffsA = [4]
		self.angMomA = 0
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.CoeffsToFullBasisFunctionForMixedCoeffExponentOpt(self.fixedExponents, self.fixedCoeffs, self.angMomA)

	@mock.patch("gen_basis_helpers.fit_cp2k_basis.basis_coeff_mappers.parseGau.GauPolyBasis")
	def testExpectedCallsToGauPoly(self, mockedGauPoly):
		expObj = mock.Mock()
		expExponents, expCoeffs = [3,4], [1,2]
		mockedGauPoly.side_effect = lambda *args,**kwargs: expObj
		actObj = self.testObjA(self.testCoeffsA)
		mockedGauPoly.assert_called_with(expExponents, [expCoeffs], label=self.angMomA)
		self.assertEqual(expObj, actObj)

class TestCoeffsToBasisFunctionForFixedExponents(unittest.TestCase):

	def setUp(self):
		self.exponentsA = [1,2]
		self.angMomA = 0
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.CoeffsToFullBasisFunctionForFixedExponents(self.exponentsA, self.angMomA)

	def testExpectedForTestInputA(self):
		testCoeffsA = [3,4]
		expGauPolyBas = parseGau.GauPolyBasis(self.exponentsA, [testCoeffsA], label=self.angMomA)
		actGauPolyBas = self.testObjA(testCoeffsA)
		self.assertEqual(expGauPolyBas, actGauPolyBas)


class TestCoeffsToFullBasisNewFunctionCoeffsOpt(unittest.TestCase):

	def setUp(self):
		self.exponentsA = [1,2]
		self.angMomA = 0
		self.fixedOrbExpansions = [mock.Mock()]
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.CoeffsToFullBasisSetForFixedExponents(self.fixedOrbExpansions, self.exponentsA, self.angMomA)

	def testExpectedForTestInput(self):
		testCoeffsA = [3,4]
		expGauPolyBas = parseGau.GauPolyBasis(self.exponentsA, [testCoeffsA], label=self.angMomA)
		expOutput = self.fixedOrbExpansions + [expGauPolyBas]
		actOutput = self.testObjA(testCoeffsA)
		self.assertEqual(expOutput, actOutput)

class TestCoeffsToNormalisedValuesForFixedExponents(unittest.TestCase):

	def setUp(self):
		self.exponentsA = [2]
		self.coeffsA = [3]
		self.angMomA = 0
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.CoeffsToNormalisedValuesFixedExponents(self.exponentsA, self.angMomA)

	def testRaisesForAngMomNotZero(self):
		self.angMomA = 1
		self.createTestObjs()
		with self.assertRaises(ValueError):
			self.testObjA(self.coeffsA)

	def testGetScaleFactorAsExpected(self):
		areaOfGaussian = (math.pi/(2*self.exponentsA[0]))**(3/2)
		expScaleFactor = 1/ math.sqrt(self.coeffsA[0]*self.coeffsA[0]*areaOfGaussian)
		actScaleFactor = self.testObjA._getScaleFactor(self.coeffsA)
		self.assertAlmostEqual(expScaleFactor,actScaleFactor)

	@mock.patch("gen_basis_helpers.fit_cp2k_basis.basis_coeff_mappers.CoeffsToNormalisedValuesFixedExponents._getScaleFactor")
	def testScaleFactorAppliedOnCall(self, mockedGetScaleFactor):
		expScaleFactor = 2
		mockedGetScaleFactor.side_effect = lambda coeffs: expScaleFactor
		expCoeffs = [x*expScaleFactor for x in self.coeffsA]
		actCoeffs = self.testObjA(self.coeffsA)
		mockedGetScaleFactor.assert_called_with(self.coeffsA)
		self.assertEqual(expCoeffs, actCoeffs)

	def testIdemptotent(self):
		""" Applying twice to any set of coeffs should give the same result as applying once """
		#Check that applying it once changes the coeffs
		coeffsB = list( self.testObjA(self.coeffsA) )
		expCoeffs = list(coeffsB)
		for cA,cB in it.zip_longest(self.coeffsA, coeffsB):
			self.assertNotAlmostEqual(cA,cB)

		#Check that applying the second time does NOT change the coeffs
		actCoeffs = self.testObjA(coeffsB)
		for exp,act in it.zip_longest(expCoeffs,actCoeffs):
			self.assertAlmostEqual(exp,act)



