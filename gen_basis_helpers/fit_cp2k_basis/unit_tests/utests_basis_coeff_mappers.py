
import math
import unittest
import unittest.mock as mock

import plato_pylib.plato.parse_gau_files as parseGau

import gen_basis_helpers.fit_cp2k_basis.basis_coeff_mappers as tCode



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
		expScaleFactor = 1/(self.coeffsA[0]*self.coeffsA[0]*areaOfGaussian)
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



