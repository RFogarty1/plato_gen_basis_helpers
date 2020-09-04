
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



