
import itertools as it
import unittest
import unittest.mock as mock

import gen_basis_helpers.shared.basis_set_objs as tCode

class TestGauSumOrbitalBasisFunction(unittest.TestCase):

	def setUp(self):
		self.nVal = 2
		self.lVal = 1
		self.gauFunct = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.GauSumOrbitalBasisFunction(self.nVal, self.lVal, self.gauFunct)

	def testCallsCorrectlyForGetRadialVals(self):
		testDistances = mock.Mock() #Values are super irrelevant
		fakeOutVals = mock.Mock()
		self.gauFunct.evalFunctAtDists.side_effect = lambda x: fakeOutVals
		actOutVals = self.testObjA.getRadialValsAtDists(testDistances)
		self.gauFunct.evalFunctAtDists.assert_called_once_with(testDistances)
		self.assertEqual(fakeOutVals,actOutVals)

	def testCorrectValsForPureRadial(self):
		self.gauFunct.evalFunctAtDists.side_effect = lambda x: x
		testXVals = [1,2]
		expVals = [2.04665341589298, 8.18661366357191]
		actVals = self.testObjA.getRadialValsAtDists(testXVals, pureRadial=True)
		for exp,act in it.zip_longest(expVals,actVals):
			self.assertAlmostEqual(exp,act)



class TestOrbitalBasisSetStandard(unittest.TestCase):

	def setUp(self):
		self.basisName = "fake_basis"
		self.lValA = 0
		self.lValB = 1
		self.lValC = 0
		self.createTestObjs()

	def createTestObjs(self):
		self.basisFunctA = mock.Mock()
		self.basisFunctB = mock.Mock()
		self.basisFunctC = mock.Mock()

		self.basisFunctA.lVal = self.lValA
		self.basisFunctB.lVal = self.lValB
		self.basisFunctC.lVal = self.lValC
	
		self.orbitalBasisFunctions = [self.basisFunctA, self.basisFunctB, self.basisFunctC]
		self.testObjA = tCode.OrbitalBasisSetStandard(self.orbitalBasisFunctions, self.basisName)

	def testExtractBasisFunctsWithL0(self):
		testLVal = 0
		expBasisFuncts = [self.basisFunctA, self.basisFunctC]
		actBasisFuncts=  self.testObjA.getBasisFunctionsWithLVal(testLVal)
		self.assertEqual(expBasisFuncts, actBasisFuncts)


