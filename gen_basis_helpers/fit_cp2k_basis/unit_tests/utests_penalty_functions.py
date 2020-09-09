

import unittest 
import unittest.mock as mock

import gen_basis_helpers.fit_cp2k_basis.penalty_functions as tCode

class TestMaximumCoeffPenaltyFunct(unittest.TestCase):

	def setUp(self):
		self.maxValue = 4
		self.penaltyIfOver = None
		self.coeffsA = [5,6,3]
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.MaxAbsValCoeffPenaltyFunctionStandard(self.maxValue, penaltyIfOver=self.penaltyIfOver)

	def testExpValGivenForCoeffsA(self):
		expVal = 3
		actVal = self.testObjA(self.coeffsA)
		self.assertEqual(expVal, actVal)

	def testExpValGivenForCustomPenaltyFunct(self):
		self.penaltyIfOver = lambda targVal,actVal: 2*(abs(targVal-actVal))
		self.createTestObjs()
		expVal = 6
		actVal = self.testObjA(self.coeffsA)
		self.assertEqual(expVal,actVal)

class TestGetAdaptedStdInptFromPenaltyFunct(unittest.TestCase):

	def setUp(self):
		self.functRetVal = 3
		self.penaltyFunct = lambda coeffs:self.functRetVal
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.getAdaptedStdInptFromPenaltyFunct(self.penaltyFunct)

	def testCorrectObjectiveFunctionValGiven(self):	
		outputObj = self.testObjA.createOutputObj()
		expObjFunctVal = self.functRetVal
		actObjFunctVal = outputObj.data[0].objFunct #Is this really the correct one??
		self.assertEqual(expObjFunctVal, actObjFunctVal)

	def testRunCommsReturnsEmptyList(self):
		expRunComms = list()
		actRunComms = self.testObjA.runComms
		self.assertEqual(expRunComms,actRunComms)
