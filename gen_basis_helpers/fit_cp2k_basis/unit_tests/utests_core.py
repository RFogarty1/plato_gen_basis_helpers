
import types
import unittest
import unittest.mock as mock

import gen_basis_helpers.fit_cp2k_basis.core as tCode

class TestCoeffUpdaterStandard(unittest.TestCase):

	def setUp(self):
		self.observerA = mock.Mock()
		self.transformerA = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjNoTransform = tCode.CoeffUpdaterStandard()
		self.testObjWithTransform = tCode.CoeffUpdaterStandard(transformer=self.transformerA)
		self.testObjNoTransform.addObserver(self.observerA)
		self.testObjWithTransform.addObserver(self.observerA)

	def testTransmitsWithNoTransformFunction(self):
		testCoeffs = [1,2,3]
		self.testObjNoTransform(testCoeffs)
		self.observerA.updateCoeffs.assert_called_with(testCoeffs)

	def testTransformation(self):
		inpCoeffs = [1,2,3]
		expCoeffs = mock.Mock()
		self.transformerA.side_effect = lambda *args: expCoeffs
		self.testObjWithTransform(inpCoeffs)
		self.transformerA.assert_called_with(inpCoeffs)
		self.observerA.updateCoeffs.assert_called_with(expCoeffs)
		

class TestObjFunctCalculator(unittest.TestCase):

	def setUp(self):
		self.objA = mock.Mock()
		self.objB = mock.Mock()
		self.coeffUpdaterA = mock.Mock()
		self.nCores = 4
		self.createTestObjs()

	def createTestObjs(self):
		self.objs = [self.objA, self.objB]
		self.testObjA = tCode.ObjFunctCalculatorStandard(self.objs, self.coeffUpdaterA, nCores=self.nCores)

	def testCoeffsUpdaterCalled(self):
		self.testCoeffs = [1,2,3]
		self.testObjA._updateCoeffs(self.testCoeffs)
		self.coeffUpdaterA.assert_called_with(self.testCoeffs)

	@mock.patch("gen_basis_helpers.fit_cp2k_basis.core.jobRunHelp.executeRunCommsParralel")
	def testPreRunShellComms(self, mockedRunner):
		runCommsA, runCommsB = [mock.Mock()], [mock.Mock()]
		expRunComms = runCommsA + runCommsB
		self.objA.runComms, self.objB.runComms = runCommsA, runCommsB
		self.testObjA._doPreRunShellComms()
		mockedRunner.assert_called_with(expRunComms, self.nCores, quiet=True, noCommsOk=True)
		
	def testExpectedObjFunctionValsReturned(self):
		expVals = [2,3]
		outputObjA, outputObjB = mock.Mock(), mock.Mock()
		outputObjA.data = [types.SimpleNamespace(objFunct=expVals[0])]
		outputObjB.data = [types.SimpleNamespace(objFunct=expVals[1])]
		self.objA.createOutputObj.side_effect = lambda: outputObjA
		self.objB.createOutputObj.side_effect = lambda: outputObjB
		actVals = self.testObjA._getObjFunctVals()
		self.assertEqual(expVals, actVals)

	def testCombineObjFunct(self):
		testInput = [3,4]
		expOutput = sum(testInput)
		actOutput = self.testObjA._combineObjFunctVals(testInput)
		self.assertEqual(expOutput, actOutput)

	@mock.patch("gen_basis_helpers.fit_cp2k_basis.core.ObjFunctCalculatorStandard._calcTotalObjFunct")
	@mock.patch("gen_basis_helpers.fit_cp2k_basis.core.ObjFunctCalculatorStandard._doPreRunShellComms")
	@mock.patch("gen_basis_helpers.fit_cp2k_basis.core.ObjFunctCalculatorStandard._updateCoeffs")
	def testExpectedCallsMade(self, mockedUpdateCoeffs, mockedPreRunComms, mockedCalcObjFunct):
		testCoeffs = [1,2]
		expOutput = mock.Mock()
		mockedCalcObjFunct.side_effect = lambda: expOutput
		actOutput = self.testObjA(testCoeffs) 
		mockedUpdateCoeffs.assert_called_with(testCoeffs)
		mockedPreRunComms.assert_called_with()
		mockedCalcObjFunct.assert_called_with()
		self.assertEqual(expOutput,actOutput)

		



