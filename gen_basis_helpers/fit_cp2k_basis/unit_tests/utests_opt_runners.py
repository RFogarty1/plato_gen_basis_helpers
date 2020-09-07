
import types
import unittest
import unittest.mock as mock

import gen_basis_helpers.fit_cp2k_basis.core as coreHelp

import gen_basis_helpers.fit_cp2k_basis.opt_runners as tCode

class TestBasicOptRunner(unittest.TestCase):

	def setUp(self):
		self.transformerA = lambda x:[a*2 for a in x]
		self.contribObjFunctA = mock.Mock()
		self.contribObjFunctA.runComms = list()
		self.startCoeffsA = [1,2,4]
		self.createTestObjs()

	def createTestObjs(self):
		self.coeffUpdaterA = coreHelp.CoeffUpdaterStandard(transformer=self.transformerA)
		self.objFunctA = coreHelp.ObjFunctCalculatorStandard( [self.contribObjFunctA],self.coeffUpdaterA)


	@mock.patch("gen_basis_helpers.fit_cp2k_basis.core.ObjFunctCalculatorStandard._calcTotalObjFunct")
	@mock.patch("gen_basis_helpers.fit_cp2k_basis.opt_runners.minimize")
	def testTransformedCoeffsAsExpectedWhenTransformerSet(self, mockedMinimize, mockedCalcTotalObjFunct):
		rawCoeffs = [1,2,3]
		expTransformedCoeffs = [2,4,6]
		mockedMinimize.side_effect = lambda *args,**kwargs: types.SimpleNamespace(x=rawCoeffs)
		actOutput = tCode.carryOutOptimisationBasicOptions(self.objFunctA, self.startCoeffsA)
		self.assertEqual(expTransformedCoeffs, actOutput.transformedCoeffs)

	@mock.patch("gen_basis_helpers.fit_cp2k_basis.core.ObjFunctCalculatorStandard._calcTotalObjFunct")
	@mock.patch("gen_basis_helpers.fit_cp2k_basis.opt_runners.minimize")
	def testTransformedCoeffsAsExpectedWhenNotSettingTransformer(self, mockedMinimize, mockedCalcTotalObjFunct):
		self.transformerA = None
		self.createTestObjs()
		rawCoeffs = [1,2,3]
		mockedMinimize.side_effect = lambda *args,**kwargs: types.SimpleNamespace(x=rawCoeffs)
		actOutput = tCode.carryOutOptimisationBasicOptions(self.objFunctA, self.startCoeffsA)
		self.assertEqual(rawCoeffs, actOutput.transformedCoeffs)


