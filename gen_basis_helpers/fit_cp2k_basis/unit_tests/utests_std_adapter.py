
import unittest
import unittest.mock as mock

import gen_basis_helpers.fit_cp2k_basis.adapter_stdinp as tCode


class TestGetAdapaterFromWorkflow(unittest.TestCase):

	def setUp(self):
		self.workflowA = mock.Mock()
		self.outputToObjFunctA = mock.Mock()
		self.labelA = None

	@mock.patch("gen_basis_helpers.fit_cp2k_basis.adapter_stdinp.calcRunners.StandardInputObj")
	@mock.patch("gen_basis_helpers.fit_cp2k_basis.adapter_stdinp.GenericMapFunctionForGettingObjFunctFromStdInp")
	def testExpectedCallsMade(self, mockedGenericMapFunct, mockedStdInp):
		expMapFunct, expStdInp = mock.Mock(), mock.Mock()
		mockedGenericMapFunct.side_effect = lambda *args: expMapFunct
		mockedStdInp.side_effect = lambda *args, **kwargs: expStdInp		

		actStdInp = tCode.getAdaptedStdInptFromWorkflow(self.workflowA, self.outputToObjFunctA, label=self.labelA)
		mockedGenericMapFunct.assert_called_with(self.outputToObjFunctA)
		mockedStdInp.assert_called_with(self.workflowA, self.labelA, mapFunction=expMapFunct)
		self.assertEqual(expStdInp,actStdInp)

