

import unittest
import unittest.mock as mock


import gen_basis_helpers.shared.label_objs as labelHelp
import gen_basis_helpers.job_helpers.eos_curves as tCode

class TestEosCurvesCreator(unittest.TestCase):

	def setUp(self):
		self.eleKey, self.methodKey, self.structKey = "fake_ele", "fake_method", "fake_struct"
		self.structStrs = ["hcp","bcc"]
		self.baseCreators = mock.Mock()
		self.geoms = [[mock.Mock()],[mock.Mock()]]
		self.eosFitFuncts = [ mock.Mock(), mock.Mock() ]
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"eleKey":self.eleKey, "methodKey":self.methodKey, "structKey":self.structKey,
		             "structStrs":self.structStrs, "baseCreators":self.baseCreators,
		             "geoms":self.geoms, "eosFitFunctions":self.eosFitFuncts}
		self.testObjA = tCode.StandardInputCreatorEosCurves(**kwargDict)

	def testGetIterFromAttr_nonIter(self):
		expIter = [self.baseCreators for x in self.structStrs]
		actIter = self.testObjA._getIterForAttr("baseCreators")
		self.assertEqual(expIter,actIter)

	def testGetIterFromAttr_iter(self):
		expIter = self.geoms
		actIter = self.testObjA._getIterForAttr("geoms")
		self.assertEqual(expIter,actIter)

	def testGetIterFromAttr_wrongLenIter(self):
		self.structStrs = ["hcp"]
		self.geoms = [ [mock.Mock()], [mock.Mock()] ]
		self.createTestObjs()

		with self.assertRaises(AssertionError):
			self.testObjA._getIterForAttr("geoms")

	def testGetCalcObjsForSingleWorkflow(self):
		structStr = "bcc"
		structIdx = 1

		self.baseCreators = [mock.Mock(), mock.Mock()]
		expOutput = [mock.Mock(), mock.Mock()]
		self.geoms[1] = [mock.Mock(), mock.Mock()]
		self.baseCreators[structIdx].create.side_effect = expOutput

		self.createTestObjs()

		actOutput = self.testObjA._getCalcObjsForStructStr(structStr)
		for x in self.geoms[1]:
			self.baseCreators[1].create.assert_any_call(geom=x)
		self.assertEqual(expOutput,actOutput)

#                gen_basis_helpers.job_helpers.eos_curves
	@mock.patch("gen_basis_helpers.job_helpers.eos_curves.StandardInputCreatorEosCurves._getCalcObjsForStructStr")
	@mock.patch("gen_basis_helpers.job_helpers.eos_curves.eosFlowHelp.EosWorkflow")
	def testExpectedWorkflowForSingleStructStr(self, mockedWorkflow, mockedGetCalcObjs):
		expStructKey = "bcc"
		expFitFunct = self.eosFitFuncts[1]
		expLabel = labelHelp.StandardLabel(eleKey=self.eleKey, structKey=expStructKey, methodKey=self.methodKey)
		expCalcObjs, expWorkflow = mock.Mock(), mock.Mock()
		mockedGetCalcObjs.side_effect = lambda *args:expCalcObjs 
		mockedWorkflow.side_effect = lambda *args,**kwargs: expWorkflow

		actOutput = self.testObjA._getWorkflowForStructStr(expStructKey)

		mockedGetCalcObjs.assert_called_with(expStructKey)
		mockedWorkflow.assert_called_with(expCalcObjs, expFitFunct, expLabel)
		self.assertEqual(expWorkflow, actOutput)

	@mock.patch("gen_basis_helpers.job_helpers.eos_curves.baseFlow.StandardLabelledWorkflowComposite")
	@mock.patch("gen_basis_helpers.job_helpers.eos_curves.StandardInputCreatorEosCurves._getWorkflowForStructStr")
	def testExpectedCompositeWorkflow(self, mockedSingleFlowGetter, mockedWorkflowComposite):
		expWorkflow = mock.Mock()
		expHcpFlow, expBccFlow = mock.Mock(), mock.Mock()

		mockedSingleFlowGetter.side_effect = [expHcpFlow, expBccFlow]
		mockedWorkflowComposite.side_effect = lambda *args, **kwargs: expWorkflow
		
		actOutput = self.testObjA._createWorkflow()
		for x in self.structStrs:
			mockedSingleFlowGetter.assert_any_call(x)
		mockedWorkflowComposite.assert_called_with([expHcpFlow, expBccFlow])
		self.assertEqual(expWorkflow,actOutput)

	@mock.patch("gen_basis_helpers.job_helpers.eos_curves.calcRunners.StandardInputObj")
	@mock.patch("gen_basis_helpers.job_helpers.eos_curves.StandardInputCreatorEosCurves._createWorkflow")
	def testExpectedCallToStandardInputCreator(self, mockedWorkflowGetter, mockedStandardInp):
		expLabel = labelHelp.StandardLabel(eleKey=self.eleKey, structKey=self.structKey, methodKey=self.methodKey)
		expStdInp, expWorkflow, expOutput = mock.Mock(), mock.Mock(), mock.Mock()

		mockedWorkflowGetter.side_effect = lambda *args,**kwargs: expWorkflow
		mockedStandardInp.side_effect = lambda *args,**kwargs: expStdInp
		expStdInp.create.side_effect = lambda *args,**kwargs: expOutput

		actOutput = self.testObjA.create()
		mockedWorkflowGetter.assert_called_with()
		mockedStandardInp.assert_called_with(expWorkflow, expLabel)
		self.assertEqual(expOutput, actOutput)	

