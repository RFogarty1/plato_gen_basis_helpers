
import itertools as it
import types
import unittest
import unittest.mock as mock

import gen_basis_helpers.shared.label_objs as labelHelp
import gen_basis_helpers.shared.calc_runners as tCode

class TestStandardInputObj(unittest.TestCase):

	def setUp(self):
		self.preRunShellCommsA = ["fake_comm_a"]
		self.eleKey = "Mg"
		self.methodKey = "fake_method_a"
		self.structKey = "fake_struct_a"
		self.mapFunction = None
		self.output = types.SimpleNamespace(convVals=[(1,2),(2,4)])
		self.createTestObjs()

	def createTestObjs(self):
		self.mockWorkflowA = mock.Mock()
		self.mockWorkflowA.preRunShellComms = self.preRunShellCommsA
		self.mockWorkflowA.output = self.output
		self.testLabelA = labelHelp.StandardLabel(eleKey=self.eleKey, structKey=self.structKey,
		                                          methodKey=self.methodKey)
		self.testObjA = tCode.StandardInputObj( self.mockWorkflowA, self.testLabelA, mapFunction=self.mapFunction )

	def testRunCommsExpected(self):
		expRunComms = self.preRunShellCommsA
		actRunComms = self.testObjA.runComms
		self.assertEqual(expRunComms, actRunComms)

	def testLabelGetterGivesExpected(self):
		expVal = [self.testLabelA]
		actVal = self.testObjA.label
		self.assertEqual(expVal, actVal)

	@mock.patch("gen_basis_helpers.shared.calc_runners.StandardOutputObj")
	def testExpectedArgsPassedToOutputObjInit(self, mockedOutputClass):
		self.testObjA.createOutputObj()
		self.mockWorkflowA.run.assert_called_with()
		mockedOutputClass.assert_called_with(self.output, self.testLabelA)

	@mock.patch("gen_basis_helpers.shared.calc_runners.StandardOutputObj")
	def testCreateOutputObjWithCustomMapFunction(self, mockedOutputClass):
		#Modify our test input object
		fakeOutput = types.SimpleNamespace(convVals=[(2*x,y) for x,y in self.output.convVals])
		self.mapFunction = lambda x: fakeOutput
		self.createTestObjs()

		#Check our hacked version acts as expected
		self.testObjA.createOutputObj()
		mockedOutputClass.assert_called_with(fakeOutput, self.testLabelA)

