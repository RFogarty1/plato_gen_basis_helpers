

import os
import unittest
import unittest.mock as mock

import gen_basis_helpers.job_helpers.conv_prop_total_energy as tCode

class TestEnergyConvergenceTemplate(unittest.TestCase):

	def setUp(self):
		self.eleKey = "MgO"
		self.methodKey = "cutoff"
		self.structKey = "exptGeom"
		self.baseWorkFolder = "fake_workfolder"
		self.convVals = [1,2,3]
		self.convKwarg = "fake_conv_kwarg"
		self.baseCreatorObj = mock.Mock()
		self.geom = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		self.expOutFolder = os.path.join(self.baseWorkFolder, self.eleKey, self.structKey, self.methodKey)
		self.mockCreators = [mock.Mock() for x in self.convVals]
		kwargDict = { "eleKey":self.eleKey,
		              "methodKey":self.methodKey,
		              "structKey":self.structKey,
		              "baseWorkFolder":self.baseWorkFolder,
		              "convVals":self.convVals,
		              "geom":self.geom,
		              "convKwarg":self.convKwarg,
		              "baseCreatorObj":self.baseCreatorObj }
		self.testObjA = tCode.CodeSpecificStandardInputCreatorTemplate(**kwargDict)

	def testOutfolderAsExpected(self):
		expPath = self.expOutFolder
		actPath = self.testObjA.outFolder
		self.assertEqual(expPath,actPath)

	def testExpectedFileNamesGenerated(self):
		expBaseFileNames = ["conv_val_{:.3f}".format(x).replace(".","pt")  for x in self.convVals]
		actBaseFileNames = self.testObjA._outBaseFileNames
		self.assertEqual(expBaseFileNames, actBaseFileNames)

	def testCreatorObjModifiedWithSharedValues(self):
		testObj = mock.Mock()
		self.testObjA._modifyCreatorObjWithSharedOptions(testObj)
		expGeom, expBaseFolder = self.geom, self.expOutFolder
		self.assertEqual(expGeom, testObj.geom)
		self.assertEqual(expBaseFolder, testObj.workFolder)


	@mock.patch("gen_basis_helpers.job_helpers.conv_prop_total_energy.CodeSpecificStandardInputCreatorTemplate._getBasicCreatorObjForEachRequiredCalc")
	@mock.patch("gen_basis_helpers.job_helpers.conv_prop_total_energy.CodeSpecificStandardInputCreatorTemplate._outBaseFileNames",new_callable=mock.PropertyMock)
	def testOutFileNamesAppliedToCreators(self, mockedFileNames, mockedGetCreators):
		expFileNames = [mock.Mock() for x in self.convVals]
		mockedFileNames.return_value = expFileNames
		mockedGetCreators.side_effect = lambda *args,**kwargs: self.mockCreators 

		outCreators = self.testObjA._getCreatorsForWorkflow()
		actFileNames = [x.fileName for x in outCreators]

		self.assertEqual(outCreators, self.mockCreators)
		self.assertEqual(expFileNames, actFileNames)

	@mock.patch("gen_basis_helpers.job_helpers.conv_prop_total_energy.CodeSpecificStandardInputCreatorTemplate._modifyCreatorObjWithSharedOptions")
	@mock.patch("gen_basis_helpers.job_helpers.conv_prop_total_energy.CodeSpecificStandardInputCreatorTemplate._getBasicCreatorObjForEachRequiredCalc")
	def testSharedOptionsAppliedToCreators(self, mockedGetCreators, mockedModWithSharedOptions):
		mockedGetCreators.side_effect = lambda *args,**kwargs: self.mockCreators
		self.testObjA._getCreatorsForWorkflow()
		for x in self.mockCreators:
			mockedModWithSharedOptions.assert_any_call(x)

	@mock.patch("gen_basis_helpers.job_helpers.conv_prop_total_energy.CodeSpecificStandardInputCreatorTemplate._getBasicCreatorObjForEachRequiredCalc")
	def testConvValsAppliedToCreators(self, mockedGetCreators):
		mockedGetCreators.side_effect = lambda *args,**kwargs:self.mockCreators
		outCreators = self.testObjA._getCreatorsForWorkflow()
		expConvVals = self.convVals
		actConvVals = [getattr(x, self.convKwarg) for x in outCreators]
		self.assertEqual(expConvVals,actConvVals)

	def testRaisesWhenConvValsNotUnique(self):
		self.convVals[1] = self.convVals[0]
		self.createTestObjs()
		with self.assertRaises(ValueError):
			self.testObjA.create()

	@mock.patch("gen_basis_helpers.job_helpers.conv_prop_total_energy.CodeSpecificStandardInputCreatorTemplate._getCalcObjs")
	@mock.patch("gen_basis_helpers.job_helpers.conv_prop_total_energy.convFlow.GridConvergenceEnergyWorkflow")
	def testWorkflowCalledWithExpectedArgs(self, mockedGridConvFlow, mockedCalcObjs):
		expCalcObjs = [mock.Mock() for x in self.convVals]
		mockedCalcObjs.side_effect = lambda *args, **kwargs: expCalcObjs
		self.testObjA.create()
		mockedGridConvFlow.assert_called_once_with( expCalcObjs, self.convVals )
		
