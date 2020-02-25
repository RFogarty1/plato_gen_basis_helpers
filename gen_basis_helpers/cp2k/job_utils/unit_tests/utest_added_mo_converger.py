
import os
import unittest
import unittest.mock as mock


import gen_basis_helpers.cp2k.method_register as methRegister
import gen_basis_helpers.cp2k.cp2k_calc_objs as cp2kCalcObjs

import gen_basis_helpers.cp2k.job_utils.added_mo_conv as tCode

class TestAddedMoConverger(unittest.TestCase):

	def setUp(self):
		self.basicObject = methRegister.createCP2KObjFromMethodStr("cp2k_test_object")
		self.objCreationFunct = lambda x: cp2kCalcObjs.CP2KCalcObj( methRegister.createCP2KObjFromMethodStr("cp2k_test_object") )
		self.convVals = [1,2,3,5]
		self.fakeLabel = "heres_a_str"
		self.baseFolder = os.path.abspath( os.getcwd() )
		self.createTestObj()

	def createTestObj(self):
		self.testObjA = tCode.AddedMoConverger(self.baseFolder, self.objCreationFunct, self.convVals, self.fakeLabel)

	@mock.patch("gen_basis_helpers.cp2k.job_utils.added_mo_conv.AddedMoConverger.createStandardInputObjFromCalcObjs")
	def testConvParamsCorrect(self, mockedStdInpFromCalcObjs):
		mockedStdInpFromCalcObjs.side_effect = lambda calcObjs: calcObjs #self argument implicitly passed i think
		outObjs = self.testObjA.createStandardInput() 
		outConvVals = [x.addedMOs for x in outObjs]
		self.assertEqual( self.convVals, outConvVals )

	@mock.patch("gen_basis_helpers.cp2k.job_utils.added_mo_conv.AddedMoConverger.createStandardInputObjFromCalcObjs")
	def testExpPathsCorrect(self, mockedStdInpFromCalcObjs):
		mockedStdInpFromCalcObjs.side_effect = lambda calcObjs: calcObjs #self argument implicitly passed i think
		outObjs = self.testObjA.createStandardInput()
		expPaths = [os.path.join(self.baseFolder, "added_mos_{}.inp".format(x)) for x in self.convVals] 
		actPaths = [x.inpPath for x in outObjs]
		self.assertEqual(expPaths,actPaths)

#	@mock.patch("gen_basis_helpers.cp2k.job_utils.added_mo_conv.convFlow.GridConvergenceEnergyWorkflow")
#	def testCorrectParamsPassedToMockWorkflow(self, mockedWorkflow):
#		outObjs = self.testObjA.createStandardInput()
#		mockedWorkflow.assert_called_once_with(
#		self.assertTrue(False)

if __name__ == '__main__':
	unittest.main()


