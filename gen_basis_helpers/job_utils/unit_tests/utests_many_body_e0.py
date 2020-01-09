#!/usr/bin/python3


import types
import unittest
import unittest.mock as mock

import gen_basis_helpers.job_utils.many_body_e0_helpers as tCode

class TestManyBodyCorrWorkflow(unittest.TestCase):

	def setUp(self):
		self.fullCalcObj = types.SimpleNamespace( e0=5, workFolder=6, writeFiles=mock.Mock() )
		self.twoBodyCalcObj = types.SimpleNamespace( e0 = 7, workFolder=6, writeFiles=mock.Mock() )
		self.createTestObj()	

	def createTestObj(self):
		self.testObj = tCode.ManyBodyCorrWorkflow(self.fullCalcObj, self.twoBodyCalcObj)

	def testExpectedManyBodyCorrXcObtained(self):
		expVal = -2 #Many body = Total minus 2-body
		self.testObj.run() #Assume jobs already been run
		actVal = self.testObj.output.manyBodyE0Xc
		self.assertEqual(expVal, actVal)

	def testErrorThrownForVaryingWorkFolders(self):
		self.fullCalcObj.workFolder = self.twoBodyCalcObj.workFolder*2
		with self.assertRaises(AssertionError):
			self.createTestObj()


if __name__ == '__main__':
	unittest.main()



