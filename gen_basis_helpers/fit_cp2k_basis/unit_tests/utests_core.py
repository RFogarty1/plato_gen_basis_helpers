
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
		

