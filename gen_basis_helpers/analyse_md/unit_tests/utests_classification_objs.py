
import types
import unittest

import gen_basis_helpers.analyse_md.classifier_objs as tCode

class TestByReferenceClassifier(unittest.TestCase):


	def setUp(self):
		self.byRefExecCount = 0
		self.fromRefExecCount = 0
		self.expIndices = [3,4,5]
		self.dudSparseMatrixCalculator = None
		self.createTestObjs()

	def createTestObjs(self):
		self.refClassifier = _StubClassifier(self.expIndices, self.fromRefExecCount)
		self.testObj = tCode.getByReferenceClassifiers([self.refClassifier],startExecCount=self.byRefExecCount)[0]

	def testExpectedCaseA(self):
		#Setup the stored result
		self.refClassifier.classify(self.dudSparseMatrixCalculator)

		#Check our "byReference" gets the correct result
		actIndices = self.testObj.classify(self.dudSparseMatrixCalculator)
		expIndices = self.expIndices
		self.assertEqual( expIndices, actIndices )

	def testRaisesWhenExecutionCountSame(self):
		# Both should now have been executed the same number of times; whereas the reference should always have been executed once more by the time the "byReference" object gets called
		self.byRefExecCount = 1
		self.createTestObjs()
		self.refClassifier.classify(self.dudSparseMatrixCalculator)

		with self.assertRaises(ValueError):
			self.testObj.classify(self.dudSparseMatrixCalculator)

	def testExpectedForMultiCycles(self):
		""" Would fail if execCounts were out of sync """
		expIndices = self.expIndices
		nCycles = 4
		for cycleIdx in range(nCycles):
			self.refClassifier.classify(self.dudSparseMatrixCalculator)
			actIndices = self.testObj.classify(self.dudSparseMatrixCalculator)

		self.assertEqual(expIndices, actIndices)


class _StubClassifier():

	def __init__(self, outputVals, startExecCount):
		self.outputVals = outputVals
		self.execCount = startExecCount

	def classify(self, sparseMatrixCalculator):
		self.execCount += 1
		self.storedClassifyResult = self.outputVals
		return self.outputVals

