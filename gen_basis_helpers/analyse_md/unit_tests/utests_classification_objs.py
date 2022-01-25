
import copy
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

	def testExpectedMultiCycles_chainOfRefObjs(self):
		""" Want to make sure ref<-byRef<-byRef works """
		expIndices = self.expIndices
		self.testObjB = tCode.getByReferenceClassifiers([self.testObj],startExecCount=self.byRefExecCount)[0]
		nCycles = 4

		for cycleIdx in range(nCycles):
			self.refClassifier.classify(self.dudSparseMatrixCalculator)
			self.testObj.classify(self.dudSparseMatrixCalculator)
			actIndices = self.testObjB.classify(self.dudSparseMatrixCalculator)

		self.assertEqual(expIndices, actIndices)



class _StubClassifier():

	def __init__(self, outputVals, startExecCount):
		self.outputVals = outputVals
		self.execCount = startExecCount

	def classify(self, sparseMatrixCalculator):
		self.execCount += 1
		self.storedClassifyResult = self.outputVals
		return self.outputVals


class TestEqualityForAtomsWithinMinDistRange(unittest.TestCase):


	def setUp(self):
		self.atomIndices = [1,3,4]
		self.distFilterIndices = [ 5,6 ]
		self.distFilterRange = [2.2, 5.6]
		self.minDistVal = -0.05
		self.execCount = 2
		self.createTestObjs()

	def createTestObjs(self):
		currArgs = [self.atomIndices, self.distFilterIndices, self.distFilterRange]
		currKwargs = {"minDistVal":self.minDistVal, "execCount":self.execCount}
		self.testObjA = tCode._AtomsWithinMinDistRangeClassifier(*currArgs, **currKwargs)

	def testCmpEqualA(self):
		objA, objB = self.testObjA, copy.deepcopy(self.testObjA)
		self.assertEqual(objA, objB)

	def testDiffIndicesCmpNotEqual(self):
		objA = copy.deepcopy(self.testObjA)
		self.atomIndices[-1] += 2
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testDiffMinDistCmpNotEqual(self):
		objA = copy.deepcopy(self.testObjA)
		self.minDistVal -= 0.4
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testDiffFilterRangeCmpNotEqual(self):
		objA = copy.deepcopy(self.testObjA)
		self.distFilterRange[-1] += 0.2
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA,objB)

