
import copy
import itertools as it
import types
import unittest
import unittest.mock as mock

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



class TestNonHyChainedClassifierAllCommon(unittest.TestCase):

	def setUp(self):
		self.nonHyA = [ [3], [6,5], [4] ]
		self.nonHyB = [ [5,6], [3] ]

		self.hyA = [ [10]   , [11,12], [13] ]
		self.hyB = [ [11,12], [10] ]

		self.dudSparseMatrixObj = 4 #Just want something that normal equalities work with

		self.createTestObjs()

	def createTestObjs(self):
		self.classA, self.classB = mock.Mock(), mock.Mock()
		self.classA.classify.side_effect = lambda x: (self.nonHyA, self.hyA)
		self.classB.classify.side_effect = lambda x: (self.nonHyB, self.hyB)

		self.testObj = tCode._NonHyAndHyChainedClassifier_allCommon([self.classA, self.classB])

	def _checkExpAndActEqual(self, expNonHy, expHy, actNonHy, actHy):
		for expHy,actHy in it.zip_longest(expHy, actHy):
			self.assertEqual( list(expHy), list(actHy) )

		for expNonHy, actNonHy in it.zip_longest(expNonHy, actNonHy):
			self.assertEqual( list(expNonHy), list(actNonHy) )

	def testExpectedCaseA(self):
		expNonHy = [ [3] , [5,6]   ]
		expHy =    [ [10], [11,12] ]

		actNonHy, actHy = self.testObj.classify(self.dudSparseMatrixObj)
		self.classA.classify.assert_called_with(self.dudSparseMatrixObj)
		self.classB.classify.assert_called_with(self.dudSparseMatrixObj)

		self._checkExpAndActEqual(expNonHy, expHy, actNonHy, actHy)

	def testExpected_noHyIndices(self):
		self.hyA = [ list(), list(), list() ]
		self.hyB = [ list(), list() ]
		self.createTestObjs()
		expNonHy = [ [3], [5,6] ]
		expHy = [ list(), list() ]

		actNonHy, actHy = self.testObj.classify(self.dudSparseMatrixObj)
		self._checkExpAndActEqual(expNonHy, expHy, actNonHy, actHy)

	def testExpected_noNonHyIndices(self):
		self.nonHyA, self.nonHyB = [ list(), list(), list()], [list(), list()]
		self.createTestObjs()
		expNonHy, expHy = [list(), list()], [ [10], [11,12] ]

		actNonHy, actHy = self.testObj.classify(self.dudSparseMatrixObj)
		self._checkExpAndActEqual(expNonHy, expHy, actNonHy, actHy)

	def testExpected_reversedOrderIndices(self):
		self.nonHyA, self.nonHyB = [x for x in reversed(self.nonHyA)], [x for x in reversed(self.nonHyB)]
		self.hyA, self.hyB = [x for x in reversed(self.hyA)], [x for x in reversed(self.hyB)]

		expNonHy = [ [5,6], [3] ]
		expHy = [ [11,12], [10] ]  
		actNonHy, actHy = self.testObj.classify(self.dudSparseMatrixObj)
		self._checkExpAndActEqual(expNonHy, expHy, actNonHy, actHy)




