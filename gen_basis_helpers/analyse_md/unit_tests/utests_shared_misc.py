
import unittest
import unittest.mock as mock

import gen_basis_helpers.analyse_md.shared_misc as tCode

class TestGetSlicesForMerging(unittest.TestCase):

	def setUp(self):
		self.stepIndicesA = [0,6]
		self.stepIndicesB = [8,16]
		self.nStepsA = 4
		self.nStepsB = 5
		self.overlapStrat = None
		self.createTestObjs()

	def createTestObjs(self):
		self.stepIndices = [self.stepIndicesA, self.stepIndicesB]
		self.nSteps = [self.nStepsA, self.nStepsB]

	def _runTestFunct(self):
		return tCode.getSlicesForMergingTrajectories(self.stepIndices, self.nSteps, overlapStrat=self.overlapStrat)

	def testExpectedForSimpleCaseA(self):
		expIndices = [ [0,4], [0,5] ]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices, actIndices)

	def testRaisesErrorForOverlapStratNone(self):
		self.stepIndicesB[0] = self.stepIndicesA[1]
		with self.assertRaises(ValueError):
			self._runTestFunct()

	def testStratSimpleExpectedOutput(self):
		self.stepIndicesB[0] = self.stepIndicesA[1]
		self.overlapStrat = "simple"
		expIndices = [ [0,3], [0,5] ]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices,actIndices)

