
import copy
import itertools as it

import unittest
import unittest.mock as mock

import gen_basis_helpers.shared.multi_plotters_base as tCode



class TestMultiPlotterStandard(unittest.TestCase):

	def setUp(self):
		self.plotterA = mock.Mock()
		self.plotterB = mock.Mock()
		self.gridCreatorA = mock.Mock()
		self.annotateStrs = ["a","b"]
		self.annotateGraphs = True
		self.annotatePositionsAbs = [0.05,0.2]
		self.testXlimA, self.testXlimB = [0,2], [0,4]
		self.testYlimA, self.testYlimB = [0,4], [0,6]


		self.createTestObjs()

	def createTestObjs(self):
		#Step 1 create mock axes
		self.mockAxes = [mock.Mock(), mock.Mock()]
		self.gridCreatorA.create.side_effect = lambda *args, **kwargs: (mock.Mock(), self.mockAxes)

		self.mockAxes[0].get_xlim.side_effect = lambda *args,**kwargs: self.testXlimA
		self.mockAxes[1].get_xlim.side_effect = lambda *args,**kwargs: self.testXlimB
		self.mockAxes[0].get_ylim.side_effect = lambda *args,**kwargs: self.testYlimA
		self.mockAxes[1].get_ylim.side_effect = lambda *args,**kwargs: self.testYlimB




		self.plotterIterA = [self.plotterA, self.plotterB]
		self.testObjA = tCode.MultiPlotterStandard(self.plotterIterA, self.gridCreatorA, annotateGraphs=self.annotateGraphs,
		                                           annotateStrs=self.annotateStrs, annotatePositionsAbs=self.annotatePositionsAbs)

	def testStandardMultiPlotterCalledWithOutputAxes(self):
		self.testObjA.create()
		for plotFact, mockAxis in it.zip_longest(self.plotterIterA, self.mockAxes):
			plotFact.createPlot.assert_called_once_with(axHandle=mockAxis)

	def testExpAnnotationStrsPassed(self):
		self.testObjA.create()
		expStrs = self.annotateStrs
		for mockAx, expStr in it.zip_longest(self.mockAxes, expStrs):
			callArgs, callKwargs = mockAx.annotate.call_args
			actStr = callArgs[0]
			self.assertEqual(expStr,actStr)

	def testExpAnnotationStrsPassedUsingInpArg(self):
		self.annotateGraphs=False
		self.createTestObjs()
		self.testObjA.create(annotateGraphs=True)
		expStrs = self.annotateStrs
		for mockAx, expStr in it.zip_longest(self.mockAxes, expStrs):
			callArgs, callKwargs = mockAx.annotate.call_args
			actStr = callArgs[0]
			self.assertEqual(expStr,actStr)

	def testThrowsErrorIfSettingPositionsRelAndAbsAtSameTime(self):
		testPosBoth = [0.1,0.1]
		with self.assertRaises(ValueError):
			tCode.MultiPlotterStandard(self.plotterIterA, self.gridCreatorA, annotatePositionsAbs=testPosBoth, annotatePositionsRel=testPosBoth)

	def testExpPassedWhenUsingRelativePos(self):
		testRelPos = [0.5,0.5]
		self.createTestObjs()

		#Figure out what absolute positions we want passed
		xlimRangeA, xlimRangeB = self.testXlimA[1]-self.testXlimA[0], self.testXlimB[1]-self.testXlimB[0] 
		ylimRangeA, ylimRangeB = self.testYlimA[1]-self.testYlimA[0], self.testYlimB[1]-self.testYlimB[0]
		expAbsPosA = [ self.testXlimA[0] + (xlimRangeA*testRelPos[0]), self.testYlimA[0] + (ylimRangeA*testRelPos[0]) ]
		expAbsPosB = [ self.testXlimB[0] + (xlimRangeB*testRelPos[0]), self.testYlimB[0] + (ylimRangeB*testRelPos[0]) ]
		expAbsPositions = [expAbsPosA, expAbsPosB]

		#Check actual match expected
		self.testObjA.create(annotatePositionsRel=testRelPos)
		for mockAx,expPos in it.zip_longest(self.mockAxes,expAbsPositions):
			callArgs, callKwargs = mockAx.annotate.call_args
			actPos = callArgs[1]
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expPos,actPos)]


	def testGetterAndSetterConsistentForAnnotatePosRel(self):
		self.createTestObjs()
		testRelPos = [0.5,0.5]
		expVal = copy.deepcopy(testRelPos)
		self.testObjA.annotatePositionsRel = testRelPos
		actVal = self.testObjA.annotatePositionsRel
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expVal,actVal)]







