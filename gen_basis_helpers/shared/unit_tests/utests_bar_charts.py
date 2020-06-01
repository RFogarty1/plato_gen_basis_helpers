

import itertools as it
import numpy as np

import unittest
import unittest.mock as mock

import gen_basis_helpers.shared.bar_chart_plotter as tCode

class TestBarChartPlotter(unittest.TestCase):

	def setUp(self):
		self.yVals = [ [1,2], [3,4] ]
		self.labels = ["labelA", "labelB"]
		self.fixDataSeriesOrLabels = "labels"
		self.barWidths = 2
		self.barSeparationWithinLabels = 1
		self.barSeparationBetweenLabels = 3
		self.xStartPos = 1
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"data":self.yVals, "fixDataSeriesOrLabels":self.fixDataSeriesOrLabels,
		             "barWidths":self.barWidths, "xStartPos":self.xStartPos,
		             "barSeparationWithinLabels":self.barSeparationWithinLabels,
		             "barSeparationBetweenLabels":self.barSeparationBetweenLabels}
		self.testObjA = tCode.LabelledBarChartPlotterStandard(**kwargDict)

	def testDataSeriesOrLabelsSetterRaisesForInvalidVal(self):
		testVal = "fake_invalid_value"
		with self.assertRaises(AttributeError):
			self.testObjA.fixDataSeriesOrLabels = testVal

	#This is the more "Native" format of data; all data for one label is together 
	def testOutputDataForFixedLabelDataCase(self):
		self.fixDataSeriesOrLabels = "labels"
		self.createTestObjs()
		expVals = self.yVals
		actVals = self.testObjA._getYValsInFixedLabelFormat()
		self.assertEqual(expVals,actVals)


	def testOutputDataForFixedDataSeriesCase(self):
		self.yVals = [ [1,2,3], [4,5,6] ] #3 labels; 2 data series
		self.fixDataSeriesOrLabels = "dataSeries"
		self.createTestObjs()
		expVals = [ [1,4], [2,5], [3,6] ]
		actVals = self.testObjA._getYValsInFixedLabelFormat()
		self.assertEqual(expVals,actVals)


	def testGetXValsForExampleSpacings(self):
		startXPos = self.xStartPos+(self.barWidths*0.5)
		expXVals = [ [startXPos, startXPos + self.barSeparationWithinLabels + self.barWidths],
		             [10, 13] ]
		actXVals = self.testObjA._getXVals()

		for expList,actList in it.zip_longest(expXVals, actXVals):
			for exp,act in it.zip_longest(expList,actList):
				self.assertAlmostEqual( exp,act )

	@mock.patch("gen_basis_helpers.shared.bar_chart_plotter.LabelledBarChartPlotterStandard._getYValsInFixedLabelFormat")
	@mock.patch("gen_basis_helpers.shared.bar_chart_plotter.LabelledBarChartPlotterStandard._getXVals")
	def testGetXVsYForDataSeries(self, xValGetter, yValGetter):
		""" Test that we can combine the label-centric x vs y data into one thats data-series centric (we plot 1 data series at a time) """
		startXVals = [ [1,2,3], [5, 6 ,7 ] ] #Each list is x-values for one LABEL (e.g. a basis set)
		startYVals = [ [7,8,9], [10,11,12] ] #Each list is y-values for one LABEL (e.g. [eqmVol, surfaceEnergy] for a basis set)
		xValGetter.side_effect = lambda *args:startXVals
		yValGetter.side_effect = lambda *args:startYVals

		expDataSeriesA = [ (1,7), (5,10) ]
		expDataSeriesB = [ (2,8), (6,11) ]
		expDataSeriesC = [ (3,9), (7,12) ]

		expXYVals = [ np.array(expDataSeriesA), np.array(expDataSeriesB), np.array(expDataSeriesC) ] 
		actXYVals = self.testObjA._getToPlotData()

		#Actually dont need to convert anything to a np array to pass this test
		for exp,act in it.zip_longest(expXYVals, actXYVals):
			self.assertTrue( np.allclose(exp,act) )


class TestGetXTickPositionsForLabels(unittest.TestCase):

	def setUp(self):
		self.testDataA = [ [1,2], [3,4] ]
		self.testDataB = [ [5,5], [7,7] ]
		self.testDataC = [ [10,1], [11,1] ]
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = [ np.array(x) for x in [self.testDataA,self.testDataB,self.testDataC] ]

	def runTestFunct(self, vals):
		return tCode._getXTickPositionsForLabels(vals)

	def testExpectedTicksGivenForTestData(self):
		expTickPositions = [2, 6, 10.5] #Halfway between the values
		actTickPositions = self.runTestFunct(self.testObjA)
		for exp,act in it.zip_longest(expTickPositions,actTickPositions):
			self.assertAlmostEqual(exp,act)


