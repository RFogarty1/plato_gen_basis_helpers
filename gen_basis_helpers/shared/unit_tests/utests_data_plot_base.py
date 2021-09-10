
import itertools as it
import os
import unittest
import unittest.mock as mock

import numpy as np

import gen_basis_helpers.shared.data_plot_base as tCode

class TestDataPlotterStandardXTickLabels(unittest.TestCase):

	def setUp(self):
		self.xTickLabels = ["labelA","labelB"]
		self.xTickValues = [1,2]
		self.xTickRotation = None
		self.removeXTicks = False
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"removeXTicks":self.removeXTicks,
		             "xTickLabels":self.xTickLabels, "xTickValues": self.xTickValues, "xTickLabelRotation":self.xTickRotation}
		self.testObjA = tCode.DataPlotterStandard(**kwargDict)

	def testXTickLabelsRaisesIfTickPositionsNotPassed(self):
		self.xTickValues = None
		self.createTestObjs()
		mockAxis = mock.Mock()
		with self.assertRaises(ValueError):
			self.testObjA._sortXTickLabels(mockAxis)

	@mock.patch("gen_basis_helpers.shared.data_plot_base.putXTickLabelsOnPlot")
	def testPutXTickLabelsOnPlotCalledWithExpectedArgs(self, mockedPutXTickLabelsOn):
		mockAxis = mock.Mock()
		self.testObjA._sortXTickLabels(mockAxis)
		mockedPutXTickLabelsOn.assert_called_once_with(mockAxis, self.xTickValues, self.xTickLabels, rotation=self.xTickRotation)


class TestSerializePlotters(unittest.TestCase):

	def setUp(self):
		self.xLimA = [1,2]
		self.xLimB = [3,4]
		self.dataA = [ [1,1], [2,4]  ]
		self.dataB = [ [3,9], [4,16] ]
		self.axHandle = None #This COULD be set to something non serializable so... 
		self.tempFilePath = "_data_plotters_dump_test.json"
		self.createTestObjs()

	def tearDown(self):
		os.remove(self.tempFilePath)

	def createTestObjs(self):
		kwargsA = {"xLim":self.xLimA, "data":self.dataA, "axHandle":self.axHandle}
		self.plotterA = tCode.DataPlotterStandard(**kwargsA)
		kwargsB = {"xLim":self.xLimB, "data":self.dataB, "axHandle":self.axHandle}
		self.plotterB = tCode.DataPlotterStandard(**kwargsB)
		self.dataPlotters = [self.plotterA, self.plotterB]
		self._dumpFile()

	def _dumpFile(self):
		tCode.dumpStandardDataPlottersToJson(self.dataPlotters, self.tempFilePath)

	def _readFile(self):
		return tCode.readStandardDataPlottersFromJson(self.tempFilePath)

	#Note that mapPlotDataFunct is not json serialiazble; so that gets tested here too
	def testExpected_noNumpyArrays(self):
		expDataPlotters = self.dataPlotters
		actDataPlotters = self._readFile()

		for expPlotter, actPlotter in it.zip_longest(expDataPlotters, actDataPlotters):
			self._checkExpAndActPlotterEqual(expPlotter, actPlotter)

	def testExpected_numpyArrays(self):
		""" These arent json serializable so need converting to a list first """
		self.dataA = np.array(self.dataA)
		self.dataB = np.array(self.dataB)
		self.createTestObjs()

		expDataPlotters = self.dataPlotters
		actDataPlotters = self._readFile()

		for expPlotter, actPlotter in it.zip_longest(expDataPlotters, actDataPlotters):
			self._checkExpAndActPlotterEqual(expPlotter, actPlotter)

	def _checkExpAndActPlotterEqual(self, expPlotter, actPlotter):
		""" Cant easily implement at general __eq__ so will just use this, knowing which props we actually set """
		expData, actData = np.array(expPlotter.data), np.array(actPlotter.data)
		self.assertTrue( np.allclose(expData,actData) )

		directCmps = ["xLim", "axHandle"]
		for attr in directCmps:
			expAttr,actAttr = getattr(expPlotter,attr), getattr(actPlotter,attr)
			self.assertEqual(expAttr,actAttr)





