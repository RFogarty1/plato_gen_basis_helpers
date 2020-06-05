
import unittest
import unittest.mock as mock

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

