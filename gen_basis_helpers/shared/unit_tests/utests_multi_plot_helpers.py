

import unittest
import unittest.mock as mock

import gen_basis_helpers.shared.multi_plot_helpers as tCode


class TestCreateSimpleMultiPlot(unittest.TestCase):

	def setUp(self):
		self.plotFactoryA = mock.Mock()
		self.plotFactoryB = mock.Mock()
		self.nCols = 4
		self.multiPlotterKwargDict = None
		self.createTestObjs()

	def createTestObjs(self):
		self.plotFactories = [self.plotFactoryA, self.plotFactoryB]

	def runFunct(self):
		return tCode.createSimpleMultiPlot(self.plotFactories, nCols=self.nCols, multiPlotterKwargDict=self.multiPlotterKwargDict)



	@mock.patch("gen_basis_helpers.shared.multi_plot_helpers.mPlotHelp.MultiPlotterStandard")
	@mock.patch("gen_basis_helpers.shared.multi_plot_helpers.multiGrids.RectangularPlotGrid")
	def testSingleColumnPassedWhenOnlyOneFactoryGiven(self, mockedPlotGrid, mockedMultiPlotter):
		self.plotFactories = [mock.Mock()]
		self.nCols = None
		self.runFunct()
		expArgs = ( len(self.plotFactories), len(self.plotFactories) ) #1 plot and therefore 1 column
		actArgs, actKwargs = mockedPlotGrid.call_args
		self.assertEqual(expArgs, actArgs)

	@mock.patch("gen_basis_helpers.shared.multi_plot_helpers.multiGrids.RectangularPlotGrid")
	@mock.patch("gen_basis_helpers.shared.multi_plot_helpers.mPlotHelp.MultiPlotterStandard")
	def testExpArgsToMultiGrid(self, mockedStdPlotter, mockedGridCreator):
		expGrid, expFactories = mock.Mock(), self.plotFactories
		mockedGridCreator.side_effect = lambda *args,**kwargs : expGrid
		expGrid.create.side_effect = lambda *args,**kwargs: mock.Mock(), mock.Mock() 
		self.runFunct()
		mockedStdPlotter.assert_called_once_with(expFactories,expGrid)

	@mock.patch("gen_basis_helpers.shared.multi_plot_helpers.multiGrids.RectangularPlotGrid")
	@mock.patch("gen_basis_helpers.shared.multi_plot_helpers.mPlotHelp.MultiPlotterStandard")
	def testKwargsPassedToMultiPlotter(self, mockedStdPlotter, mockedGridCreator):
		self.multiPlotterKwargDict = {"annotateGraphs":True}
		expKwargs = self.multiPlotterKwargDict
		self.runFunct()
		actArgs,actKwargs = mockedStdPlotter.call_args
		for key in expKwargs:
			self.assertEqual( expKwargs[key], actKwargs[key] )




