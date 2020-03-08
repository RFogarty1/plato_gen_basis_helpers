

import unittest
import unittest.mock as mock

import gen_basis_helpers.shared.multi_plot_grids as tCode

class TestRectangularGrid(unittest.TestCase):

	def setUp(self):
		self.nCols = 4
		self.nPlots = 9
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.RectangularPlotGrid(self.nCols,self.nPlots)

	def testCorrectNumbRowsAndCols(self):
		expRowsAndCols = (3, self.nCols)
		actRowsAndCols = self.testObjA.dims
		self.assertEqual(expRowsAndCols,actRowsAndCols)	

	@mock.patch("gen_basis_helpers.shared.multi_plot_grids.RectangularPlotGrid._addSubPlotsToFigure")
	@mock.patch("gen_basis_helpers.shared.multi_plot_grids.plt")
	def testExpectedCallsToPyPlot(self, mockedPlt, mockedAddSubplots):
		retFigMock = mock.Mock()
		mockedPlt.figure.side_effect = lambda *args,**kwargs: retFigMock
		self.testObjA.create()
		mockedPlt.figure.assert_called_once_with(constrained_layout=True)
		retFigMock.add_gridspec.assert_called_once_with(*self.testObjA.dims)

