


import unittest
import unittest.mock as mock

import gen_basis_helpers.shared.multi_axis_plotter as tCode

class TestMultiAxisPlotterStandard(unittest.TestCase):

	def setUp(self):
		self.plotterA = mock.Mock()
		self.plotterB = mock.Mock()
		self.plotterFactories = [self.plotterA, self.plotterB]
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.MultiAxisPlotterStandard(self.plotterFactories)

	def testErrorInCreateIfMoreThanTwoFactoriesGiven(self):
		self.plotterFactories.append( mock.Mock() )
		self.createTestObjs()
		with self.assertRaises(ValueError):
			self.testObjA.create()

	def testExpCallForSingleFactoryCreateNoAxHandleSet(self):
		expOutfig = mock.Mock()
		expAxHandle = mock.Mock()
		self.plotterA.axHandle = None
		self.plotterA.createPlot.side_effect = [expOutfig]
		expOutfig.get_axes.side_effect = lambda : [None] #Only element 0 should be accesed
		actOutAx, actOutfig = self.testObjA._getAxisHandleAndOutFigHandle()
		self.plotterA.createPlot.assert_called_with()
		expOutfig.get_axes.assert_called_once_with() #This is how we get the axis in this case
		self.assertEqual(expOutfig, actOutfig)

	def testExpCallForSingleFactoryWithAxHandleSet(self):
		expOutfig = None
		expAxHandle = mock.Mock()
		self.plotterA.axHandle = expAxHandle
		actAx, actFig = self.testObjA._getAxisHandleAndOutFigHandle()
		self.assertEqual(expOutfig, actFig)
		self.assertEqual(expAxHandle, actAx)

	@mock.patch("gen_basis_helpers.shared.multi_axis_plotter.MultiAxisPlotterStandard._getAxisHandleAndOutFigHandle")
	def testExpCallToFirstPlotFactory(self, mockedGetAxAndFigure):
		expFig, expAx = mock.Mock(), mock.Mock()
		mockedGetAxAndFigure.side_effect = lambda : [expAx, expFig]
		self.testObjA.create()
		self.plotterA.createPlot.assert_called_once_with(axHandle=expAx)

	@mock.patch("gen_basis_helpers.shared.multi_axis_plotter.MultiAxisPlotterStandard._getSecondYIndependentAxis")
	@mock.patch("gen_basis_helpers.shared.multi_axis_plotter.MultiAxisPlotterStandard._getAxisHandleAndOutFigHandle")
	def testExpCallToSecondPlotFactory(self, mockedGetAxAndFigure, mockedGetSecondAx):
		expFig, expAx = mock.Mock(), mock.Mock()
		expSecondAxis = mock.Mock()
		mockedGetSecondAx.side_effect = lambda arg: expSecondAxis
		mockedGetAxAndFigure.side_effect = lambda : [expAx, expFig]
		self.testObjA.create()
		self.plotterB.createPlot.assert_called_once_with(axHandle=expSecondAxis)


