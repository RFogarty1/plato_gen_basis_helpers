

import itertools as it

import unittest
import unittest.mock as mock

import gen_basis_helpers.shared.multi_plotters_base as tCode



class TestMultiPlotterStandard(unittest.TestCase):

	def setUp(self):
		self.plotterA = mock.Mock()
		self.plotterB = mock.Mock()
		self.gridCreatorA = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		self.plotterIterA = [self.plotterA, self.plotterB]
		self.testObjA = tCode.MultiPlotterStandard(self.plotterIterA, self.gridCreatorA)

	def testStandardMultiPlotterCalledWithOutputAxes(self):
		mockAxes =  [mock.Mock(), mock.Mock()]
		self.gridCreatorA.create.side_effect = lambda *args, **kwargs: (mock.Mock(), mockAxes)
		self.testObjA.create()
		for plotFact, mockAxis in it.zip_longest(self.plotterIterA, mockAxes):
			plotFact.createPlot.assert_called_once_with(axHandle=mockAxis)



