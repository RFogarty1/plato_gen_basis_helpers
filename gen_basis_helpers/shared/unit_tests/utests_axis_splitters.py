
import itertools as it

import unittest
import unittest.mock as mock

import gen_basis_helpers.shared.axis_splitters as tCode


class TestSingleAxisSplitterTemplate(unittest.TestCase):

	def setUp(self):
		self.positions = [ [0.1,0.4] ] #The FRACTIONAL [start,end] positions of the new axes
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.SingleAxisSplitterTemplate.__new__(tCode.SingleAxisSplitterTemplate)
		self.testObjA.positions = self.positions
		self.testObjA.splitAxes = False

	#Testing that we convert the fractional values of positions into actual values
	@mock.patch("gen_basis_helpers.shared.axis_splitters.SingleAxisSplitterTemplate._getOrigAxisReducedPositionsFromInpHandle")
	def testGetNewAxisPositions(self, mockedGetOrigAxisPosition):
		fakeStartPos, fakeLength = 1.0, 10.0
		fakePositions = [fakeStartPos, fakeLength]
		mockedAxHandle = mock.Mock()
		mockedGetOrigAxisPosition.side_effect = lambda inpHandle: fakePositions

		#Figure out the expected positions
		expStartPosA = fakeStartPos + (fakeLength*self.positions[0][0])
		expLengthA = fakeLength * (self.positions[0][1] - self.positions[0][0])
		expOutPositions = [ [expStartPosA, expLengthA] ]
		actOutPositions = self.testObjA._getNewAxesPositions(mock.Mock())

		for exp,act in it.zip_longest(expOutPositions, actOutPositions):
			self.assertAlmostEqual(exp[0],act[0])
			self.assertAlmostEqual(exp[1],act[1])

	def testGetParentFigure(self):
		expFigureHandle = mock.Mock()
		mockAxHandle = mock.Mock()
		mockAxHandle.figure = expFigureHandle
		actFigureHandle = self.testObjA._getFigureHandleFromAxHandle(mockAxHandle)
		self.assertEqual(expFigureHandle, actFigureHandle)

	@mock.patch("gen_basis_helpers.shared.axis_splitters.SingleAxisSplitterTemplate._getFullNewAxesPositionsFromOrigAx")
	@mock.patch("gen_basis_helpers.shared.axis_splitters.SingleAxisSplitterTemplate._getFigureHandleFromAxHandle")
	def testCreateNewAxesFromOrigHandle(self, mockGetFigHandle, mockGetAxPositions):
		self.positions = [ [0.1,0.3], [0.5,0.8] ]
		self.fakePositions = [2,3]
		self.createTestObjs()
		mockInpAx, mockFigureHandle = mock.Mock(), mock.Mock()
		expNewAxes = [mock.Mock(), mock.Mock()]

		mockGetFigHandle.side_effect = lambda inpAx: mockFigureHandle
		mockGetAxPositions.side_effect = lambda *args: self.fakePositions
		mockFigureHandle.add_axes.side_effect = expNewAxes
		actNewAxes = self.testObjA._createNewAxesFromOrigHandle(mockInpAx)

		mockGetFigHandle.assert_called_once_with(mockInpAx)
		mockFigureHandle.add_axes.assert_any_call(self.fakePositions[0])
		mockFigureHandle.add_axes.assert_any_call(self.fakePositions[1])
		self.assertTrue(self.testObjA.splitAxes)
		self.assertEqual(expNewAxes,actNewAxes)

	@mock.patch("gen_basis_helpers.shared.axis_splitters.SingleAxisSplitterTemplate._getSharedKwargDictForNewAxes")
	@mock.patch("gen_basis_helpers.shared.axis_splitters.copy.deepcopy")
	def testPlotterInstanceCalledOnNewAxes(self, mockCopy, mockGetNewAxisDict):
		expNewAxDict = dict()
		mockCopy.side_effect = lambda x:x
		mockGetNewAxisDict.side_effect = lambda : expNewAxDict
		plotterInstance, mockAxis = mock.Mock(), mock.Mock()
		mockAxes = [mockAxis]
		self.testObjA._addPlotToAxes(mockAxes, plotterInstance)
		plotterInstance.createPlot.assert_any_call( axHandle=mockAxis, **expNewAxDict ) 


