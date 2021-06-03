
import copy
import itertools as it
import unittest
import unittest.mock as mock

import gen_basis_helpers.analyse_md.analyse_metadyn_hills as aMetaHillsHelp
import gen_basis_helpers.analyse_md.fit_hills_impl as tCode

class TestGetStartParamsSimple(unittest.TestCase):

	def setUp(self):
		self.inpX =  [1,2,3,4,5]
		self.targY = [5,6,7,6,5]
		self.nGaus = 2
		self.scales = [2,2]

	def _runTestFunct(self):
		return tCode.getStartParamsToFitSimple(self.inpX, self.targY, nGaus=self.nGaus, scales=self.scales)

	def testExpParamsA(self):
		expPos, expHeight = 3, 7/2
		expScale = 2
		expParams = [expPos, expScale, expHeight, expPos, expScale, expHeight]
		actParams = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expParams, actParams)]

class TestGetSimpleFitHillsObjFunct(unittest.TestCase):

	def setUp(self):
		self.inpXVals = [2,3]
		self.targYVals = [6,7]
		self.outYVals = [9,10]
		self.inpParams = [4,5]
		self.expOutput = 4
		self.createTestObjs()

	def createTestObjs(self):
		self.objFunct = mock.Mock()
		self.hillsInfoObj = mock.Mock()
		self.groupedHillsObj = mock.Mock()

		self.objFunct.side_effect = lambda *args, **kwargs: self.expOutput
		self.hillsInfoObj.createGroupedHills.side_effect = lambda *args, **kwargs: self.groupedHillsObj
		self.groupedHillsObj.evalFunctAtVals.side_effect = lambda *args, **kwargs: self.outYVals

	def _runTestFunct(self):
		outFunct = tCode.getFitHillsObjFunctSimple(self.inpXVals, self.targYVals, objFunct=self.objFunct)
		outVal = outFunct(self.inpParams)
		return outVal

	@mock.patch("gen_basis_helpers.analyse_md.fit_hills_impl._getHillsInfoObjFromOneDimParams")
	def testExpectedResultA(self, mockedGetHillsInfoObj):
		#Setup mocks
		mockedGetHillsInfoObj.side_effect = lambda *args,**kwargs: self.hillsInfoObj

		#run Code
		actVal = self._runTestFunct()

		#Check expected calls made + expected value returned 
		mockedGetHillsInfoObj.assert_called_with(self.inpParams)
		self.hillsInfoObj.createGroupedHills.assert_called_with()
		self.groupedHillsObj.evalFunctAtVals.assert_called_with( [[x] for x in self.inpXVals] )
		self.objFunct.assert_called_with( self.targYVals, self.outYVals )
		self.assertEqual(self.expOutput, actVal)


class TestHillsInfoObjToAndFromParams(unittest.TestCase):

	def setUp(self):
		self.scales = [1,2]
		self.heights = [4,5]
		self.positions = [6,7]
		self.createTestObjs()

	def createTestObjs(self):
		self.times = [0 for x in range(len(self.scales))]
		outKwargs = {"times":self.times, "positions":[[x] for x in self.positions],
		             "scales":[[x] for x in self.scales], "heights": [ [x] for x in self.heights]}
		self.hillsObjA = aMetaHillsHelp.MetadynHillsInfo(**outKwargs)
		self.expParams = [ self.positions[0], self.scales[0], self.heights[0],
		                   self.positions[1], self.scales[1], self.heights[1] ]

	def testExpectedParamsFromHillObj(self):
		expParams = self.expParams
		actParams = tCode._getParamsFromOneDimHillsInfoObj( self.hillsObjA )
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expParams,actParams)]

	def testExpectedHillsObjFromParams(self):
		expHills = self.hillsObjA
		actHills = tCode._getHillsInfoObjFromOneDimParams(self.expParams)
		self.assertEqual(expHills,actHills)


class TestPrepareFitDataSimple(unittest.TestCase):

	def setUp(self):
		#Note that rawInpData isnt in ascending order
		self.rawInpData = [ [0,0], [2,4] ,[1,2], [3,9], [4,25] ] 

	def testSetToConstantValueOutsideRange_bothPartsSet(self):
		dataRange = [1.9,3.1]
		expVals = [ [0,4], [2,4], [1,4], [3,9], [4,9] ]
		tCode._setDataToConstantValueOutsideOfRange( self.rawInpData, dataRange )
		actVals = self.rawInpData
		for exp,act in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp, act)]

	def testSetToConstantOutsideRange_onePartSet(self):
		dataRange = [1.9,None]
		expVals = [ [0,4], [2,4], [1,4], [3,9], [4,25] ]
		tCode._setDataToConstantValueOutsideOfRange( self.rawInpData, dataRange )
		actVals = self.rawInpData
		for exp,act in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp, act)]

	def testSetToConstantOutsideRange_neitherSet(self):
		dataRange = [None, None]
		expVals = copy.deepcopy( self.rawInpData )
		tCode._setDataToConstantValueOutsideOfRange( self.rawInpData, dataRange )
		actVals = self.rawInpData
		for exp,act in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp, act)]

	def testGetMirroredData_firstPointNonZero(self):
		self.rawInpData = [ [2,4], [1,2] ]
		expOutData = [ [2,4], [1,2], [3,2] ]
		actOutData = tCode._getPesDataMirroredAroundFirstIdx(self.rawInpData)
		self.assertEqual(expOutData, actOutData)

	#Note this flips the other direction to normal
	def testGetFlippedData(self):
		expVals = [ [0,25], [2,21], [1,23], [3,16], [4,0] ]
		actVals = tCode._getPesDataFlippedWithMinvalZero(self.rawInpData)
		for exp,act in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp, act)]

	def testAllCombinedStandard(self):
		dataRange = [1.9,3.1]
		expValsPreFlip = [ [0,4], [ 2,4], [ 1,4], [ 3,9], [ 4,9], #Pre-mirror values
		                          [-2,4], [-1,4], [-3,9], [-4,9] ] #Post mirror values
		expVals = tCode._getPesDataFlippedWithMinvalZero(expValsPreFlip) #Just messy to figure out...
		actVals = tCode.getStandardPreparedPesDataForFitToFlatten(self.rawInpData, dataRange)
		self.assertEqual(expVals, actVals)


