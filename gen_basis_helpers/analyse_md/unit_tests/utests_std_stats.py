

import unittest
import unittest.mock as mock

import gen_basis_helpers.analyse_md.std_stats as tCode

class TestBlockingFunction(unittest.TestCase):

	def setUp(self):
		#These inpVals were generated from a random distribution; so really correlation=0 which makes this a kind of awful test
		#But i couldnt be bothered to generate correlated random numbers since its a PITA
		#import random; random.seed(0); [random.randint(-10,10) for x in range(16)]
		self.inpVals = [2, 3, -9, -2, 6, 5, 2, -1, 5, 1, 8, -4, 6, -6, -1, -6]
		self.maxBlockOrder = 2

	def _runTestFunct(self):
		return tCode.getStatsFromBlockingDataUpToMaxOrder(self.inpVals, self.maxBlockOrder)

	def _loadExpAnswerA(self):
		#Here for future reference mainly
		expValsZero = [2, 3, -9, -2, 6, 5, 2, -1, 5, 1, 8, -4, 6, -6, -1, -6]
		expValsOne  = [2.5, -5.5, 5.5, 0.5, 3, 2, 0, -3.5]
		expValsTwo  = [-1.5, 3, 2.5, -1.75]

		#The actual expceted outputs
		orderZeroDict = {"order":0, "mean":0.5625,  "mean_std_dev": 1.248227910546254 , "mean_std_dev_std_dev": 0.22789419450457463}
		orderOneDict  = {"order":1, "mean":0.5625,  "mean_std_dev": 1.2657490018618565, "mean_std_dev_std_dev": 0.33828565018701134}
		orderTwoDict  = {"order":2, "mean":0.5625,  "mean_std_dev": 1.2680981494611   , "mean_std_dev_std_dev": 0.5176989016578825}

		return [orderZeroDict, orderOneDict, orderTwoDict]

	def testExpectedForSimpleRandomCaseA(self):
		expOutput = self._loadExpAnswerA()
		actOutput = self._runTestFunct()
		self.assertTrue( _doExpAndActBlockingDictsMatch(expOutput,actOutput) )

	def testBlockingOfData_oddNumberEntries(self):
		inpData = [1,2,3,4,5]
		expOutput = [1.5, 3.5]
		actOutput = tCode._getDataDividedIntoTwoBlocks(inpData)

		self.assertEqual(len(expOutput), len(actOutput))
		for exp,act in zip(expOutput,actOutput):
			self.assertAlmostEqual(exp,act)

def _doExpAndActBlockingDictsMatch(expDicts, actDicts, eTol=1e-5):
	if len(expDicts) != len(actDicts):
		return False

	numCmpKeys = ["order", "mean", "mean_std_dev", "mean_std_dev_std_dev"]
	for expDict, actDict in zip(expDicts, actDicts):
		for key in numCmpKeys:
			aDiff = abs( expDict[key] - actDict[key] )
			if aDiff > eTol:
				return False

	return True	



class TestGetStandardErrorOfMeanForUncorrelatedData(unittest.TestCase):

	def setUp(self):
		self.nSamples = 100
		self.stdDev = 20

	def _runTestFunct(self):
		return tCode.calcStandardErrorOfMeanForUncorrelatedData(self.stdDev, self.nSamples)

	def testForValA(self):
		expVal = 2
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal)

	def testForValB(self):
		self.nSamples = 100**2
		expVal = 2/10
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal)

