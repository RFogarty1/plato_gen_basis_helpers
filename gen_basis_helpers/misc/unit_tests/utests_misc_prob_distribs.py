

import unittest

import gen_basis_helpers.misc.misc_prob_distribs as tCode

class TestHexagonalNAdsAdjacentDistrib(unittest.TestCase):

	def setUp(self):
		self.nStartSites = 35
		self.nAdded = 3
		self.nPopulateAdjacent = 0
		self.nAdjacentSites = 6

	def _runTestFunct(self):
		currArgs = [self.nStartSites, self.nAdded, self.nPopulateAdjacent]
		return tCode.getHexagonalNAdsAdjacentProbability(*currArgs, nAdjacentSites=self.nAdjacentSites)

	#No combinatorics needed here
	def testZeroAdsCase(self):
		expVal = ( 29/35 )*( 28/34 )*( 27/33 )
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal)

	def testOneAdsCase(self):
		self.nPopulateAdjacent = 1

		#All combos are same probability; but writing it out like this to make it clearer
		expValComboA = (29/35) * (28/34) * (6/33)
		expValComboB = (29/35) * (6/34 ) * (28/33)
		expValComboC = (6/35)  * (29/34) * (28/33)

		expVal = sum( [expValComboA,expValComboC,expValComboC] )
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal,actVal)

	def testTwoAdsCase(self):
		self.nPopulateAdjacent = 2

		expValCombo = (6/35)*(5/34)*(29/33)
		expNCombos = 3
		expVal = expNCombos*expValCombo
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal,actVal)

	def testProbabilitiesAddToOne_5Added(self):
		self.nAdded = 5
		probSum = 0
		for nPopulate in range(self.nAdded+1):
			self.nPopulateAdjacent = nPopulate
			probSum += self._runTestFunct()

		expProbSum = 1
		self.assertAlmostEqual(expProbSum,probSum)
	
	def testZeroProbEdgeCase_nAddedLessThanPopulatedAdjacent(self):
		self.nAdded = 2
		self.nPopulateAdjacent = 5

		expProb = 0
		actProb = self._runTestFunct()
		self.assertAlmostEqual( expProb, actProb )

	def testZeroProbEdgeCase_nAddedMoreThanTotalSites(self):
		self.nAdded = self.nStartSites + 2
		with self.assertRaises(ValueError):
			actProb = self._runTestFunct()

	def testZeroProbEdgeCase_nPopulatedAdjacentGreaterThanNAdjacentSites(self):
		self.nPopulateAdjacent = self.nAdjacentSites+2
		expProb = 0
		actProb = self._runTestFunct()
		self.assertAlmostEqual(expProb,actProb)

