
import itertools as it
import unittest
import unittest.mock as mock

import gen_basis_helpers.adsorption.site_occupiers as tCode


class TestOccupyInOrder(unittest.TestCase):

	def setUp(self):
		self._loadDefaultFourSitePositions()
		self.createTestObjs()

	def _loadDefaultFourSitePositions(self):
		sitePositions = _returnDefaultFourSitePositions()
		self.sitePosA, self.sitePosB, self.sitePosC, self.sitePosD = sitePositions

	def createTestObjs(self):
		self.allSitePosA = [self.sitePosA, self.sitePosB, self.sitePosC, self.sitePosD]
		self.testObjA = tCode.OccupyInOrderSiteOccupier()

	def testFullOccupancy(self):
		testAdsorbates = [mock.Mock() for x in self.allSitePosA]
		expAdsorbates = testAdsorbates
		actAdsorbates = self.testObjA(self.allSitePosA, testAdsorbates)
		self.assertEqual(expAdsorbates,actAdsorbates)

	def testNoOccupancy(self):
		testAdsorbates = list()
		expAdsorbates = [None for x in self.allSitePosA]
		actAdsorbates = self.testObjA(self.allSitePosA, testAdsorbates)
		self.assertEqual(expAdsorbates, actAdsorbates)

	def testHalfOccupancy(self):
		testAdsorbates = [mock.Mock(), mock.Mock()]
		expAdsorbates = testAdsorbates + [None,None]
		actAdsorbates = self.testObjA(self.allSitePosA, testAdsorbates)
		self.assertEqual(expAdsorbates, actAdsorbates)

class TestOccupyClosestIgnoringGeoms(unittest.TestCase):

	def setUp(self):
		self._loadDefaultFourSitePositions()
		self.createTestObjs()

	def _loadDefaultFourSitePositions(self):
		sitePositions = _returnDefaultFourSitePositions()
		self.sitePosA, self.sitePosB, self.sitePosC, self.sitePosD = sitePositions

	def createTestObjs(self):
		self.allSitePosA = [self.sitePosA, self.sitePosB, self.sitePosC, self.sitePosD]
		self.testObjA = tCode.OccupyClosestIgnoringAdsorbateGeomsSiteOccupier()
	
	def testNoCoverageCase(self):
		testAdsorbates = list() #signifies no adsorbates to be added
		expOutput = [None for x in self.allSitePosA]
		actOutput = self.testObjA(self.allSitePosA, testAdsorbates)
		self.assertEqual( len(expOutput), len(actOutput) )
		for exp,act in it.zip_longest(expOutput, actOutput):
			self.assertEqual(exp,act)

	def testHalfCoverageCase(self):
		adsA,adsB = mock.Mock(), mock.Mock()
		testAdsorbates = [adsA, adsB]
		expOutput = [adsA, None, None, adsB]
		actOutput = self.testObjA(self.allSitePosA, testAdsorbates)
		self.assertEqual( len(expOutput),len(actOutput) )

		for exp,act in zip(expOutput, actOutput):
			self.assertEqual(exp,act)

	def testFullCoverageCase(self):
		testAdsorbates = [mock.Mock() for x in self.allSitePosA]
		expOutput = [ testAdsorbates[0], testAdsorbates[3] , testAdsorbates[2], testAdsorbates[1]]
		actOutput = self.testObjA(self.allSitePosA, testAdsorbates)
		for exp,act in it.zip_longest(expOutput, actOutput):
			self.assertEqual(exp,act)

def _returnDefaultFourSitePositions():
	siteA = [0.0, 0.0, 0.0]
	siteB = [1.1, 0.0, 0.0]
	siteC = [0.0, 1.1, 0.0]
	siteD = [0.0, 0.5, 0.0]
	return [siteA, siteB, siteC, siteD]



class TestMinimizeAverageDistToOccupiedSiteOccupier(unittest.TestCase):

	def setUp(self):
		self.createTestObjs()

	def createTestObjs(self):
		self.sitePosA = self._loadFlat2DimGrid()
		self.testObjA = tCode.MinimizeAverageOfNonPeriodicOccSiteDistsStandardSiteOccupier()

	def _loadFlat2DimGrid(self):
		allXVals = [0,1,2, 0,1,2, 0,1,2]
		allYVals = [0,0,0, 1,1,1, 2,2,2]
		allZ = 0.0
		allVals = list()
		for xVal,yVal in it.zip_longest(allXVals, allYVals):
			currVals = [xVal,yVal,allZ]
			allVals.append(currVals)
		return allVals

	def testForHalfishOccupancyA(self):
		adsA, adsB, adsC, adsD = [mock.Mock() for x in range(4)]
		testAdsorbates = [adsA, adsB, adsC, adsD]
		expOccIndices = [0, 1, 3, 4]
		expOutput = [None for x in range(len(self.sitePosA))]
		expOutput[0], expOutput[1], expOutput[3], expOutput[4] = adsA, adsB, adsC, adsD
		actOutput = self.testObjA(self.sitePosA, testAdsorbates)
		self.assertEqual(expOutput, actOutput)

