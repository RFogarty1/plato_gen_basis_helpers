
import math
import unittest
import unittest.mock as mock

import gen_basis_helpers.adsorption.add_adsorbates as tCode

#TODO:Mock copy.deepcopy or mocks get screwed
class TestSingleTypeGetAdsorbates(unittest.TestCase):

	def setUp(self):
		self.adsorbateObjA = mock.Mock()
		self.fractCoverA = 1
		self.nSites = 4
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.SingleTypeGetAdsorbatesForSites(self.adsorbateObjA, self.fractCoverA)
		self.inpSitesA = [mock.Mock() for x in range(self.nSites)]

	def testGetNumberOfAdsorbates_fullCoverage(self):
		expVal = self.nSites
		actVal = self.testObjA.getNumberOfAdsorbateObjs(self.nSites)
		self.assertEqual(expVal, actVal)

	@mock.patch("gen_basis_helpers.adsorption.add_adsorbates.copy.deepcopy")
	def testFullCoverageLeadsToExpectedAdsorbateObjs(self, mockedCopy):
		mockedCopy.side_effect = lambda x:x
		expObjs = [self.adsorbateObjA for x in range(self.nSites)]
		actObjs = self.testObjA.getInputAdsorbateObjs(self.nSites)
		self.assertEqual(expObjs,actObjs)

	def testThrowsValueErrorWhenFractInconsistentWithInpSites(self):
		self.fractCoverA = 0.6
		self.createTestObjs()
		with self.assertRaises(ValueError):
			self.testObjA(self.inpSitesA)

	@mock.patch("gen_basis_helpers.adsorption.add_adsorbates.copy.deepcopy")
	def testExpectedAdsorbateObjsForHalfCoverage(self, mockedCopy):
		mockedCopy.side_effect = lambda x:x
		self.fractCoverA = 0.5
		self.createTestObjs()
		expObjs = [self.adsorbateObjA for x in range(2)] + [None for x in range(2)]
		actObjs = self.testObjA(self.inpSitesA)
		self.assertEqual(expObjs,actObjs)

class TestAddWaterAdsorbatesToBilayer(unittest.TestCase):

	def setUp(self):
		self.waterA, self.waterB = 1, 2
		self.distA, self.distB = 1, 2
		self.createTestObjs()

	def createTestObjs(self):
		argsA = [self.waterA, self.waterB, self.distA, self.distB]
		self.testObjA = tCode.AddWaterAdsorbatesToBilayerSitesStandard(*argsA)

	def testForHexagonOfSitesA(self):
		inpSites = self._loadHexagonOfSites()
		expAdsList = [self.waterA, self.waterB, self.waterB, self.waterA, self.waterA, self.waterB]
		actAdsList = self.testObjA(inpSites)
		self.assertEqual(expAdsList,actAdsList)

	def testDistsForHexagonOfSites(self):
		inpSites = self._loadHexagonOfSites()
		expDists = [self.distA, self.distB, self.distB, self.distA, self.distA, self.distB]
		actDists = self.testObjA.getDistances(inpSites)
		self.assertEqual(expDists, actDists)

	def _loadHexagonOfSites(self):
		outCoords = list()
		dist = 1
		xCompA = dist*math.cos(math.radians(60)) #Two atoms nearest origin have this x-val
		yCompA = math.sqrt( 1-(xCompA**2) ) #two top-most atoms have this y-val

		outCoords.append( [0.0,0.0,0.0] )
		outCoords.append( [xCompA, yCompA, 0.0] )
		outCoords.append( [xCompA, -1*yCompA, 0.0] )
		outCoords.append( [xCompA+dist, yCompA, 0.0] )
		outCoords.append( [xCompA+dist,-1*yCompA, 0.0] )
		outCoords.append( [(2*dist), 0.0, 0.0] )
		return outCoords


