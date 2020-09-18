
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


