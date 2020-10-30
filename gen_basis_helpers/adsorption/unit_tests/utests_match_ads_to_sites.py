
import itertools as it
import types
import unittest
import unittest.mock as mock

import gen_basis_helpers.adsorption.adsorbate_rep_objs as adsRepObjs
import gen_basis_helpers.adsorption.match_ads_to_sites as tCode


class TestH2GeomEnforcerFilterFunct(unittest.TestCase):

	def setUp(self):
		self.maxDist = 4
		self.surfOutwardsVect = [1,1,1]
		self.sitePos = [1,1,1]
		self.adsObj = types.SimpleNamespace( geom=[[0,0,0,"H"],[0,0,1,"H"]] )
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.AdsObjIsH2FilterFunct(maxDist=self.maxDist)

	def _runTestFunctStd(self):
		return self.testObjA(self.sitePos, self.surfOutwardsVect, self.adsObj)

	def testTrueForH2WithShorterBondThanMaxDist(self):
		expVal = True
		actVal = self._runTestFunctStd()
		self.assertEqual(expVal,actVal)

	def testFalseForH3(self):
		currGeom = [ [0,0,0,"H"], [0,0,1,"H"], [0,0,0.5,"H"] ]
		self.adsObj = types.SimpleNamespace( geom=currGeom )
		expVal = False
		actVal = self._runTestFunctStd()
		self.assertEqual(expVal, actVal)

	def testFalseForHF(self):
		currGeom = [ [0,0,0,"H"], [0,0,1,"F"] ]
		self.adsObj = types.SimpleNamespace( geom=currGeom )
		expVal = False
		actVal = self._runTestFunctStd()
		self.assertEqual(expVal,actVal)

	def testFalseForH2WithBondLengthGreaterThanMax(self):
		currGeom = [ [0,0,0,"H"], [0,0,self.maxDist+1,"H"] ]
		self.adsObj = types.SimpleNamespace( geom=currGeom )
		expVal = False
		actVal = self._runTestFunctStd()
		self.assertEqual(expVal,actVal)


class TestAdsObjCentreWithinMaxDistOfSiteFilterFunct(unittest.TestCase):

	def setUp(self):
		self.maxDistA = 2
		self.surfOutwardsVect = [0,0,1]
		self.sitePos = [0,0,0]
		self.geomA = [ [0,0,1,"H"], [0,0,2.2,"H"] ]
		self.createTestObjs()

	def createTestObjs(self):
		self.adsObj = types.SimpleNamespace( geom=self.geomA )
		self.testObjA = tCode.AdsObjCentroidDistFromSiteFilterFunct(self.maxDistA)

	def _runTestFunct(self):
		return self.testObjA(self.sitePos, self.surfOutwardsVect, self.adsObj)

	def testTrueWhenCentreWithinDist(self):
		expVal = True
		actVal = self._runTestFunct()	
		self.assertEqual(expVal,actVal)

	def testFalseWhenCentreBeyondDist(self):
		self.maxDistA = 1
		self.createTestObjs()
		expVal = False
		actVal = self._runTestFunct()
		self.assertEqual(expVal,actVal)


class TestAdsObjCentroidWithinHorizontalDist(unittest.TestCase):

	def setUp(self):
		self.maxDistA = 1.05
		self.surfOutwardsVect = [0,0,1]
		self.sitePos = [0,0,0]
		self.geomA = [ [0.9,0,2,"H"], [1.1,0,2,"H"] ]
		self.createTestObjs()

	def createTestObjs(self):
		self.adsObj = types.SimpleNamespace( geom=self.geomA )
		self.testObjA = tCode.AdsObjCentroidHorizontalDistFromSiteFilterFunct(self.maxDistA)

	def _runTestFunct(self):
		return self.testObjA(self.sitePos, self.surfOutwardsVect, self.adsObj)

	def testTrueWhenCentreWithinDist(self):
		expVal = True
		actVal = self._runTestFunct()
		self.assertEqual(expVal,actVal)		

	def testTrueWhenCentroidWithinDistWithFlippedSurfVector(self):
		self.surfOutwardsVect = [-1*x for x in self.surfOutwardsVect]
		expVal = True
		actVal = self._runTestFunct()
		self.assertEqual(expVal, actVal)

	def testFalseWhenCentreBeyondDist(self):
		self.maxDistA = 0.99
		self.createTestObjs()
		expVal = False
		actVal = self._runTestFunct()
		self.assertEqual(expVal, actVal)


class TestMatchSitesToH2AdsObjs(unittest.TestCase):

	def setUp(self):
		self.maxDistTotal = 2
		self.maxDistHoz = 1
		self.maxHH = 2
		self.surfOutwardsVector = [0,0,1]
		self.geomA = [ [0,0,0,"H"], [0,0,1,"H"] ]
		self.sitesA = [ [0,0,0] ]
		self.createTestObjs()

	def createTestObjs(self):
		self.adsObjsA = [ types.SimpleNamespace(geom=self.geomA) ]
		kwargDict = {"maxDistHoz":self.maxDistHoz, "maxHHBond":self.maxHH}
		self.testObjA = tCode.MatchSitesToH2AdsObjsStandard(self.maxDistTotal, **kwargDict)

	def _runTestFunct(self):
		return self.testObjA(self.sitesA, self.surfOutwardsVector, self.adsObjsA)

	def testExpectedObjReturnedWhenWithinRange(self):
		expObjs = [ [self.adsObjsA[0]] ]
		actObjs = self._runTestFunct()
		self.checkExpListMatchesActual(expObjs, actObjs)

	def testEmptyListReturnedWhenOutsideHozRange(self):
		self.geomA[1][0] +=2.2
		expObjs = [ list() ]
		actObjs = self._runTestFunct()
		self.checkExpListMatchesActual(expObjs, actObjs)

	def checkExpListMatchesActual(self,expLists,actLists):
		for expList,actList in it.zip_longest(expLists,actLists):
			for e,a in it.zip_longest(expList,actList):
				self.assertTrue( adsRepObjs.adsorbatesSameWithinError(e,a) )

