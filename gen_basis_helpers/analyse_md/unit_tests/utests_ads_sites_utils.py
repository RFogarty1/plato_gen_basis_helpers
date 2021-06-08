
import copy
import itertools as it
import math
import unittest
import unittest.mock as mock

import numpy as np

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.ads_sites_impl as adsSiteImplHelp
import gen_basis_helpers.analyse_md.traj_core as trajCoreHelp
import gen_basis_helpers.analyse_md.ads_sites_utils as tCode



class TestGetPlotDataFromAssignedIndices_adsAtomsAdsorbedOverTime(unittest.TestCase):

	def setUp(self):
		self.inpTimes = [1,2]
		self.indicesTimeA = [ [] , [20], [3]  ] #3 adsorption sites, and 4 adsorbates
		self.indicesTimeB = [ [2], [5]  , [20] ]
		self.createTestObjs()

	def createTestObjs(self):
		self.assignedIndices = [self.indicesTimeA, self.indicesTimeB]

	def testMapToAdsorptionCentricA(self):
		expIndicesTimeA = [ [] , [2], [] , [1] ] #Order is by idx. So [2,3,5,20]
		expIndicesTimeB = [ [0], [] , [1], [2] ] 
		expAssignedIndices = [expIndicesTimeA, expIndicesTimeB]
		actAssignedIndices = tCode._getAssignedAdsIndicesInAdsorbingAtomCentricForm(self.assignedIndices)
		self.assertEqual(expAssignedIndices, actAssignedIndices)

	def testExpectedPlotDataReturned(self):
		timeA, timeB = self.inpTimes
		expData = [ [[0, np.nan], [0, timeB]] ,
		            [[1, timeA] , [1, np.nan]],
		            [[2, np.nan], [2, timeB]],
		            [[3, timeA] , [3, timeB]] ]

		actData = tCode.getPlotDataFromAssignedIndices_adsorbatesAdsorbedOverTime(self.inpTimes, self.assignedIndices)
		self._checkExpAndActDataEqual(expData, actData)

	def _checkExpAndActDataEqual(self, expData, actData):
		self.assertEqual( len(expData), len(actData) )
		for expSeries, actSeries in it.zip_longest(expData, actData):
			self.assertEqual( len(expSeries), len(actSeries) )
			for exp, act in it.zip_longest(expSeries, actSeries):
				self.assertAlmostEqual(exp[0],act[0])
				if (exp[1] is np.nan) and (act[1] is np.nan):
					pass
				else:
					self.assertAlmostEqual(exp[1],act[1])



class TestGetPlotDataFromAssignedIndices_adsSiteOccupiedOverTime(unittest.TestCase):

	def setUp(self):
		self.inpTimes = [1,2]
		self.indicesTimeA = [ [] , [2], [] ] 
		self.indicesTimeB = [ [1], [] , [] ]
		self.createTestObjs()

	def createTestObjs(self):
		self.assignedIndices = [self.indicesTimeA, self.indicesTimeB]

	def _runTestFunct(self):
		args = [self.inpTimes, self.assignedIndices]
		return tCode.getPlotDataFromAssignedIndices_adsSiteOccupiedOverTime(*args)

	def testExpectedCaseA(self):
		timeA, timeB = self.inpTimes

		expData = [  [[0,np.nan], [0, timeB]], #0th adsorption site
		             [[1, timeA], [1,np.nan]],
		             [[2,np.nan], [2,np.nan]] ]

		actData = self._runTestFunct()
		self._checkExpAndActDataEqual(expData, actData)


	def testRaisesIfTimesLengthDoesntMatchAssignedLength(self):
		self.inpTimes.append(4)
		with self.assertRaises(AssertionError):
			self._runTestFunct()

	def _checkExpAndActDataEqual(self, expData, actData):
		self.assertEqual( len(expData), len(actData) )
		for expSeries, actSeries in it.zip_longest(expData, actData):
			self.assertEqual( len(expSeries), len(actSeries) )
			for exp, act in it.zip_longest(expSeries, actSeries):
				self.assertAlmostEqual(exp[0],act[0])
				if (exp[1] is np.nan) and (act[1] is np.nan):
					pass
				else:
					self.assertAlmostEqual(exp[1],act[1])


class TestGetAssignedAdsIndicesForTrajectory(unittest.TestCase):

	def setUp(self):
		self.steps = [1,2]
		self.uCells = [mock.Mock(),mock.Mock()]
		self.adsSiteObjs = mock.Mock()
		self.inpIndices = mock.Mock()
		self.maxHozDist = 5
		self.maxTotDist = 10
		self.createTestObjs()

	def createTestObjs(self):
		self.trajStepA = trajCoreHelp.TrajStepBase(unitCell=self.uCells[0], step=self.steps[0])
		self.trajStepB = trajCoreHelp.TrajStepBase(unitCell=self.uCells[1], step=self.steps[1])
		self.trajObjA = trajCoreHelp.TrajectoryInMemory([self.trajStepA, self.trajStepB])

	def _runTestFunct(self):
		args = [self.trajObjA, self.adsSiteObjs, self.inpIndices]
		kwargs = {"maxHozDist":self.maxHozDist, "maxTotDist":self.maxTotDist}
		return tCode.getAssignedAdsIndiceForTrajectory(*args, **kwargs)

	@mock.patch("gen_basis_helpers.analyse_md.ads_sites_utils.assignAdsIndicesToIndividualAdsorptionSites")
	def testExpectedCallsA(self, mockAssignOneGeom):
		expOutput = [mock.Mock(), mock.Mock()]
		mockAssignOneGeom.side_effect = expOutput

		actOutput = self._runTestFunct()

		#Check expected called
		sharedKwargs = {"maxHozDist":self.maxHozDist, "maxTotDist":self.maxTotDist}
		mockAssignOneGeom.assert_any_call(self.uCells[0], self.adsSiteObjs, self.inpIndices, **sharedKwargs)
		mockAssignOneGeom.assert_called_with(self.uCells[1], self.adsSiteObjs, self.inpIndices, **sharedKwargs)

		#Check expected output
		self.assertEqual(expOutput, actOutput)



class TestAssignAdsIndicesToIndividualAdsorptionSites(unittest.TestCase):

	def setUp(self):
		self.adsSiteIndices = [1,2]
		self.adsIndices = [0]
		self.maxHozDist = 4 #If somethings further than this its considered to be "None"
		self.siteName = "site_a"
		self.maxTotDist = None
		self.createTestObjs()

	def createTestObjs(self):
		self.adsObjs = [adsSiteImplHelp.TopStandard(x, siteName=self.siteName) for x in self.adsSiteIndices]
		self.cellA = loadTestCellA()
		cartCoords = [ [9,9,9,"X"],
		               [1,1,1,"Y"],
		               [6,6,5,"Z"] ]
		self.cellA.cartCoords = cartCoords	

	def _runTestFunct(self):
		args = [self.cellA, self.adsObjs, self.adsIndices]
		kwargs = {"maxHozDist":self.maxHozDist, "maxTotDist":self.maxTotDist}
		return tCode.assignAdsIndicesToIndividualAdsorptionSites(*args, **kwargs)

	def testSimpleCaseA(self):
		expOutput = [ [0], [] ]
		actOutput = self._runTestFunct()
		self.assertEqual(expOutput, actOutput)

	#Just flip the indices around; will make sure we dont just use indices within adsIndices rather than cartCoords for example
	def testSimpleCaseB(self):
		self.adsSiteIndices = [2,0]
		self.adsIndices = [1]
		self.createTestObjs()
		expOutput = [ [], [1] ]
		actOutput = self._runTestFunct()
		self.assertEqual(expOutput, actOutput)

	def testExpectedWhenNoAdsorbateInRange(self):
		self.maxHozDist = 1
		expOutput = [ [], [] ]
		actOutput = self._runTestFunct()
		self.assertEqual( expOutput, actOutput )


class TestAssignAdsIndicesToAdsorptionSites(unittest.TestCase):

	def setUp(self):
		self.topIndices = [1,2]
		self.adsIndices = [0]
		self.maxHozDist = 4 #If somethings further than this its considered to be "None"
		self.siteName = "site_a"
		self.maxTotDist = None
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = loadTestCellA()
		self.adsObjs = [adsSiteImplHelp.TopStandard(x, siteName=self.siteName) for x in self.topIndices]

	def _runTestFunct(self):
		args = [self.cellA, self.adsObjs, self.adsIndices]
		kwargs = {"maxHozDist":self.maxHozDist, "maxTotDist":self.maxTotDist}
		return tCode.assignAdsIndicesToAdsorptionSites(*args, **kwargs)

	def testExpectedCaseA(self):
		expDict = {"None":list(), "site_a":[0]}
		actDict = self._runTestFunct()
		self.assertEqual(expDict, actDict)

	def testExpectedCaseB(self):
		self.maxHozDist = 2 #sqrt(8) is the hoz dist for the None site
		self.topIndices = [0]
		self.adsIndices = [1,2]
		self.createTestObjs()
		expDict = {"None":[2], "site_a":[1]}
		actDict = self._runTestFunct()
		self.assertEqual(expDict, actDict)

	def testExpectedTooFarFromHozDist(self):
		self.maxHozDist = 1
		expDict = {"None":[0], "site_a":list()}
		actDict = self._runTestFunct()
		self.assertEqual(expDict, actDict)

	def testExpectedWithMaxTotalDist_closeEnoughToBoth(self):
		self.maxTotDist = 10
		expDict = {"None":list(), "site_a":[0]}
		actDict = self._runTestFunct()
		self.assertEqual(expDict, actDict)

	def testExpectedWithMaxTotalDist_tooFarFromBoth(self):
		self.maxTotDist = math.sqrt(2) #Closest is sqrt(3)
		expDict = {"None":[0], "site_a":list()}
		actDict = self._runTestFunct()
		self.assertEqual(expDict, actDict)

	def testExpectedWhenPbcsMatter(self):
		cartCoords = [ [9,9,9,"X"],
		               [1,1,1,"Y"],
		               [6,6,5,"Z"] ]
		self.cellA.cartCoords = cartCoords
		self.maxHozDist = math.sqrt(8) + 0.1 #Index 1 should work here, if PBCs accounted for
		expDict = {"None":list(), "site_a":[0]}
		actDict = self._runTestFunct()
		self.assertEqual(expDict, actDict)


class TestAddAdsSitesStandard(unittest.TestCase):

	def setUp(self):
		self.topIndices = [1,2] 
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = loadTestCellA()
		self.topSites = [adsSiteImplHelp.TopStandard(idx,siteName="top") for idx in self.topIndices]
		self.testObjA = tCode.AddAdsSitesToGeomsStandard(self.topSites)

	def _runTestFunct(self):
		args = self.cellA
		return self.testObjA.addToUnitCell(self.cellA)

	def testAddToUnitCellSimple(self):
		expCoords = [ [4,4,6,"X"],
		              [5,5,7,"Y"],
		              [6,6,8,"Z"],
		              [5,5,7,"top"],
		              [6,6,8,"top"] ]
		self._runTestFunct()
		actCoords = self.cellA.cartCoords
		testVal = areExpAndActCoordsEqualUsingTemplateCell(expCoords, actCoords, self.cellA)
		self.assertTrue(testVal)


def areExpAndActCoordsEqualUsingTemplateCell(expCoords, actCoords, templCell):
	expCell, actCell = copy.deepcopy(templCell), copy.deepcopy(templCell)
	expCell.cartCoords = expCoords
	actCell.cartCoords = actCoords
	if expCell == actCell:
		return True
	else:
		return False

def loadTestCellA():
	lattParams, lattAngles = [10,10,10], [90,90,90]
	cartCoords = [ [4,4,6,"X"],
	               [5,5,7,"Y"],
	               [6,6,8,"Z"] ]

	outCell = uCellHelp.UnitCell(lattParams=lattParams, lattAngles=lattAngles)
	outCell.cartCoords = cartCoords
	return outCell


