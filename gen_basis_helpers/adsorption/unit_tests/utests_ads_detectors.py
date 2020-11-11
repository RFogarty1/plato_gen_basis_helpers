
import itertools as it
import types
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.adsorption.adsorbate_rep_objs as adsRepObjs
import gen_basis_helpers.adsorption.adsorbate_detectors as tCode


class TestDetectH2OAdsorbates(unittest.TestCase):

	def setUp(self):
		self.minBondLength = 0.1
		self.maxBondLength = 2.5
		self.minAngle = 0
		self.maxAngle = 180

		self.lattParamsA = [10,10,10]
		self.lattAnglesA = [90,90,90]

		#Simple 90 degrees case (since its easier to build @90 degrees)
		self.cartCoordsA = [ [0,0,3,"X"],
		                     [0,5,5,"H"],
		                     [0,5,7,"H"],
		                     [1,2,3,"X"],
		                     [0,4,6,"O"] ]

		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"minBondLength":self.minBondLength, "maxBondLength":self.maxBondLength}
		self.detectorA = tCode.DetectH2OAdsorbatesFromInpGeomStandard(**kwargDict)
		self._createTestCellA()

	def _createTestCellA(self):
		kwargDict = {"lattParams":self.lattParamsA, "lattAngles":self.lattAnglesA}
		self.testCellA = uCellHelp.UnitCell(**kwargDict)
		self.testCellA.cartCoords = self.cartCoordsA

	def testAdsBothInSameCellA_distCriteriaOnly(self):
		expAds = [types.SimpleNamespace(geom=[ [0,4,6,"O"], [0,5,5,"H"], [0,5,7,"H"] ])]
		actAds = self.detectorA( self.testCellA )

		for exp,act in it.zip_longest(expAds,actAds):
			self.assertTrue( adsRepObjs.adsorbatesSameWithinError(exp,act) )

	def testEmptyListWhenBondLengthTooLong(self):
		self.maxBondLength=1
		self.createTestObjs()
		expAds = list()
		actAds = self.detectorA( self.testCellA )
		self.assertEqual(expAds, actAds)

	def testWhenBothHInImageCell(self):
		#Shift to separate H and O into different cells
		self.cartCoordsA[1][1] = 0.9 #H
		self.cartCoordsA[2][1] = 0.9 #H
		self.cartCoordsA[4][1] = 9.9 #O
		self.createTestObjs()

		expAds = [ types.SimpleNamespace(geom=[ [0, 9.9,6,"O"], [0,10.9,5,"H"], [0,10.9,7,"H"] ]),
		           types.SimpleNamespace(geom=[ [0, 0.9,5,"H"], [0, 0.9,7,"H"], [0,-0.1,6,"O"] ]) ]
		actAds = self.detectorA( self.testCellA )

		for exp,act in it.zip_longest(expAds,actAds):
			self.assertTrue( adsRepObjs.adsorbatesSameWithinError(exp,act) )

	def testWhenOneHInImageCell(self):
		self.cartCoordsA[1] = [0,5, 9.5,"H"]
		self.cartCoordsA[4] = [0,4, 0.5,"O"]
		self.cartCoordsA[2] = [0,5, 1.5,"H"]
		self.createTestObjs()

		expAds = [ types.SimpleNamespace(geom=[ self.cartCoordsA[4], self.cartCoordsA[2], [0,5,-0.5,"H"] ]), 
		           types.SimpleNamespace(geom=[ self.cartCoordsA[1], [0,4,10.5,"O"], [0,5,11.5,"H"] ]) ]
		actAds = self.detectorA( self.testCellA )

		for exp,act in it.zip_longest(expAds, actAds):
			self.assertTrue( adsRepObjs.adsorbatesSameWithinError(exp,act) )


class TestDetectH2Adsorbates(unittest.TestCase):

	def setUp(self):
		self.minBondLength = 0.1
		self.maxBondLength = 1.5
		self.minDistOtherNebs = 1.3

		self.lattParamsA = [10,10,10]
		self.lattAnglesA = [90,90,90]

		self.cartCoordsA = [ [0,0,3,"X"],
		                     [0,0,5,"X"],
		                     [5,5,7,"H"],
		                     [5,5,8,"H"] ]

		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"minBondLength":self.minBondLength, "maxBondLength":self.maxBondLength,
		             "minDistOtherNebs":self.minDistOtherNebs}
		self.detectorA = tCode.DetectH2AdsorbatesFromInpGeomStandard(**kwargDict)
		self._createTestCellA()

	def _createTestCellA(self):
		kwargDict = {"lattParams":self.lattParamsA, "lattAngles":self.lattAnglesA}
		self.testCellA = uCellHelp.UnitCell(**kwargDict)
		self.testCellA.cartCoords = self.cartCoordsA


	def testAdsFoundWhenBothInSameCellA(self):
		expAds = [ types.SimpleNamespace(geom =[[5,5,7,"H"], [5,5,8,"H"]]) ]
		actAds = self.detectorA( self.testCellA )
		for exp,act in it.zip_longest(expAds,actAds):
			adsRepObjs.adsorbatesSameWithinError(exp,act)

	def testAdsFoundWhenAtomsInNeighbouringCells(self):
		self.cartCoordsA[2] = [5,5,0.5,"H"]
		self.cartCoordsA[3] = [5,5,9.5,"H"]
		self.createTestObjs()
		expAds = [ types.SimpleNamespace(geom=[ [5,5,0.5,"H"], [5,5,-0.5,"H"] ]),
		           types.SimpleNamespace(geom=[ [5,5,9.5,"H"], [5,5,10.5,"H"] ]) ]
		actAds = self.detectorA( self.testCellA )
		for exp,act in it.zip_longest(expAds, actAds):
			adsRepObjs.adsorbatesSameWithinError(exp,act)

	def testH3NotDetected(self):
		self.cartCoordsA.append( [5,5,6,"H"] ) 
		self.createTestObjs()
		expAds = list()
		actAds = self.detectorA( self.testCellA )
		self.assertEqual(expAds,actAds)

	def testH3NotDetectedForImageH(self):
		self.cartCoordsA[2] = [5,5,0.5,"H"]
		self.cartCoordsA[3] = [5,5,9.5,"H"]
		self.cartCoordsA.append( [5,5,8.5,"H"] )
		self.createTestObjs()
		expAds = list()
		actAds = self.detectorA( self.testCellA )
		self.assertEqual(expAds,actAds)


