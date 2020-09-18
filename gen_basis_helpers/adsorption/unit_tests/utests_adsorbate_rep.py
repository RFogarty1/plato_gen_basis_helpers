
import copy
import itertools as it
import types
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.adsorption.adsorbate_rep_objs as tCode

class TestAdsorbateLayerGeoms(unittest.TestCase):

	def setUp(self):
		self.positionsA = [[1,2,3]]
		self.adsorbateGeomA = [ [0.0,1.0,0.0,"Mg"] ]
		self.distanceA = 2.0
		self.surfVector = [0,0,1]
		self.createTestObjs()

	def createTestObjs(self):
		adsorbateA = adsorbateObjFromCartCoords(self.adsorbateGeomA)
		adsorbateListA = [adsorbateA]
		self.testObjA = tCode.AdsorbateLayer(self.positionsA, adsorbateListA, [self.distanceA], self.surfVector)

	def testExpectedGeomNoSymbolsA(self):
		expGeom = [ [1,2+1.0,3 + 0.0 + self.distanceA] ]
		actGeom = self.testObjA.cartCoordsNoAtomSymbols
		self.assertEqual( len(expGeom), len(actGeom) )
		for exp,act in zip(expGeom,actGeom):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp,act)]

	def testExpectedGeomA(self):
		expGeom = [ [1,2+1.0,3 + 0.0 + self.distanceA,"Mg"] ]
		actGeom = self.testObjA.cartCoords
		self.assertEqual( len(expGeom), len(actGeom) )
		for exp,act in zip(expGeom,actGeom):
			[self.assertAlmostEqual(e,a) for e,a in zip(exp[:3],act[:3])]
			self.assertEqual( exp[-1], act[-1] )

class TestAdsorbateObjsEquality(unittest.TestCase):

	def setUp(self):
		self.posA = [0.0, 0.0, 0.0]
		self.posB = [0.0, 0.4, 0.5]
		self.posC = [0.4, 0.0, 0.2]
		self.posD = [0.3, 0.2, 0.6]
		self.eleListA = ["Mg", "Zr", "Ge"]
		self.eleListB = ["Mg", "Zr", "Po"]
		self.createTestObjs()

	def createTestObjs(self):
		self.createCartCoords()
		self.adsorbateA = adsorbateObjFromCartCoords(self.cartCoordsA)
		self.adsorbateB = adsorbateObjFromCartCoords(self.cartCoordsB)
		self.adsorbateC = adsorbateObjFromCartCoords(self.cartCoordsC)

	def createCartCoords(self):
		self.cartCoordsA = self._createCartCoordsA()
		self.cartCoordsB = self._createCartCoordsB()
		self.cartCoordsC = self._createCartCoordsC()

	def _createCartCoordsA(self):
		outCoords = list()
		currPositions = [list(x) for x in [self.posA, self.posB, self.posC]]
		for idx,x in enumerate(currPositions):
			currCoords = list(x) + [self.eleListA[idx]]
			outCoords.append(currCoords)
		return outCoords

	def _createCartCoordsB(self):
		outCoords = list()
		currPositions = [list(x) for x in [self.posA, self.posB, self.posD]]
		for idx,x in enumerate(currPositions):
			currCoords = list(x) + [self.eleListA[idx]]
			outCoords.append(currCoords)
		return outCoords

	def _createCartCoordsC(self):
		outCoords = list()
		currPositions = [list(x) for x in [self.posA, self.posB, self.posC]]
		for idx,x in enumerate(currPositions):
			currCoords = list(x) + [self.eleListB[idx]]
			outCoords.append(currCoords)
		return outCoords

	def testEqualAdsorbateObjsCompareEqual(self):
		copiedAdsorbate = copy.deepcopy(self.adsorbateA)
		objsEqual = tCode.adsorbatesSameWithinError(self.adsorbateA, copiedAdsorbate)
		self.assertTrue(objsEqual)

	def testUnequalAdsorbateObjsCompareUnequal_diffEleLists(self):
		adsA, adsB = self.adsorbateA, self.adsorbateC
		expVal = False
		actVal = tCode.adsorbatesSameWithinError(adsA, adsB)
		self.assertEqual(expVal,actVal)

	def testUnequalAdsorbateObjsCompareUnequal_diffCoords(self):
		adsA, adsB = self.adsorbateA, self.adsorbateB
		expVal = False #Not equal
		actVal = tCode.adsorbatesSameWithinError(adsA, adsB)
		self.assertEqual( expVal,actVal )


class TestAdsorbateLayerEquality(unittest.TestCase):

	def setUp(self):
		self.sitePositionsA = [ [1,2,3], [4,5,6] ]
		self.sitePositionsB = [ [7,8,9], [1,4,7] ]
		self.distsA = [1,2]
		self.distsB = [3,4]
		self.surfVectorA = [-1,-2,-3]
		self.surfVectorB = [-4,-5,-6]
		self.createAdsorbateObjs()
		self.createTestObjs()

	def createAdsorbateObjs(self):
		self.adsorbateA = self._createAdsorbateA()
		self.adsorbateB = self._createAdsorbateB()

	def _createAdsorbateA(self):
		cartCoords = [ [0.0, 0.4, 0.5, "Mg"] ]
		return adsorbateObjFromCartCoords(cartCoords)

	def _createAdsorbateB(self):
		cartCoords = [ [2.0, 4.4, 1.5, "Zr"] ]
		return adsorbateObjFromCartCoords(cartCoords)

	def createTestObjs(self):
		argsA = [ self.sitePositionsA, [self.adsorbateA], self.distsA, self.surfVectorA]
		self.adsorbateLayerA = tCode.AdsorbateLayer(*argsA)

	def testEqualObjsCompareEqual(self):
		copiedObj = copy.deepcopy(self.adsorbateLayerA)
		self.assertEqual(self.adsorbateLayerA, copiedObj)

	def testUnequalIfSitePositionsUnequal(self):
		objA = copy.deepcopy(self.adsorbateLayerA)
		self.sitePositionsA[0][0] += 1
		self.createTestObjs()
		objB = self.adsorbateLayerA
		self.assertNotEqual( objA,objB )

	def testUnequalIfAdsorbateUnequal(self):
		objA = copy.deepcopy(self.adsorbateLayerA)
		self.adsorbateA = self._createAdsorbateB()
		self.createTestObjs()
		objB = self.adsorbateLayerA
		self.assertNotEqual(objA,objB)

	def testUnequalIfDistancesUnequal(self):
		objA = copy.deepcopy(self.adsorbateLayerA)
		self.distsA[0] += 2
		self.createTestObjs()
		objB = self.adsorbateLayerA
		self.assertNotEqual(objA, objB)

	def testUnequalIfSurfVectorsUnequal(self):
		objA = copy.deepcopy(self.adsorbateLayerA)
		self.surfVectorA[0] += 1
		self.createTestObjs()
		objB = self.adsorbateLayerA
		self.assertNotEqual(objA,objB)

def adsorbateObjFromCartCoords(cartCoords):
	return types.SimpleNamespace(geom=cartCoords)


class TestSurfaceWithAdsorbateLayersUsingMocks(unittest.TestCase):


	def setUp(self):
		self.cellLengthA = 3
		self.surfCartCoordsWithSymbols = [ [ 1.5,1.5,1.5,"Mg" ] ]
		self.adsorbateACartCoordsWithSymbols = [ [ 2.0, 2.0, 2.0, "O" ] ] 
		self.createTestObjs()

	def createTestObjs(self):
		self.surfGeomA = self._createCubicSurfaceGeom()
		self.cleanSurfaceStubA = types.SimpleNamespace( unitCell=self.surfGeomA )
		self.adsorbateLayerStubA = types.SimpleNamespace( cartCoords=self.adsorbateACartCoordsWithSymbols )
		self.testObjA = tCode.SurfaceWithAdsorbatesStandard( self.cleanSurfaceStubA, [self.adsorbateLayerStubA] )

	def _createCubicSurfaceGeom(self):
		kwargs = {"lattParams":[self.cellLengthA for x in range(3)],
		          "lattAngles":[90 for x in range(3)]}
		outCell = uCellHelp.UnitCell(**kwargs)
		outCell.cartCoords = self.surfCartCoordsWithSymbols
		return outCell

	def testExpecetedForSingleLayer_fixedAddedVac(self):
		expCartCoords = self.surfCartCoordsWithSymbols + self.adsorbateACartCoordsWithSymbols
		expUCell = self._createCubicSurfaceGeom()
		expUCell.cartCoords = expCartCoords
		actUCell = self.testObjA.unitCell
		self.assertEqual(expUCell,actUCell)

