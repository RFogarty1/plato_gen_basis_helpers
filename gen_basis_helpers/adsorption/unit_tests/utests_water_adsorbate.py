
import itertools as it
import math
import unittest
import unittest.mock as mock

import gen_basis_helpers.adsorption.water_adsorbate as tCode

class TestStandardRefPosGetter(unittest.TestCase):

	def setUp(self):
		self.ohDists = [1,2]
		self.angle = 90 #Allows us to have 45 degrees around the bisect; which makes math easier
		self.createTestObjs()

	def createTestObjs(self):
		self.testObj = tCode.WaterRefPosGetterStandard()

	def testVsExpA(self):
		expGeom = self._loadExpGeomA()
		actGeom = self.testObj(self.ohDists, self.angle)
		self._checkGeomsEqual(expGeom,actGeom)

	def _loadExpGeomA(self):
		contribA = (self.ohDists[0]) / math.sqrt(2)
		contribB = (self.ohDists[1]) / math.sqrt(2) #MAAAAAAAAYBE

		outGeom = [ [ 0.0, 0.0, 0.0, "O"],
		            [ contribA, contribA, 0.0, "H"],
		            [ contribB, -1*contribB, 0.0, "H"] ]
		return outGeom

	def _checkGeomsEqual(self, expGeom, actGeom):
		for exp,act in it.zip_longest(expGeom, actGeom):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp[:3],act[:3])]
			self.assertEqual(exp[-1],act[-1])

class TestGeometryFromWaterAdsorbate(unittest.TestCase):

	def setUp(self):
		self.ohDists = [1,2]
		self.angle = 90
		self.pitch = 0
		self.roll = 0
		self.azimuthal = 0
		self.refPosGetter = tCode.WaterRefPosGetterStandard()
		self.createTestObjs()

	def createTestObjs(self):
		kwargs = {"refPosGetter":self.refPosGetter, "pitch":self.pitch, "roll":self.roll,
		          "azimuthal":self.azimuthal}
		self.testObjA = tCode.WaterAdsorbateStandard(self.ohDists, self.angle, **kwargs)

	def testWithOnlyInternalCoords(self):
		expGeom = self._loadExpGeomA()
		actGeom = self.testObjA.geom
		self._checkGeomsEqual(expGeom,actGeom)

	def testWith90DegRoll(self):
		self.roll = 90
		self.createTestObjs()

		contribA = (self.ohDists[0]) / math.sqrt(2)
		contribB = (self.ohDists[1]) / math.sqrt(2)
		expGeom = [ [0.0, 0.0, 0.0, "O"],
		            [contribA, 0.0, contribA, "H"],
		            [contribB, 0.0, -1*contribB, "H"] ]
		actGeom = self.testObjA.geom
		self._checkGeomsEqual(expGeom,actGeom)

	def testWith90DegPitch(self):
		self.pitch = 90
		self.createTestObjs()

		contribA = (self.ohDists[0]) / math.sqrt(2)
		contribB = (self.ohDists[1]) / math.sqrt(2)
		expGeom = [ [0.0,0.0,0.0,"O"],
					[0.0, contribA, contribA, "H"],
					[0.0, -1*contribB, contribB, "H"] ]
		actGeom = self.testObjA.geom
		self._checkGeomsEqual(expGeom,actGeom)

	def testWithMinus90DegAzimuthal(self):
		self.azimuthal = -90
		self.createTestObjs()

		contribA = (self.ohDists[0]) / math.sqrt(2)
		contribB = (self.ohDists[1]) / math.sqrt(2)
		expGeom = [ [0.0, 0.0, 0.0, "O"],
		            [contribA,-1*contribA, 0.0, "H"],
		            [-1*contribB, -1*contribB, 0.0, "H"] ] 
		actGeom = self.testObjA.geom
		self._checkGeomsEqual(expGeom,actGeom)

	def testWith90DegAzimuthal(self):
		self.azimuthal = 90
		self.createTestObjs()

		contribA = (self.ohDists[0]) / math.sqrt(2)
		contribB = (self.ohDists[1]) / math.sqrt(2)
		expGeom = [ [0.0, 0.0, 0.0, "O"],
		            [-1*contribA, contribA, 0.0, "H"],
		            [contribB, contribB, 0.0, "H"] ] 
		actGeom = self.testObjA.geom
		self._checkGeomsEqual(expGeom,actGeom)

	def testWith270DegAzimuthal(self):
		self.azimuthal = 270
		self.azimuthal = -90
		self.createTestObjs()

		contribA = (self.ohDists[0]) / math.sqrt(2)
		contribB = (self.ohDists[1]) / math.sqrt(2)
		expGeom  = [ [0.0, 0.0, 0.0, "O"],
		             [contribA, -1*contribA, 0.0, "H"],
		             [-1*contribB, -1*contribB, 0.0, "H"] ]

		actGeom = self.testObjA.geom
		self._checkGeomsEqual(expGeom,actGeom)

	def _loadExpGeomA(self):
		contribA = (self.ohDists[0]) / math.sqrt(2)
		contribB = (self.ohDists[1]) / math.sqrt(2)

		outGeom = [ [ 0.0, 0.0, 0.0, "O"],
		            [ contribA, contribA, 0.0, "H"],
		            [ contribB, -1*contribB, 0.0, "H"] ]
		return outGeom

	def _checkGeomsEqual(self, expGeom, actGeom):
		for exp,act in it.zip_longest(expGeom, actGeom):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp[:3],act[:3])]
			self.assertEqual(exp[-1],act[-1])



class TestGetInternalCoordsFromGeom(unittest.TestCase):

	def setUp(self):
		self.eleListA = ["O","H","H"]
		self.coordO = [0.0, 0.0, 0.0]
		self.coordH_one = [-1.0, 0.0, 1.0]
		self.coordH_two = [ 1.0,-0.0, 1.1]
		self.allCoords = [self.coordO, self.coordH_one, self.coordH_two]
		self.createTestObjs()

	def createTestObjs(self):
		self.testCoordsA = [ coords + [symbol] for coords,symbol in it.zip_longest(self.allCoords, self.eleListA)]

	def testCheckWaterRaisesWhenLenWrong(self):
		self.eleListA.append(None)
		self.allCoords.append([0,0,0])
		self.createTestObjs()
		with self.assertRaises(ValueError):
			tCode.getWaterInternalCoordsFromXYZ(self.testCoordsA)

	def testCheckWaterRaisesWhenEleWrong(self):
		self.eleListA = ["O","H","O"]
		self.createTestObjs()
		with self.assertRaises(ValueError):
			tCode.getWaterInternalCoordsFromXYZ(self.testCoordsA)

	def testDDistsAndInternalAngleMatchesExpectedA(self):
		expDists = [ math.sqrt(2), math.sqrt(1 + (1.1*1.1)) ]
		expAngle = 87.27368900609373
		actDists, actAngle = tCode.getWaterInternalCoordsFromXYZ(self.testCoordsA)
		self.assertAlmostEqual(expAngle,actAngle)
		for exp,act in it.zip_longest(expDists,actDists):
			self.assertAlmostEqual(exp,act)
