
import copy
import itertools as it
import math
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

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
		self.translationVector = [0,0,0]
		self.refPosGetter = tCode.WaterRefPosGetterStandard()
		self.createTestObjs()

	def createTestObjs(self):
		kwargs = {"refPosGetter":self.refPosGetter, "pitch":self.pitch, "roll":self.roll,
		          "azimuthal":self.azimuthal, "translationVector":self.translationVector}
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

	def testTranslationVectorA(self):
		self.translationVector = [1,2,3]
		self.createTestObjs()
		contribA = (self.ohDists[0]) / math.sqrt(2)
		contribB = (self.ohDists[1]) / math.sqrt(2)

		expGeom = [ [ 1.0, 2.0, 3.0, "O"],
		            [ contribA+1, contribA+2     , 3.0, "H"],
		            [ contribB+1, (-1*contribB)+2, 3.0, "H"] ]
		
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


class TestGetAxisRotationsFromInpXyz(unittest.TestCase):

	def setUp(self):
		self.ohDists = [1.1,1.1]
		self.angle = 105
		self.roll, self.pitch, self.azimuthal = 0,0,0
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"roll":self.roll, "pitch":self.pitch, "azimuthal":self.azimuthal}
		self.waterObjA = tCode.WaterAdsorbateStandard(self.ohDists,self.angle, **kwargDict)

	def testAllZeroWhenInRefGeom(self):
		expVals = [0.0, 0.0, 0.0]
		actVals = tCode.getStandardAxisRotationsFromXyz(self.waterObjA.geom)
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expVals,actVals)]

	def testOutputConsistentWith45Roll(self):
		self.waterObjA.roll = 45
		expVals = [45, 0, 0]
		actVals = tCode.getStandardAxisRotationsFromXyz(self.waterObjA.geom)
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expVals,actVals)]

	def testOutputConsistentWith45Pitch(self):
		self.waterObjA.pitch = 45
		expVals = [0, 45, 0]
		actVals = tCode.getStandardAxisRotationsFromXyz(self.waterObjA.geom)
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expVals,actVals)]

	def testOutputConsistentWith45Azimuthal(self):
		self.waterObjA.azimuthal = 45
		expVals = [0, 0, 45]
		actVals = tCode.getStandardAxisRotationsFromXyz(self.waterObjA.geom)
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expVals,actVals)]

	def testOutputConsitentWithAllNeg45(self):
		self.waterObjA.azimuthal, self.waterObjA.pitch, self.waterObjA.roll = -45, -45, -45
		expVals = [-45, -45, -45]
		actVals = tCode.getStandardAxisRotationsFromXyz(self.waterObjA.geom)
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expVals,actVals)]

	def testOutputConsistentWithRangeOfAngleCombos(self):
		inpCombos =[ [ 30,  40, 50],
		             [-20, -10, 12],
		             [14, 26, 19] ]
		for combo in inpCombos:
			self._runStandardTestForInpCombo(combo)

	def _runStandardTestForInpCombo(self, inpCombo):
		self.waterObjA.roll = inpCombo[0]
		self.waterObjA.pitch = inpCombo[1]
		self.waterObjA.azimuthal = inpCombo[2]
		actVals = tCode.getStandardAxisRotationsFromXyz(self.waterObjA.geom)
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(inpCombo,actVals)]

	#Badly behaved for angles close to 90; best is to just throw an error for now
	def testRaisesIfAngleNear90(self):
		self.waterObjA.roll = -89.999
		with self.assertRaises(ValueError):
			tCode.getStandardAxisRotationsFromXyz(self.waterObjA.geom)



class CheckExpAndActEqualMixin():

	def _checkExpAndActCoordsEqual(self, expCoords, actCoords):
		self.assertEqual( len(expCoords), len(actCoords) )
		for exp,act in zip(expCoords,actCoords):
			self.assertEqual( len(act), 4 )
			self.assertEqual(exp[-1],act[-1])
			[self.assertAlmostEqual(e,a) for e,a in zip(exp[:3],act[:3])]



class TestGetStandardWaterAdsObjFromXYZCoords(unittest.TestCase, CheckExpAndActEqualMixin):

	def setUp(self):
		self.ohDists = [1,2] #Different ohDists must slightly alter the azimuth value
		self.angle = 90
		self.pitch = 12
		self.roll = 25
		self.azimuthal = 18
		self.translationVector = [1,1,1]
		self.refPosGetter = tCode.WaterRefPosGetterStandard()
		self.createTestObjs()

	def createTestObjs(self):
		kwargs = {"refPosGetter":self.refPosGetter, "pitch":self.pitch, "roll":self.roll,
		          "azimuthal":self.azimuthal, "translationVector":self.translationVector}
		self.testObjA = tCode.WaterAdsorbateStandard(self.ohDists, self.angle, **kwargs)
		self.testXyzA = self.testObjA.geom

	def testExpWaterObjA(self):
		expXyz = self.testXyzA
		actObj = tCode.getStandardWaterAdsObjFromXyz(expXyz)
		self._checkExpAndActCoordsEqual(expXyz, actObj.geom)


#This means messing with WaterRefPosGetter return co-ords for the OXYGEN
class TestGetWaterAdsorptionObjsFromCellAndIndices(unittest.TestCase, CheckExpAndActEqualMixin):

	def setUp(self):
		self.lattParamsA = [10,10,10]
		self.lattAnglesA = [90,90,90]

		self.indicesA = [ [1,2,3] ]
		self.coordsA = [ [3,3,3,"X"], [4,4,9,"O"], [3.5,3.5,9,"H"], [4,3.5,9,"H"] ] 
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParamsA, lattAngles=self.lattAnglesA)
		self.cellA.cartCoords = self.coordsA

	def _runTestFunct(self):
		return tCode.getWaterAdsorptionObjsFromInpCellAndWaterIndices(self.cellA, self.indicesA)

	def testExpA_pbcsIrrelevant(self):
		expObjs = [ tCode.getStandardWaterAdsObjFromXyz(self.coordsA[1:]) ]
		expGeom = expObjs[0].geom
		actObjs = self._runTestFunct()
		actGeom = actObjs[0].geom
		self._checkExpAndActCoordsEqual(expGeom, actGeom)

	def testExpB_pbcsMatter(self):
		expObjs = [ tCode.getStandardWaterAdsObjFromXyz(self.coordsA[1:]) ]
		self.coordsA[2][-2] += self.lattParamsA[-1]
		self.createTestObjs()

		expGeom = expObjs[0].geom
		actObjs = self._runTestFunct()
		actGeom = actObjs[0].geom
		self._checkExpAndActCoordsEqual(expGeom, actGeom)


class TestWaterAdsorbateEquality(unittest.TestCase):

	def setUp(self):
		self.ohDists = [1,2] #Different ohDists must slightly alter the azimuth value
		self.angle = 90
		self.pitch = 12
		self.roll = 25
		self.azimuthal = 18
		self.translationVector = [1,1,1]
		self.refPosGetter = tCode.WaterRefPosGetterStandard()
		self.createTestObjs()

	def createTestObjs(self):
		kwargs = {"refPosGetter":self.refPosGetter, "pitch":self.pitch, "roll":self.roll,
		          "azimuthal":self.azimuthal, "translationVector":self.translationVector}
		self.testObjA = tCode.WaterAdsorbateStandard(self.ohDists, self.angle, **kwargs)

	def testEqualObjsCompareEqual(self):
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = self.testObjA
		self.assertEqual(objA, objB)

	def testUnequalObjsCompareUnequal_diffAngle(self):
		objA = copy.deepcopy(self.testObjA)
		self.angle += 12
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testUnequalObjsCompareUnequal_difftranslationVector(self):
		objA = copy.deepcopy(self.testObjA)
		self.translationVector[-1] += 2
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)


class TestGetAdsorbatesForNextWaterBilayer(unittest.TestCase):

	def setUp(self):
		#The "flat" one
		self.ohDistsFlat, self.angleFlat = [1.2, 1.2], 110
		self.pitchFlat, self.rollFlat, self.azimuthFlat = 10, 12, 20
		self.tVectFlat =  [2,3,4]

		#The "Hup" one
		self.ohDistsHup, self.angleHup = [1.1, 1.1], 105
		self.pitchHup, self.rollHup, self.azimuthHup = 30, 80, 30
		self.tVectHup = [3,4,5]

		#Args for the constructor object in essence
		self.bilayerSpacing = 2
		self.top = True
		self.surfNormal = [0,0,1]
		self.layerTol = 0.5

		self.createTestObjs()

	def createTestObjs(self):
		kwargsFlat = {"pitch": self.pitchFlat, "roll": self.rollFlat,
		              "azimuthal": self.azimuthFlat, "translationVector": self.tVectFlat}
		kwargsHup = {"pitch": self.pitchHup, "roll": self.rollHup,
		             "azimuthal": self.azimuthHup, "translationVector": self.tVectHup}
		self.flatAdsObj = tCode.WaterAdsorbateStandard(self.ohDistsFlat, self.angleFlat, **kwargsFlat)
		self.hupAdsObj  = tCode.WaterAdsorbateStandard(self.ohDistsHup , self.angleHup , **kwargsHup)
		self.adsObjs = [self.flatAdsObj, self.hupAdsObj]

	def _runTestFunct(self):
		args = [self.adsObjs, self.bilayerSpacing]
		kwargDict = {"layerTolerance":self.layerTol , "top": self.top, "surfNormal": self.surfNormal}
		return tCode.getAdsorbateObjsForNextWaterBilayerBasic(*args, **kwargDict)

	def _loadExpObjsA(self):
		layerHeight = abs( self.tVectHup[-1] - self.tVectFlat[-1] )

		#First is on top of the "flat" water molecule
		expAdsObjA = copy.deepcopy(self.hupAdsObj)
		expAdsObjA.translationVector = copy.deepcopy(self.tVectFlat)
		expAdsObjA.translationVector[-1] += self.bilayerSpacing + (2*layerHeight)
		expAdsObjA.azimuthal += 180

		#Second is on top of the "Hup" water molecule
		expAdsObjB = copy.deepcopy(self.flatAdsObj)
		expAdsObjB.translationVector = copy.deepcopy(self.tVectHup)
		expAdsObjB.translationVector[-1] += self.bilayerSpacing
		expAdsObjB.azimuthal += 180

		return [expAdsObjB, expAdsObjA]

	def testExpectedValuesA(self):
		expObjs = self._loadExpObjsA()
		actObjs = self._runTestFunct()
		self.assertEqual(expObjs,actObjs)

	def testRaisesIfTopAndBottomOfBilayerIndistinghuisable(self):
		self.layerTol = 2
		with self.assertRaises(ValueError):
			self._runTestFunct()

	def testRaisesIfMoreThanTwoLayers(self):
		extraAds = copy.deepcopy(self.flatAdsObj)
		extraAds.translationVector[-1] += 5

		self.adsObjs += [extraAds]
		with self.assertRaises(ValueError):
			self._runTestFunct()



