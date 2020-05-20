#!/usr/bin/python3

""" Test specific implementations of surface objects """

import copy
import itertools as it
import math
import unittest
import unittest.mock as mock

import numpy as np
import plato_pylib.shared.ucell_class as UCell
import plato_pylib.utils.supercell as supCell

import gen_basis_helpers.shared.surfaces as tCode

class TestHcp001Surface(unittest.TestCase):

	def setUp(self):
		baseLattAngles = [90.0, 90.0, 120.0]
		baseLengths = [3.21,3.21,5.91]
		baseFractCoords = None
		self.baseUnitCell = UCell.UnitCell(lattParams=baseLengths, lattAngles=baseLattAngles, fractCoords = baseFractCoords)
		self.vacuumRegion = 4.0
		self.nLayers = 1
		self.createSurfaceObj()

	def createSurfaceObj(self):
		self.testObj = tCode.Hcp0001Surface(self.baseUnitCell, self.nLayers, self.vacuumRegion)

	def testWrongAnglesLeadsToValueError(self):
		lattAngles = self.baseUnitCell.lattAngles
		lattAngles["gamma"] = 110.0 #Not actually hcp if this is the angle, so should break
		self.baseUnitCell.lattAngles = lattAngles
		with self.assertRaises(ValueError):
			self.createSurfaceObj()

	def testExpectedSurfaceAreaGiven(self):
		lattParams = self.baseUnitCell.lattParams
		lattParams["b"] = 2*lattParams["b"] #As it would be in a supercell
		self.baseUnitCell.lattParams = lattParams
		self.createSurfaceObj()
		expSurfaceArea = 17.847224726270476
		actSurfaceArea = self.testObj.surfaceArea
		self.assertAlmostEqual(expSurfaceArea, actSurfaceArea)

	def testExpCartCoordsGiven(self):
		""" For hcp0001 the z-coords of each atom should be shifted by 0.5* vacuum length. This is easy to check
		for a single layer """
		self.baseUnitCell.fractCoords = [ [0.0,0.0,0.0,"Mg"], [0.33,0.66,0.5,"Mg"] ]
		startCartCoords = copy.deepcopy(self.baseUnitCell.cartCoords)
		self.createSurfaceObj()

		#Form the expected co-ordinates
		expCoords = np.array( copy.deepcopy([x[:3] for x in startCartCoords]) )
		expCoords[0,2] += 0.5*self.vacuumRegion
		expCoords[1,2] += 0.5*self.vacuumRegion
		actCoords = np.array( [x[:3] for x in self.testObj.unitCell.cartCoords] )

		self.assertTrue( np.allclose(expCoords,actCoords) )

	def testExpNumbAtomsInMultiLayerCell(self):
		""" Test we get the expected number of atoms when we apply multiple layers """
		self.baseUnitCell.fractCoords = [ [0.0,0.0,0.0,"Mg"], [0.33,0.66,0.5,"Mg"] ]
		self.nLayers = 4
		self.createSurfaceObj()

		#Test the initially constructed
		expectedNumbAtoms = 8
		actNumbAtoms = len( self.testObj.unitCell.cartCoords )
		self.assertEqual(expectedNumbAtoms, actNumbAtoms)

		#Test it after i've changed nLayers post construction
		self.testObj._nLayers = 3
		expectedNumbAtoms = 6
		actNumbAtoms = len( self.testObj.unitCell.cartCoords )
		self.assertEqual(expectedNumbAtoms, actNumbAtoms)



class TestBuildRocksalt001SingleLayerFromPrimitive(unittest.TestCase):

	def setUp(self):
		self.a = 4
		self.b = 4
		self.c = 4
		self.alpha = 60
		self.beta = 60
		self.gamma = 60
		self.eleA = "Mg"
		self.eleB = "O"
		self.createTestObjs()

	def createTestObjs(self):
		self.fractCoordsA = [ [0.0,0.0,0.0,self.eleA], [0.5,0.5,0.5,self.eleB] ]
		self.testCellA = UCell.UnitCell( lattParams=[self.a,self.b,self.c], lattAngles=[self.alpha,self.beta,self.gamma])
		self.testCellA.fractCoords = self.fractCoordsA #Unfortunately cant set this way upon initialization; in that case coords and eles need to be in separate iterables

	def _runFunct(self):
		return tCode.getSingleLayerRocksalt001FromPrimitiveCell(self.testCellA)

	@mock.patch("gen_basis_helpers.shared.surfaces.uCell.UnitCell")
	def testExpectedUCellCallForTestCellA(self, mockedUCellClass):
		expFractCoords = [ [0,0,0,self.eleA], [0,0,0.5,self.eleB], [0.5,0.5,0,self.eleB], [0.5,0.5,0.5,self.eleA] ]
		expLattParamAlongZ = (1/math.sqrt(0.5))*self.a #2x as many atoms along z, hence this increases
		expLattVects = [ [ (1/math.sqrt(2))*self.a   , (1/math.sqrt(2))*self.a, 0.0   ],
		                 [ -1*(1/math.sqrt(2))*self.a, (1/math.sqrt(2))*self.a, 0.0   ],
		                 [ 0                         , 0                      , expLattParamAlongZ] ]

		mockOutput = mock.Mock()
		mockedUCellClass.fromLattVects.side_effect = [mockOutput]

		actOutput = self._runFunct()
		inpArgs,kwargs = mockedUCellClass.fromLattVects.call_args
		actLattVects,actFractCoords = inpArgs

		#Check expected lattice vectors passed
		for exp,act in it.zip_longest(expLattVects, actLattVects):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp,act)]

		#Check expected fractional coordinates passed
		for exp,act in it.zip_longest(expFractCoords, actFractCoords):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp[:3],act[:3])]
			self.assertEqual(exp[-1], act[-1])

		self.assertEqual(mockOutput,actOutput)

	def testRaisesForWrongNumberAtoms(self):
		testCoords = self.testCellA.fractCoords
		testCoords.append( [0.0,0.0,0.0,"O"] )
		self.testCellA.fractCoords = testCoords
		self.assertTrue( len(self.testCellA.fractCoords)!=2 )
		with self.assertRaises(AssertionError):
			self._runFunct()


	def testRaisesForWrongAngles(self):
		self.alpha, self.beta, self.gamma = 90,90,90
		self.createTestObjs()
		with self.assertRaises(AssertionError):
			self._runFunct()

	def testRaisesForNonCubicLattParams(self):
		self.c = self.b*2 #Should all be equal for primitive rocksalt unit cell
		self.createTestObjs()
		with self.assertRaises(AssertionError):
			self._runFunct()



class TestRocksalt001Surface(unittest.TestCase):

	def setUp(self):
		self.eleA = "Mg"
		self.eleB = "O"
		self.a = 5.0 #Only need 1 lattice vector to start since its cubic
		self.lenVac = 10.0
		self.nLayers = 3
		self.createTestObjs()

	def createTestObjs(self):
		self.fractCoordsA = [ [0,0,0,self.eleA], [0,0,0.5,self.eleB], [0.5,0.5,0,self.eleB], [0.5,0.5,0.5,self.eleA] ]
		self.lattVectsA =  [ [ (1/math.sqrt(2))*self.a   , (1/math.sqrt(2))*self.a, 0.0   ],
		                     [ -1*(1/math.sqrt(2))*self.a, (1/math.sqrt(2))*self.a, 0.0   ],
		                     [ 0                         , 0                      , self.a] ]
		self.singleLayerCellA = UCell.UnitCell.fromLattVects(self.lattVectsA, self.fractCoordsA)

		self.testObjA = tCode.Rocksalt001Surface( self.singleLayerCellA, self.nLayers, self.lenVac )


	def testExpectedLattParams(self):
		expLattAParam, expLattBParam = self.a, self.a
		expLattCParam = (self.nLayers*self.a) + self.lenVac
		expLattParams = {"a":expLattAParam, "b":expLattBParam, "c":expLattCParam}

		outUCell = self.testObjA.unitCell
		actLattParams = outUCell.lattParams

		for key in expLattParams.keys():
			exp,act = expLattParams[key], actLattParams[key]
			self.assertAlmostEqual(exp,act)

	def testExpSurfaceArea(self):
		expSurfArea = self.singleLayerCellA.volume / self.a #Surface area is independent of nLayers and lenVacuum
		actSurfArea = self.testObjA.surfaceArea
		self.assertAlmostEqual(expSurfArea, actSurfArea)

class TestHcp1010FromPrimCell(unittest.TestCase):

	def setUp(self):
		self.inpFractCoords = [ [0.0,0.0,0.0],
		                        [1/3, 2/3, 0.5] ]
		self.inpAtoms = ["Mg" for x in self.inpFractCoords]

		self.lattParams = [5.0, 5.0, 7.5]
		self.lattAngles = [90.0, 90.0, 120.0]
		self.createTestObjs()

	def createTestObjs(self):
		fractCoords = [ x+[y] for x,y in it.zip_longest(self.inpFractCoords, self.inpAtoms) ] 
		self.testUCellA = UCell.UnitCell(lattAngles=self.lattAngles, lattParams = self.lattParams)
		self.testUCellA.fractCoords = fractCoords

	@mock.patch("gen_basis_helpers.shared.surfaces._centreCFractCoordsForInpCell")
	def testExpectedSingleLayer(self, mockedCentreFractCoords):
		a,b,c = self.lattParams
		expLattParams = [c,b,a]
		expLattAngles = [60,90,90]
		expFractCoords = copy.deepcopy(self.inpFractCoords)
		expFractCoords[1][0], expFractCoords[1][2] = self.inpFractCoords[1][2], -1*self.inpFractCoords[1][0]
		expUCell = UCell.UnitCell(lattParams=expLattParams, lattAngles=expLattAngles, fractCoords=expFractCoords, elementList=self.inpAtoms)
		actUCell = tCode.getSingleLayerHcp1010FromPrimitiveCell(self.testUCellA)
		self.assertEqual(expUCell,actUCell)
		mockedCentreFractCoords.assert_called_once_with(actUCell)


	def testRaisesIfABLattParamsNotEqual(self):
		self.lattParams[0] = 4
		self.lattParams[1] = 3
		self.createTestObjs()
		with self.assertRaises(AssertionError):
			tCode.getSingleLayerHcp1010FromPrimitiveCell(self.testUCellA)

	def testRaisesIfNumbAtomsNotEqualToTwo(self):
		self.inpFractCoords = [ [ 0.0,0.0,0.0 ] ]
		self.inpAtoms = ["Mg"]
		self.createTestObjs()
		with self.assertRaises(AssertionError):
			tCode.getSingleLayerHcp1010FromPrimitiveCell(self.testUCellA)

	def testRaisesIfAnglesAll90(self):
		self.lattAngles = [90,90,90]
		self.createTestObjs()
		with self.assertRaises(AssertionError):
			tCode.getSingleLayerHcp1010FromPrimitiveCell(self.testUCellA)


class TestGenericSurface(unittest.TestCase):

	def setUp(self):
		self.lattParams = [6,4,20]
		self.lattAngles = [60, 90, 90]
		self.nLayers = 1
		self.lenVac = 5
		self.createTestObjs()

	def createTestObjs(self):
		self.uCellA = UCell.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.surfA = tCode.GenericSurface(self.uCellA, self.nLayers, self.lenVac)

	def testExpectedAreaForSquareLikeSurface(self):
		expArea = self.lattParams[0]*self.lattParams[1] #Since gamma is 90 degrees, these simply form a square in this case
		actArea = self.surfA.surfaceArea
		self.assertAlmostEqual(expArea,actArea)

	def testExpectedAreaForParralelogramSurface(self):
		self.lattAngles = [90,90,60]
		self.createTestObjs()
		expArea = self.uCellA.volume/self.lattParams[-1] #This works only when c is 90 degrees to both a and b
		actArea = self.surfA.surfaceArea
		self.assertAlmostEqual(expArea,actArea)

class TestAddingVacuumRegion(unittest.TestCase):

	def setUp(self):
		self.lenVac = 20
		self.lattParams = [7.5,5.0,5.0]
		self.lattAngles = [60,90,90]
		self.createTestObjs()

	def createTestObjs(self):
		self.testCellA = UCell.UnitCell( lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.testCellA.fractCoords = [[0,0,0,"Mg"]]

#putCAlongZ=True

	def testExpectedDistanceBetweenAtomsForCubicCell(self):
		self.lattParams = [5.0,5.0,5.0]
		self.lattAngles = [90,90,90]
		self.createTestObjs()
		self._checkExpectedDistanceBetweenImages()

	def testExpectedDistanceBetweenAtomsForHexagonalCell(self):
		self._checkExpectedDistanceBetweenImages()

	def _checkExpectedDistanceBetweenImages(self):
		noVacSupercell = supCell.superCellFromUCell( self.testCellA, [1,1,2] )
		noVacDist = self._getDistTwoPoints(noVacSupercell.cartCoords[0][:3] ,noVacSupercell.cartCoords[1][:3])
		expDistance = self.lenVac + noVacDist
		tCode.addVacuumToUnitCellAlongC(self.testCellA,self.lenVac)
		doubleCell = supCell.superCellFromUCell( self.testCellA, [1,1,2] )
		posA, posB = doubleCell.cartCoords[0][:3], doubleCell.cartCoords[1][:3]
		actDistance = self._getDistTwoPoints(posA,posB)
		self.assertAlmostEqual(expDistance,actDistance) 

	def _getDistTwoPoints(self, pointA, pointB):
		return math.sqrt( sum( [x**2 for x in [b-a for b,a in it.zip_longest(pointA,pointB)]] ) )


class TestCentreFractCoordsAlongC(unittest.TestCase):

	def setUp(self):
		self.lattParams = [2,2,2]
		self.lattAngles = [90,90,90]
		self.fractCoords = [ [0.5,0.5,0.3],
		                     [0.6,0.6,0.5] ]
		self.atomList = ["Mg" for x in self.fractCoords]
		self.createTestObjs()

	def createTestObjs(self):
		self.testCellA = UCell.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles,
		                                fractCoords=self.fractCoords, elementList=self.atomList)

	def testForSimpleSetA(self):
		expCoords = [0.4, 0.6]
		self._testExpMatchesAct(expCoords)

	def testForFractCoordsAboveOne(self):
		self.fractCoords[0][-1] = 1.3
		self.createTestObjs()
		expCoords = [0.9,0.1]
		self._testExpMatchesAct(expCoords)

	def testForFractCoordsBelowZero(self):
		self.fractCoords[0][-1] = -0.4
		self.createTestObjs()
		expCoords = [0.05,0.95]
		self._testExpMatchesAct(expCoords)

	def testForBothAboveOne(self):
		self.fractCoords[0][-1] = 1.4
		self.fractCoords[1][-1] = 1.5
		self.createTestObjs()
		expCoords = [0.45, 0.55]
		self._testExpMatchesAct(expCoords)

	def testForBothBelowZero(self):
		self.fractCoords[0][-1] = -0.5
		self.fractCoords[1][-1] = -0.6
		self.createTestObjs()
		expCoords = [0.55,0.45]
		self._testExpMatchesAct(expCoords)

	def _testExpMatchesAct(self, expCoords):
		tCode._centreCFractCoordsForInpCell(self.testCellA)
		actCoords = self._getZFractValsFromUnitCell(self.testCellA)
		for exp,act in it.zip_longest(expCoords,actCoords):
			self.assertAlmostEqual(exp,act)

	def _getZFractValsFromUnitCell(self, inpCell):
		fractCoords = inpCell.fractCoords
		zCoords = [x[2] for x in fractCoords]
		return zCoords



if __name__ == '__main__':
	unittest.main()
