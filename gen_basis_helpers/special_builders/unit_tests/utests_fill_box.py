
import itertools as it

import unittest
import unittest.mock as mock

import plato_pylib.utils.supercell as supCellHelp
import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.special_builders.fill_box as tCode



class TestGetCommonTripletOfFactors(unittest.TestCase):

	def setUp(self):
		self.testVal = 17

	def _runTestFunct(self):
		return tCode.getCommonTripletOfFactors(self.testVal)

	def testCaseForPrimeNumberA(self):
		expVals = [ [17,1,1], [1,17,1], [1,1,17] ]
		actVals = self._runTestFunct()
		self.assertEqual( sorted(expVals), sorted(actVals) )

	def testCaseForVal4(self):
		self.testVal = 4
		expVals = [ [4,1,1], [1,4,1], [1,1,4],
		            [2,2,1], [1,2,2], [2,1,2] ]
		actVals = self._runTestFunct()
		self.assertEqual( sorted(expVals), sorted(actVals) )

	def testForLargerCase(self):
		self.testVal = 12
		expVals = [ [12,1,1], [1,12,1], [1,1,12],
		            [6,2,1], [6,1,2], [2,6,1], [2,1,6], [1,6,2], [1,2,6],
		            [4,3,1], [4,1,3], [3,4,1], [3,1,4], [1,3,4], [1,4,3],
		            [3,2,2], [2,3,2], [2,2,3] ]
		actVals = self._runTestFunct()
		self.assertEqual( sorted(expVals), sorted(actVals) )

	def testFor24(self):
		self.testVal = 24
		expVals = [ [24,1,1], [1,24,1], [1,1,24],
		            [12,2,1], [12,1,2], [1,12,2], [1,2,12], [2,1,12], [2,12,1],
		            [8,3,1], [8,1,3], [3,8,1], [3,1,8], [1,3,8], [1,8,3],
		            [6,4,1], [6,1,4], [1,6,4], [1,4,6], [4,6,1], [4,1,6],
		            [4,3,2], [4,2,3], [3,2,4], [3,4,2], [2,3,4], [2,4,3], #breakdown of [12,2,1]
		            [6,2,2], [2,6,2], [2,2,6] ] #Breakdown of [12,2,1]

		actVals = self._runTestFunct()
		self.assertEqual( sorted(expVals), sorted(actVals) )

	def testForProductThreePrimes(self):
		self.testVal = 105
		self.assertEqual( self.testVal, 3*5*7 )
		expVals = [ [105,1,1], [1,105,1], [1,1,105],
		            [35,3,1], [35,1,3], [3,35,1], [3,1,35], [1,35,3], [1,3,35],
		            [21,5,1], [21,1,5], [5,1,21], [5,21,1], [1,5,21], [1,21,5],
		            [15,7,1], [15,1,7], [7,15,1], [7,1,15], [1,7,15], [1,15,7],
		            [7,5,3], [7,3,5], [5,7,3], [5,3,7], [3,5,7], [3,7,5] ]  #Breakdown of [35,3,1]/[21,5,1]/[15,7,1]

		actVals = self._runTestFunct()
		self.assertEqual( sorted(expVals), sorted(actVals) )

class TestMapPrimToInpCellMinimisingAverageLatticeParamDev(unittest.TestCase):

	def setUp(self):
		self.primLattParamsA = [1,1,1]
		self.primAnglesA = [90,90,90]
		self.inpCellLattParamsA = [2,2,12]
		self.inpCellLattAnglesA = [90,90,90]
		self.primFractCoordsA = [[0.5,0.5,0.5,"X"],
		                         [0.7,0.7,0.7,"X"]]
		self.targNAtoms = 120
		self.createTestObjs()

	def createTestObjs(self):
		currKwargs = {"lattParams":self.primLattParamsA, "lattAngles":self.primAnglesA}
		self.primCellA = uCellHelp.UnitCell(**currKwargs)
		self.primCellA.fractCoords = self.primFractCoordsA
		currKwargs = {"lattParams":self.inpCellLattParamsA, "lattAngles":self.inpCellLattAnglesA}
		self.inpCellA = uCellHelp.UnitCell(**currKwargs)
		self.inpCellA.fractCoords = self.primFractCoordsA
		self.testObjA = tCode.MapPrimToInpCellToMinimiseAverageLatticeParamDeviation()

	#Example comes from a real problem found with old algorithm, whereby a 
	#5x1x12 cell was constructued which led to ridic. small z-gaps
	def testGetNumbImagesC_perfectCellMatchGivesOddAbMatch(self):
		expDims = [2,2,15] # [5,1,12] would be given based purely on old algorithm; which considered match to c lattice vector in isolation
		actDims = self.testObjA(self.targNAtoms, self.primCellA, self.inpCellA)
		self.assertEqual(expDims, actDims)

	def testGetNumbImages_multipleDegenerateOptions(self):
		self.inpCellLattParamsA = [4,2,4]
		self.targNAtoms = (2)*(4*2*2)
		self.createTestObjs()
		expDims = [2,2,4] #
		actDims = self.testObjA(self.targNAtoms, self.primCellA, self.inpCellA)
		self.assertEqual(expDims, actDims)


class TestStandardGetCartCoordsToFillUpToNAtoms(unittest.TestCase):

	def setUp(self):
		self.primLattParamsA = [3,3,3]
		self.primAnglesA = [90,90,90]
		self.inpCellLattParamsA = [2*x for x in self.primLattParamsA]
		self.inpCellLattAnglesA = [90,90,90]

		self.primFractCoordsA = [[0.3,0.3,0.3,"X"],
		                          [0.6,0.6,0.6,"X"]]
		self.targNAtoms = 16
		self.createTestObjs()

	def createTestObjs(self):
		currKwargs = {"lattParams":self.primLattParamsA, "lattAngles":self.primAnglesA}
		self.primCellA = uCellHelp.UnitCell(**currKwargs)
		self.primCellA.fractCoords = self.primFractCoordsA
		currKwargs = {"lattParams":self.inpCellLattParamsA, "lattAngles":self.inpCellLattAnglesA}
		self.inpCellA = uCellHelp.UnitCell(**currKwargs)
		self.inpCellA.fractCoords = self.primFractCoordsA
		self.testObjA = tCode.CellFillerStandard(self.primCellA)

	def testGetNLayers_perfectMatchForPossibleVal(self):
		expLayers = 2
		actLayers = self.testObjA._getNLayersInC(self.targNAtoms, self.inpCellA)
		self.assertEqual(expLayers,actLayers)

	def testGetNLayers_noPerfectMatch_twoEquallyGoodMatches(self):
		self.targNAtoms = 18
		self.primFractCoordsA = [ [0,0,0,"X"] for x in range(6) ] #6 atom prim cell
		self.createTestObjs()

		expLayers = 3 #only 3/1 are allowed here and both are equally good. We pick the larger one though
		actLayers = self.testObjA._getNLayersInC(self.targNAtoms, self.inpCellA)
		self.assertEqual(expLayers,actLayers)

	def testGetNumbImagesAB_perfect2x2Match(self):
		nLayers = 2
		expVals = [2,2]
		actVals = self.testObjA._getNumbImagesAB(self.targNAtoms, nLayers, self.inpCellA)
		self.assertEqual(expVals,actVals)

	#na=2,nb=1 and na=1,nb=2 are equally good here
	def testGetNumbImagesAB_noPerfectMatch_twoEquallyGoodMatches(self):
		nLayers = 2
		dimProduct = 2 #na*nb
		self.targNAtoms = nLayers*dimProduct*len(self.primFractCoordsA)
		self.createTestObjs()
		expVal = [2,1]
		actVal = self.testObjA._getNumbImagesAB(self.targNAtoms, nLayers, self.inpCellA)
		self.assertEqual(expVal,actVal)

	#In this case the target cell is simply the 2x2x2 supercell. So sorta trivial to deal with
	def testExpectedCoordsForSimpleCell(self):
		outCell = supCellHelp.superCellFromUCell(self.primCellA,[2,2,2])
		expCartCoords = outCell.cartCoords
		actCartCoords = self.testObjA(self.targNAtoms, self.inpCellA)
		for exp,act in it.zip_longest(expCartCoords,actCartCoords):
			self.assertEqual( len(exp), len(act) )
			[self.assertAlmostEqual(e,a) for e,a in zip(exp[:3],act[:3])]
			self.assertEqual( exp[-1], act[-1] )

	def testErrorThrownForInvalidCellAngles(self):
		self.inpCellLattAnglesA = [90,90,120]
		self.createTestObjs()
		with self.assertRaises(ValueError):
			self.testObjA.getCartCoordsToFillUpCell(self.targNAtoms, self.inpCellA)

	def testErrorThrownForImpossibleNAtomsValue(self):
		self.targNAtoms = 5
		self.createTestObjs()
		with self.assertRaises(ValueError):
			self.testObjA.getCartCoordsToFillUpCell(self.targNAtoms, self.inpCellA)


