
import itertools as it

import unittest
import unittest.mock as mock

import plato_pylib.utils.supercell as supCellHelp
import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.special_builders.fill_box as tCode

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


