
import itertools as it
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp
import gen_basis_helpers.shared.surfaces as surfHelp

import gen_basis_helpers.adsorption.hcp0001 as tCode

class TestSurfaceToSite_hcp(unittest.TestCase):

	def setUp(self):
		self.lattParamsA = [2,2,3]
		self.lattAnglesA = [90,90,120]
		self.fractCoordsA = [ [0.0,0.0,0.0],
		                     [1/3, 2/3, 0.5] ]
		self.eleListA = ["Mg", "Mg"]
		self.nLayersA = 2
		self.absVacLengthA = 10
		self.top = True
		self.createTestObjs()

	def createTestObjs(self):
		self._createUnitCells()
		self.surfA = surfHelp.GenericSurface(self.cellA, self.nLayersA, lenAbsoluteVacuum=self.absVacLengthA)
		self.surfToSiteObj = tCode.HcpSurfaceToHcpSites(top=self.top)

	def _createUnitCells(self):
		cellKwargDict = {"lattParams":self.lattParamsA, "lattAngles":self.lattAnglesA,
		                 "fractCoords":self.fractCoordsA, "elementList":self.eleListA}
		self.cellA = uCellHelp.UnitCell(**cellKwargDict)
	
	def testGetSurfVectorTop(self):
		expVector = [0,0,1]
		actVector = self.surfToSiteObj.getOutwardsSurfaceVectorFromSurface(self.surfA)
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expVector,actVector)]
	
	#use max(z) as the top surface layer; trick that only really works for this type of cell but..
	def testForSimpleCell(self):
		expSitePosition = self._getExpectedPosForSimpleCellTopA()
		actPositions = self.surfToSiteObj(self.surfA)
		self.assertTrue( len(actPositions)==1 )
		for exp,act in it.zip_longest(expSitePosition, actPositions[0]):
			self.assertAlmostEqual(exp,act)

	def testForSimpleCellBottomSurface(self):
		self.top=False
		self.createTestObjs()
		expSitePosition  = self._getExpectedPosForSimpleCellBottomA()
		actPositions = self.surfToSiteObj(self.surfA)
		for exp,act in it.zip_longest(expSitePosition, actPositions[0]):
			self.assertAlmostEqual(exp,act)

	def testFor2x2Cell(self):
		self.fractCoordsA = [ [0.0, 0.0, 0.0],
		                      [1/6, 1/3, 0.5],
		                      [0.5, 0.0, 0.0],
		                      [2/3, 1/3, 0.5],
		                      [0.0, 0.5, 0.0],
		                      [1/6, 5/6, 0.5],
		                      [0.5, 0.5, 0.0],
		                      [2/3, 5/6, 0.5] ]
		self.eleListA = ["Mg" for x in self.fractCoordsA]
		self.nLayers = 1
		self.createTestObjs()
		secondLayerIndices = [0,2,4,6] 
		expXYVals = [self.surfA.unitCell.cartCoords[idx][:2] for idx in secondLayerIndices]
		expZVal = max([x[2] for x in self.surfA.unitCell.cartCoords])
		expPositions = [list(x)+[expZVal] for x in expXYVals]
		actPositions = self.surfToSiteObj(self.surfA)

		for expPos, actPos in it.zip_longest(expPositions,actPositions):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expPos,actPos)]

	def _getExpectedPosForSimpleCellTopA(self):
		secondLayerAtomIdx = 2
		expZVal = max([x[2] for x in self.surfA.unitCell.cartCoords]) #only works because alpha/beta are 90 degrees
		expXYVals = self.surfA.unitCell.cartCoords[secondLayerAtomIdx][:2]
		expSitePosition = expXYVals + [expZVal]
		return expSitePosition

	def _getExpectedPosForSimpleCellBottomA(self):
		secondLayerAtomIdx = 1
		expZVal = min([x[2] for x in self.surfA.unitCell.cartCoords])
		expXYVals = self.surfA.unitCell.cartCoords[secondLayerAtomIdx][:2]
		expSitePosition = expXYVals + [expZVal]
		return expSitePosition

class TestSurfaceToSite_fccHollow(unittest.TestCase):


	def setUp(self):
		self.lattParamsA =  [12.00001,12,9]
		self.lattAnglesA = [90,90,120]
		self.fractCoordsA = [ [0.0, 0.0, 0.0],
		                      [1/6, 1/3, 0.5],
		                      [0.5, 0.0, 0.0],
		                      [2/3, 1/3, 0.5],
		                      [0.0, 0.5, 0.0],
		                      [1/6, 5/6, 0.5],
		                      [0.5, 0.5, 0.0],
		                      [2/3, 5/6, 0.5] ]
		self.eleListA = ["Mg" for x in self.fractCoordsA]
		self.nLayersA = 1
		self.absVacLengthA = 10
		self.top = True
		self.createTestObjs()

	def createTestObjs(self):
		self._createUnitCells()
		self.surfA = surfHelp.GenericSurface(self.cellA, self.nLayersA, lenAbsoluteVacuum=self.absVacLengthA)
		self.surfToSiteObj = tCode.HcpSurfaceToFccHollowSites(top=self.top)

	def _createUnitCells(self):
		cellKwargDict = {"lattParams":self.lattParamsA, "lattAngles":self.lattAnglesA,
		                 "fractCoords":self.fractCoordsA, "elementList":self.eleListA}
		self.cellA = uCellHelp.UnitCell(**cellKwargDict)


	#Note i got the expected values by doing it near-"manually" and visualising results
	def testFor2x2Cell(self):

		expSites = [ [6.0, 6.9282032302755105, 7.25],
		             [0.0, 6.9282032302755105, 7.25],
		             [9.0, 1.7320508075688785, 7.25],
		             [3.0, 1.7320508075688785, 7.25] ]

		actSites = self.surfToSiteObj(self.surfA)
		self.assertEqual( len(expSites), len(actSites) )

		for exp,act in it.zip_longest(expSites, actSites):
			[self.assertAlmostEqual(e,a,places=4) for e,a in it.zip_longest(exp,act)]


class TestSurfaceToSites_atop(unittest.TestCase):

	def setUp(self):
		self.lattParamsA = [2,2,3]
		self.lattAnglesA = [90,90,120]
		self.fractCoordsA = [ [0.0,0.0,0.0],
		                     [1/3, 2/3, 0.5] ]
		self.eleListA = ["Mg", "Mg"]
		self.nLayersA = 2
		self.absVacLengthA = 10
		self.top = True
		self.createTestObjs()

	def createTestObjs(self):
		self._createUnitCells()
		self.surfA = surfHelp.GenericSurface(self.cellA, self.nLayersA, lenAbsoluteVacuum=self.absVacLengthA)
		self.surfToSiteObj = tCode.HcpSurfaceToAtopSites(top=self.top)

	def _createUnitCells(self):
		cellKwargDict = {"lattParams":self.lattParamsA, "lattAngles":self.lattAnglesA,
		                 "fractCoords":self.fractCoordsA, "elementList":self.eleListA}
		self.cellA = uCellHelp.UnitCell(**cellKwargDict)

	def testForSimpleCell(self):
		expSitePosition = self._getExpectedPosForSimpleCellTopA()
		actPositions = self.surfToSiteObj(self.surfA)
		self.assertTrue( len(actPositions)==1 )
		for exp,act in it.zip_longest(expSitePosition, actPositions[0]):
			self.assertAlmostEqual(exp,act)

	def testFor2x2Cell(self):
		self.fractCoordsA = [ [0.0, 0.0, 0.0],
		                      [1/6, 1/3, 0.5],
		                      [0.5, 0.0, 0.0],
		                      [2/3, 1/3, 0.5],
		                      [0.0, 0.5, 0.0],
		                      [1/6, 5/6, 0.5],
		                      [0.5, 0.5, 0.0],
		                      [2/3, 5/6, 0.5] ]
		self.eleListA = ["Mg" for x in self.fractCoordsA]
		self.nLayersA = 1
		self.createTestObjs()
		firstLayerIndices = [1,3,5,7]
		expPositions = [list(self.surfA.unitCell.cartCoords[idx][:3]) for idx in firstLayerIndices]
		actPositions = self.surfToSiteObj(self.surfA)

		for expPos, actPos in it.zip_longest(expPositions,actPositions):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expPos,actPos)]


	def _getExpectedPosForSimpleCellTopA(self):
		firstLayerAtomIdx = -1 
		expSitePosition = list( self.surfA.unitCell.cartCoords[firstLayerAtomIdx][:3] )
		return expSitePosition

