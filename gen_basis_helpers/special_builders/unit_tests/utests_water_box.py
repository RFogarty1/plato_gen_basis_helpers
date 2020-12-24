
import copy
import types
import unittest
import unittest.mock as mock

import plato_pylib.utils.supercell as supCellHelp
import plato_pylib.shared.ucell_class as uCellHelp


import gen_basis_helpers.adsorption.adsorbate_detectors as adsDetectHelp
import gen_basis_helpers.adsorption.surface_detectors as surfDetectHelp
import gen_basis_helpers.special_builders.water_box as tCode


class TestGetLatticeParameterForTargetDensity(unittest.TestCase):

	def setUp(self):
		self.lattParamsA = [1,1,2]
		self.lattAnglesA = [90,90,90]
		self.targDensity = 10
		self.nWater = 1
		self.varyLattParam = "c"
		self.massDictA = {"O":1,"H":0}
		self.createTestObjs()

	def createTestObjs(self):
		self.testCellA = uCellHelp.UnitCell(lattParams=self.lattParamsA, lattAngles=self.lattAnglesA)

	def runTestFunct(self):
		args = [self.testCellA, self.nWater, self.targDensity]
		kwargs = {"lattParam":self.varyLattParam, "massDict":self.massDictA}
		return tCode.findLatticeParameterToGetTargetDensityForNWater(*args, **kwargs)

	@mock.patch("gen_basis_helpers.special_builders.water_box.uConvHelp")
	def testExpValForSimpleCell(self, mockedUConv):
		mockedUConv.AVOGADRO_NUMBER = 1
		expVal = 0.1
		actVal = self.runTestFunct()
		self.assertAlmostEqual(expVal, actVal)


class TestGetBulkWaterEmptyBoxForSurfWithAdsBothSides(unittest.TestCase):

	def setUp(self):
		self.lattParamsA = [10,10,10]
		self.lattAnglesA = [90,90,90]
		self.adsSymbol = "ads"
		self.surfAtomSymbol = "surf"
		self.cartCoordsA = [ [0,0,8,self.adsSymbol],
		                     [0,0,7,self.adsSymbol],
		                     [0,0,6,self.surfAtomSymbol],
		                     [0,0,5,self.surfAtomSymbol],
		                     [0,0,4,self.adsSymbol],
		                     [0,0,3,self.adsSymbol] ]

		self.densityA = mock.Mock()
		self.massDict = None
		self.extraHeight = 0
		self.targDensity = 2
		self.nWaterTotal = 30
		self.createTestObjs()

	def createTestObjs(self):
		self.testCellA = uCellHelp.UnitCell(lattParams=self.lattParamsA, lattAngles=self.lattAnglesA)
		self.testCellA.cartCoords = self.cartCoordsA
		self._createAdsorbateDetectors()
		currArgs = [self.adsDetectorTopA, self.adsDetectorBotA, self.densityA]
		self.testObjA = tCode.GetBulkWaterEmptyBoxFromSurfWithAdsorbatesBothSidesStandard(*currArgs)

	def _createAdsorbateDetectors(self):
		surfDetector = surfDetectHelp.DetectSurfaceBasedOnElementsPresent([self.surfAtomSymbol])
		topFilterFunct = adsDetectHelp.FilterToAtomsAboveOrBelowSurf(surfDetector, top=True)
		botFilterFunct = adsDetectHelp.FilterToAtomsAboveOrBelowSurf(surfDetector, top=False)
		self.adsDetectorTopA = adsDetectHelp.DetectSimpleAtomicAdsorbateFromInpGeom(self.adsSymbol,postFilterFuncts=[topFilterFunct])
		self.adsDetectorBotA = adsDetectHelp.DetectSimpleAtomicAdsorbateFromInpGeom(self.adsSymbol,postFilterFuncts=[botFilterFunct])

	@mock.patch("gen_basis_helpers.special_builders.water_box.findLatticeParameterToGetTargetDensityForNWater")
	def testForSimpleCellA(self, mockedGetLatticeParam):
		mockHeightTotal = 20
		mockedGetLatticeParam.side_effect = lambda *args,**kwargs: mockHeightTotal
		expNWaterTotal = self.nWaterTotal - 4 #4 is the number of self.adsSymbol
		totalHeightAdsLayers = 1+1 #Theres a gap of 1 between ads above and below the cell
		expOutCell = uCellHelp.UnitCell(lattParams=[10,10,mockHeightTotal-totalHeightAdsLayers], lattAngles=[90,90,90])
		actOutCell,unused = self.testObjA(self.testCellA, self.nWaterTotal, targDensity=self.targDensity)
		mockedGetLatticeParam.assert_called_with(self.testCellA, expNWaterTotal, self.targDensity, lattParam="c", massDict=self.massDict)
		expOutCell.fractCoords = list()
		self.assertEqual(expOutCell, actOutCell)

	@mock.patch("gen_basis_helpers.special_builders.water_box.findLatticeParameterToGetTargetDensityForNWater")
	def testForSimpleCellA_plusFourHeight(self, mockedGetLatticeParam):
		mockHeightTotal = 20
		mockedGetLatticeParam.side_effect = lambda *args,**kwargs:mockHeightTotal
		totalHeightAdsLayers = 1+1
		extraHeight = 4
		expOutCell = uCellHelp.UnitCell(lattParams=[10,10,mockHeightTotal-totalHeightAdsLayers+extraHeight], lattAngles=[90,90,90])
		actOutCell,unused = self.testObjA(self.testCellA, self.nWaterTotal, targDensity=self.targDensity,extraHeight=extraHeight)
		expOutCell.fractCoords = list()
		self.assertEqual(expOutCell, actOutCell)		

	@mock.patch("gen_basis_helpers.special_builders.water_box.findLatticeParameterToGetTargetDensityForNWater")
	def testExpectedNumberOfWaterReturned(self, mockedGetLatticeParam):
		mockHeightTotal = 20
		mockedGetLatticeParam.side_effect = lambda *args,**kwargs:mockHeightTotal
		expNumberWater = self.nWaterTotal - 4
		unused, actNumberWater = self.testObjA(self.testCellA, self.nWaterTotal, targDensity=self.targDensity)
		self.assertEqual(expNumberWater,actNumberWater)



class TestGetWaterBoxForMD(unittest.TestCase):

	def setUp(self):
		self.adsObjA = types.SimpleNamespace( geom=[ [0,0,0,"ads"] ] )
		self.lattParamsA, self.lattAnglesA = [5,5,5], [90,90,90]
		self.cartCoordsA = [[3,3,3,"X"]]
		self.waterAdsGapA = 0
		self.nWater = 2
		self.createTestObjs()

	def createTestObjs(self):
		self.singleCellA = uCellHelp.UnitCell(lattParams=self.lattParamsA, lattAngles=self.lattAnglesA)
		self.singleCellA.cartCoords = self.cartCoordsA
		self.adsObjsA = [self.adsObjA]
		self.testObjA = tCode.GetWaterBoxForMDFromEmptyBoxStandard(self.singleCellA, self.adsObjsA,waterAdsGap=self.waterAdsGapA)
		self.emptyCellA = supCellHelp.superCellFromUCell(self.singleCellA,[2,1,1])
		self.emptyCellA.cartCoords = list()

	def testExpectedForSimpleOneLatticeSiteCell(self):
		expCartCoords = [ [3,3,3,"ads"], [8,3,3,"ads"] ]
		expCell = copy.deepcopy(self.emptyCellA)
		expCell.cartCoords = expCartCoords
		actCell = self.testObjA(self.emptyCellA, self.nWater)
		self.assertEqual(expCell,actCell)

	@mock.patch("gen_basis_helpers.special_builders.water_box.boxFillHelp.CellFillerImproved")
	def testExpectedCallsForSimpleWithAdsGap(self, mockedCellFiller):
		self.waterAdsGapA = 1
		self.createTestObjs()
		mockFillObj = mock.Mock()
		mockedCellFiller.side_effect = lambda *args, **kwargs: mockFillObj
		mockFillObj.side_effect = lambda *args,**kwargs: list()
		expCellToFill = uCellHelp.UnitCell( lattParams=[10,5,3],lattAngles=[90,90,90] )
		expCellToFill.cartCoords = list()
		self.testObjA(self.emptyCellA, self.nWater)
		mockFillObj.assert_called_with(self.nWater, expCellToFill)

	def testExpectedLattParamsWhenAdsGapDefined(self):
		self.waterAdsGapA = 0.1
		self.createTestObjs()
		expCell = uCellHelp.UnitCell(lattParams=[10,5,5-(2*self.waterAdsGapA)],lattAngles=[90,90,90])
		actCell = self.testObjA(self.emptyCellA,self.nWater)
		expCell.cartCoords = list()
		actCell.cartCoords = list()
		self.assertEqual(expCell,actCell)



class testGetBulkWaterFromSurfWithAdsorbed(unittest.TestCase):

	def setUp(self):
		self.mockInputCell = mock.Mock()
		self.expEmptyBox = mock.Mock()
		self.expBoxForMD = mock.Mock()
		self.expWaterAdsGap = 4
		self.nWaterTotal = 4
		self.nWaterForBox = self.nWaterTotal-2
		self.targDensity, self.extraHeight = 5, 7
		self.createTestObjs()

	def createTestObjs(self):
		self.mockGetEmptyBox = mock.Mock()
		self.mockFillEmptyBox = mock.Mock()
		self.mockGetEmptyBox.side_effect = lambda *args, **kwargs: (self.expEmptyBox, self.nWaterForBox)
		self.mockFillEmptyBox.side_effect = lambda *args, **kwargs: self.expBoxForMD
		self.mockFillEmptyBox.waterAdsGap = self.expWaterAdsGap
		currArgs = [self.mockGetEmptyBox, self.mockFillEmptyBox]
		self.testObjA = tCode.GetBulkWaterBoxesFromSurfWithAdsorbatesBothSidesStandard(*currArgs)


	def testExpectedOutputWithEverythingMocked(self):
		currKwargs = {"cellForMD":self.expBoxForMD, "emptyCell":self.expEmptyBox, "waterAdsGap":self.expWaterAdsGap}
		expObj = types.SimpleNamespace( **currKwargs )
		actObj = self.testObjA(self.mockInputCell, self.nWaterTotal,targDensity=self.targDensity,
		                       extraHeight=self.extraHeight, waterAdsGap=self.expWaterAdsGap)

		self.mockGetEmptyBox.assert_called_with(self.mockInputCell, self.nWaterTotal, targDensity=self.targDensity, extraHeight=self.extraHeight)
		self.mockFillEmptyBox.assert_called_with(self.expEmptyBox, self.nWaterForBox, waterAdsGap=self.expWaterAdsGap)
		self.assertEqual(expObj,actObj)




