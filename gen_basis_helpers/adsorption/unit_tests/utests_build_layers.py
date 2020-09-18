
import types
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp
import gen_basis_helpers.shared.surfaces as surfHelp

import gen_basis_helpers.adsorption.hcp0001 as hcpSurfHelp
import gen_basis_helpers.adsorption.add_adsorbates as addAdsorbateHelp
import gen_basis_helpers.adsorption.adsorbate_rep_objs as adsorbRepHelp

import gen_basis_helpers.adsorption.build_adsorbate_layers as tCode

class TestBuilderWithSimpleFullCoverageCase(unittest.TestCase):

	def setUp(self):
		self.createCellParamsA()
		self.createCell()
		self.surfA = surfHelp.GenericSurface(self.cellA, self.nLayersA, lenAbsoluteVacuum=self.absVacLengthA)

		self.surfToSites = hcpSurfHelp.HcpTopSurfaceToSites()
		self.adsorbateGeom = [[0.0,0.0,0.0,"O"]]
		self.fractCoverage = 1
		self.distance = 2
#		self.getAdsorbatesFunct = addAdsorbateHelp.SingleTypeGetAdsorbatesForSites

		self.createTestObjs()

	def createCellParamsA(self):
		self.lattParamsA = [2,2,3]
		self.lattAnglesA = [90,90,120]
		self.fractCoordsA = [ [0.0,0.0,0.0],
		                     [1/3, 2/3, 0.5] ]
		self.eleListA = ["Mg", "Mg"]
		self.nLayersA = 2
		self.absVacLengthA = 10
		self.top = True

	def createCell(self):
		cellKwargDict = {"lattParams":self.lattParamsA, "lattAngles":self.lattAnglesA,
		                 "fractCoords":self.fractCoordsA, "elementList":self.eleListA}
		self.cellA = uCellHelp.UnitCell(**cellKwargDict)

	def createTestObjs(self):
		self.adsorbateObj = createAdsorbateObjFromGeom(self.adsorbateGeom)
		self.getAdsorbatesFunct = addAdsorbateHelp.SingleTypeGetAdsorbatesForSites(self.adsorbateObj,self.fractCoverage, distance=self.distance)
		self.testObjA = tCode.AdsorbateLayerBuilder( self.surfToSites,self.getAdsorbatesFunct )

	def _getExpectedAdsorptionLayer(self):
		expPositions = [[1.8369701987210317e-16, 3.1817257161747195e-16, 8.75]] #Determined by the surface to sites
		adsorbates = [self.adsorbateObj]
		distances = [self.distance]
		surfVector = [0,0,1]
		return adsorbRepHelp.AdsorbateLayer(expPositions, adsorbates, distances, surfVector)

	def testExpectedLayerA(self):
		expObj = self._getExpectedAdsorptionLayer()
		actObj = self.testObjA.build(self.surfA)
		self.assertEqual(expObj,actObj)



def createAdsorbateObjFromGeom(inpGeom):
	return types.SimpleNamespace(geom=inpGeom)
