
import types
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp
import plato_pylib.shared.unit_convs as uConvHelp

import gen_basis_helpers.adsorption.detect_water_layer as tCode



class TestDetectTopWaterLayer(unittest.TestCase):

	def setUp(self):
		self.waterIdxDetector = mock.Mock()
		self.indicesA = [ [1,2,3] ]
		self.mockCell = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		self.waterIdxDetector.getIndicesFromInpGeom.side_effect = lambda *args,**kwargs: self.indicesA

		self.adsObjDetector = tCode.DetectOuterAdsorbedWaterLayer( self.waterIdxDetector )

	@mock.patch("gen_basis_helpers.adsorption.detect_water_layer.waterAdsHelp.getWaterAdsorptionObjsFromInpCellAndWaterIndices")
	def testExpectedCallsMadeCaseA(self, mockedGetAdsObjs):
		expOutObjs = mock.Mock()
		mockedGetAdsObjs.side_effect = lambda *args,**kwargs: expOutObjs

		actOutObj = self.adsObjDetector(self.mockCell)
		self.waterIdxDetector.getIndicesFromInpGeom.assert_called_with(self.mockCell)
		mockedGetAdsObjs.assert_called_with(self.mockCell, self.indicesA)
		self.assertEqual(expOutObjs, actOutObj)



#class TestDetectTopWaterBilayerStandard(unittest.TestCase):
#
#	def setUp(self):
#		self.waterIndicesDetector = None
#		self.top = True
#		self.maxHeightLayer = uConvHelp.ANG_TO_BOHR
#		self.createTestObjs()
#
#	def createTestObjs(self):
#		self.cellA = self._loadEmptyCellA()
#		coordsLower = self._loadBottomBilayerCartCoordsA()
#		coordsUpper = self._loadTopBilayerCartCoordsA()
#		coordsTotal = coordsLower + coordsUpper
#		self.cellA.cartCoords = coordsTotal
#
#	def _loadBottomBilayerCartCoordsA(self):
#		outCoords = [ [2.6119438170999993, -0.0343183798, 3.3663400374, 'O'],
#		              [2.3525372318, -0.0682797642, 4.3329796712, 'H'],
#		              [3.5926436985, -0.029955222900000005, 3.4091552648, 'H'],
#		              [6.483337342699999, -0.0591543147, 3.8786787253, 'O'],
#		              [6.8815437646, 0.7310404550000001, 3.4426481421, 'H'],
#		              [6.873130732000001, -0.8678739177999999, 3.4648117865000003, 'H'],
#		              [-2.2025775523000006, 2.7455932859999996, 3.3662972877, 'O'],
#		              [-2.4620360876, 2.711736769, 4.3329148404, 'H'],
#		              [-1.2218789034, 2.7499814100999997, 3.4091807255, 'H'],
#		              [-3.146861774599999, 5.5007036701, 3.8783281097, 'O'],
#		              [-2.7489202096000005, 6.2910100616, 3.4422333971000003, 'H'],
#		              [-2.7570449152000003, 4.6920706056, 3.4643292603000004, 'H'],
#		              [1.6678022701, 2.7207101356, 3.8780437935000003, 'O'],
#		              [2.0658094066999992, 3.5109435116000003, 3.4418965538000004, 'H'],
#		              [2.0576421117999995, 1.9120415708999998, 3.4641209280000003, 'H'],
#		              [2.6117366417999994, 5.5256505674, 3.3664091411000006, 'O'],
#		              [2.3522223678000005, 5.4913579954, 4.3330049633, 'H'],
#		              [3.5924353393000006, 5.5298175345, 3.4093499376, 'H'] ]
#		return outCoords
#
#	def _loadTopBilayerCartCoordsA(self):
#		outCoords = [ [2.6065907768, 1.8862477485, 22.9972621702, 'O'],
#		              [3.5870881677000006, 1.8876339681, 22.9492694116, 'H'],
#		              [2.3416029165, 1.9248998169, 22.0323093395, 'H'],
#		              [6.4796642789, 1.9155287628, 22.4882271666, 'O'],
#		              [6.870516332, 2.7250233098, 22.900046375, 'H'],
#		              [6.8751966951000005, 1.126084546, 22.9276770618, 'H'],
#		              [-2.2083679347, 4.666180903500001, 22.9972175667, 'O'],
#		              [-1.2278639243000002, 4.6676898889, 22.949346858, 'H'],
#		              [-2.4732552171, 4.705012570800001, 22.0322438545, 'H'],
#		              [-3.1507231300000003, 7.4754875135, 22.4886327349, 'O'],
#		              [-2.7598455673, 8.2849647083, 22.9004685052, 'H'],
#		              [-2.7552341034, 6.686050545900001, 22.9281265429, 'H'],
#		              [1.6644064100999998, 4.6954567087, 22.4885682286, 'O'],
#		              [2.0552790458, 5.5049303438, 22.9003967747, 'H'],
#		              [2.0600303593, 3.9060353599, 22.9279729401, 'H'],
#		              [2.6065467197, 7.446153251100001, 22.9971970741, 'O'],
#		              [3.5870533016, 7.4475744128, 22.9492496434, 'H'],
#		              [2.3416342786, 7.4847065013, 22.0322210852, 'H'] ]
#		return outCoords
#
#	def _loadEmptyCellA(self):
#		lattParams =[9.63, 9.63, 28.98]
#		lattAngles = [90, 90, 120]
#		outCell = uCellHelp.UnitCell(lattParams=lattParams, lattAngles=lattAngles)
#		return outCell
#
#	def _loadUpperWaterAdsorbates(self):
#		allCoords = self._loadTopBilayerCartCoordsA()
#		outAdsObjs = list()
#		idx = 0
#		while idx<len(allCoords):
#			currCoords = allCoords[idx:idx+3]
#			currAdsObj = types.SimpleNamespace(geom=currCoords)
#			outAdsObjs.append(currAdsObj)
#			idx += 3
#		return outAdsObjs
#
#	#TODO: First step will be to force adsorbate detector to ONLY get those which are in the central cell
#	#TODO: Probably only look at the O co-ords for the ads objs. Permuatations for the H are too messy
#	@unittest.skip("")
#	def testExpectedForTop(self):
#		expAds = self._loadUpperWaterAdsorbates()
#		self.assertTrue(False)
#
#













