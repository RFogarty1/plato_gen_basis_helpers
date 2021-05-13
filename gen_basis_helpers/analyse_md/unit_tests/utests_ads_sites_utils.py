
import copy
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.ads_sites_impl as adsSiteImplHelp
import gen_basis_helpers.analyse_md.ads_sites_utils as tCode

class TestAddAdsSitesStandard(unittest.TestCase):

	def setUp(self):
		self.topIndices = [1,2] 
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = loadTestCellA()
		self.topSites = [adsSiteImplHelp.TopStandard(idx,siteName="top") for idx in self.topIndices]
		self.testObjA = tCode.AddAdsSitesToGeomsStandard(self.topSites)

	def _runTestFunct(self):
		args = self.cellA
		return self.testObjA.addToUnitCell(self.cellA)

	def testAddToUnitCellSimple(self):
		expCoords = [ [4,4,6,"X"],
		              [5,5,7,"Y"],
		              [6,6,8,"Z"],
		              [5,5,7,"top"],
		              [6,6,8,"top"] ]
		self._runTestFunct()
		actCoords = self.cellA.cartCoords
		testVal = areExpAndActCoordsEqualUsingTemplateCell(expCoords, actCoords, self.cellA)
		self.assertTrue(testVal)


def areExpAndActCoordsEqualUsingTemplateCell(expCoords, actCoords, templCell):
	expCell, actCell = copy.deepcopy(templCell), copy.deepcopy(templCell)
	expCell.cartCoords = expCoords
	actCell.cartCoords = actCoords
	if expCell == actCell:
		return True
	else:
		return False

def loadTestCellA():
	lattParams, lattAngles = [10,10,10], [90,90,90]
	cartCoords = [ [4,4,6,"X"],
	               [5,5,7,"Y"],
	               [6,6,8,"Z"] ]

	outCell = uCellHelp.UnitCell(lattParams=lattParams, lattAngles=lattAngles)
	outCell.cartCoords = cartCoords
	return outCell


