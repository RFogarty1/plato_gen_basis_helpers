
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.get_indices_from_geom_impl as tCode

class TestGetSurfaceIndicesStandard(unittest.TestCase):

	def setUp(self):
		self.coordsA = [ [0,1,2,"Y"], [1,2,3,"X"], [2,3,4,"X"], [3,4,5,"X"] ]
		self.surfEles = ["X"]
		self.top = True
		self.bottom = True
		self.distTol = 1e-1
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = self._createCellA()
		currKwargs = {"top":self.top, "bottom":self.bottom, "distTol":self.distTol}
		self.testObjA = tCode.GetSurfaceIndicesFromGeomStandard(self.surfEles, **currKwargs)

	def _createCellA(self):
		cellA = uCellHelp.UnitCell(lattParams=[10,10,10], lattAngles=[90,90,90])
		cellA.cartCoords = self.coordsA
		return cellA

	def testTopAndBottomSimple(self):
		expIndices = [1,3]
		actIndices = self.testObjA.getIndicesFromInpGeom(self.cellA)
		self.assertEqual(expIndices, actIndices)


