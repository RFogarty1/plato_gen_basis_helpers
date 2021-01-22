
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



class TestGetWaterMoleculeIndicesStandard(unittest.TestCase):

	def setUp(self):
		self.lattParams = [10,10,10]
		self.lattAngles = [90,90,90]
		self.coordsA = [ [5,5,5,"O"],
		                 [5,6,6,"H"],
		                 [5,4,6,"H"],
		                 [1,1,8,"O"],
		                 [9,9,1,"O"] ]

		#Parameters for the "molecule finder"
		self.minOH = 0.5
		self.maxOH = 2.0
		self.minAngle = 40
		self.maxAngle = 50

		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=[10,10,10], lattAngles=[90,90,90])
		self.cellA.cartCoords = self.coordsA

		kwargDict = {"minOH":self.minOH, "maxOH":self.maxOH, "minAngle":self.minAngle, "maxAngle":self.maxAngle}
		self.testObjA = tCode.GetWaterMoleculeIndicesFromGeomStandard(**kwargDict)

	def _runTestFunct(self):
		return self.testObjA.getIndicesFromInpGeom(self.cellA)

	def testExpected_pbcsNotImportant(self):
		expOutput = [ [0,1,2] ]
		actOutput = self._runTestFunct()
		self.assertEqual(expOutput, actOutput)

	def testExpected_pbcsImportant(self):
		self.coordsA[0][-2] += self.lattParams[-1]
		self.createTestObjs()
		expOutput = [ [0,1,2] ]
		actOutput = self._runTestFunct()
		self.assertEqual(expOutput, actOutput)

	def testExpected_angleWrong(self):
		self.minAngle = 90
		self.maxAngle = 150
		self.createTestObjs()
		expOutput = list()
		actOutput = self._runTestFunct()
		self.assertEqual(expOutput,actOutput)








