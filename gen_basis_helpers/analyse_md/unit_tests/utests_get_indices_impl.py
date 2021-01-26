
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
		self.minAngle = 80
		self.maxAngle = 100

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
		self.minAngle = 10
		self.maxAngle = 20
		self.createTestObjs()
		expOutput = list()
		actOutput = self._runTestFunct()
		self.assertEqual(expOutput,actOutput)


#Simply mock the water index getter and stack 3 water molecules on top of each other
class TestGetTopWaterLayerIndicesA(unittest.TestCase):

	def setUp(self):
		self.lattParams = [10,10,10]
		self.lattAngles = [90,90,90]
		self.waterA = [ [5, 5, 5  ,"O"],
		                [5, 5, 5.5,"H"],
		                [5, 5, 4.5,"H"] ]

		self.waterB = [ [3, 3, 4  , "O"],
		                [3, 3, 3.5, "H"],
		                [3, 3, 3.5, "H"] ]

		self.waterC = [ [8, 8, 6  , "O"],
		                [8, 8, 6.5, "H"],
		                [8, 8, 5.5, "H"] ]

		self.waterIndices = [ [0,1,2], [3,4,5] , [6,7,8] ]

		self.maxLayerHeight = 0.5
		self.top = True

		self.createTestObjs()

	def createTestObjs(self):
		self.coordsA = list()
		self.coordsA.extend(self.waterA)
		self.coordsA.extend(self.waterB)
		self.coordsA.extend(self.waterC)

		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coordsA

		self.waterDetector = mock.Mock()
		self.waterDetector.getIndicesFromInpGeom.side_effect = lambda *args,**kwargs: self.waterIndices

		self.testObjA = tCode.GetTopWaterLayerIndices(waterDetector=self.waterDetector, maxHeightLayer=self.maxLayerHeight, top=self.top)

	def _runTestFunct(self):
		return self.testObjA.getIndicesFromInpGeom(self.cellA)

	def testExpRes_pbcsDontMatter(self):
		expIndices = [ [6,7,8] ]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices,actIndices)

	def testExpRes_pbcsImportant(self):
		self.waterC[0][2] -= self.lattParams[-1]
		self.createTestObjs()

		expIndices = [ [6,7,8] ]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices,actIndices)

	def testExpRes_twoLayersExpected(self):
		self.maxLayerHeight = 1.5
		self.createTestObjs()
		expIndices = [ [0,1,2], [6,7,8] ]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices,actIndices)

	def testExpRes_bottom(self):
		self.top = False
		self.createTestObjs()
		expIndices = [ [3,4,5] ]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices,actIndices)



