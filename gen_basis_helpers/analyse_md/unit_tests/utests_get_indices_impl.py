
import unittest
import unittest.mock as mock

import plato_pylib.shared.unit_convs as uConvHelp
import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.get_indices_from_geom_impl as tCode

class TestGetSurfaceIndicesStandard(unittest.TestCase):

	def setUp(self):
		self.coordsA = [ [0,1,2,"Y"], [1,2,3,"X"], [2,3,4,"X"], [3,4,5,"X"] ]
		self.surfEles = ["X"]
		self.top = True
		self.bottom = True
		self.distTol = 1e-1
		self.nLayers = 1
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = self._createCellA()
		currKwargs = {"top":self.top, "bottom":self.bottom, "distTol":self.distTol, "nLayers":self.nLayers}
		self.testObjA = tCode.GetSurfaceIndicesFromGeomStandard(self.surfEles, **currKwargs)

	def _createCellA(self):
		cellA = uCellHelp.UnitCell(lattParams=[10,10,10], lattAngles=[90,90,90])
		cellA.cartCoords = self.coordsA
		return cellA

	def _runTestFunct(self):
		return self.testObjA.getIndicesFromInpGeom(self.cellA)

	def testTopAndBottomSimple(self):
		expIndices = [1,3]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices, actIndices)

	def testTwoLayers_top(self):
		self.nLayers = 2
		self.top, self.bottom = True, False
		self.createTestObjs()
		expIndices = [2,3]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices, actIndices)

	def testTwoLayers_topAndBottom(self):
		self.nLayers = 2
		self.surfEles = ["X","Y"]
		self.createTestObjs()
		expIndices = [0,1,2,3]
		actIndices = self._runTestFunct()
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

class TestGetIndicesForVaryTypesOfSurfAtom_waterBilayersSimple(unittest.TestCase):

	def setUp(self):
		self.cellA = _loadGeomTestVaryTypesOfSurfAtomA()
		self.surfEles = ["Mg"]
		self.surfDetector = tCode.GetSurfaceIndicesFromGeomStandard(["Mg"], top=True, bottom=False, distTol=2.0)
		self.maxWaterPlaneDist = 3.5*uConvHelp.ANG_TO_BOHR
		self.maxCloseDist = 2.5*uConvHelp.ANG_TO_BOHR
		self.maxHozDist = 1*uConvHelp.ANG_TO_BOHR
		self.createTestObjs()

	def createTestObjs(self):
		args = [self.surfDetector, self.surfEles, self.maxWaterPlaneDist, self.maxCloseDist, self.maxHozDist]
		self.testObjA = tCode.GetIndicesForVaryTypesOfSurfAtom_waterBilayerAdsorbedSimple(*args)

	def _runTestFunct(self):
		return self.testObjA.getIndicesFromInpGeom(self.cellA)

	def testExpCaseA(self):
		#Figured out by looking in VESTA...
		expIndices = {"free":[1,15,16], "close":[5,9,14], "far":[4,8,17]} 
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices, actIndices)


def _loadGeomTestVaryTypesOfSurfAtomA():
	lattVects = [ [18.19806, 0, 0],
	              [-9.099031, 15.75998, 0],
	              [0, 0, 63.85384] ]

	fractCoords = [ [-0.027588, -0.002342, 0.587224, "Mg"],
	                [0.079539, 0.217771, 0.667214, "Mg"],
	                [0.300281, -0.004826, 0.587257, "Mg"],
	                [0.636574, -0.001089, 0.586709, "Mg"],
	                [0.414056, 0.220636, 0.661339, "Mg"],
	                [0.746687, 0.219854, 0.684787, "Mg"],
	                [-0.033513, 0.326386, 0.587223, "Mg"],
	                [-0.031496, 0.665902, 0.588190, "Mg"],
	                [0.076254, 0.548272, 0.661229, "Mg"],
	                [0.076750, 0.884194, 0.680762, "Mg"],
	                [0.301607, 0.332359, 0.587536, "Mg"],
	                [0.303572, 0.662784, 0.587134, "Mg"],
	                [0.639368, 0.330226, 0.587770, "Mg"],
	                [0.633157, 0.660160, 0.588158, "Mg"],
	                [0.413768, 0.553394, 0.675142, "Mg"],
	                [0.415103, 0.887354, 0.662140, "Mg"],
	                [0.747242, 0.550726, 0.667078, "Mg"],
	                [0.750572, 0.888935, 0.661048, "Mg"],
	                [0.393481, 0.145364, 0.763020, "O"],
	                [0.692287, 0.144888, 0.750652, "O"],
	                [0.119679, 0.576804, 0.759279, "O"],
	                [0.112905, 0.868882, 0.745192, "O"],
	                [0.421279, 0.590298, 0.742401, "O"],
	                [0.703571, 0.863535, 0.760277, "O"],
	                [0.468150, 0.229714, 0.848836, "O"],
	                [0.230949, 0.661875, 0.833297, "O"],
	                [0.739291, 0.885323, 0.839228, "O"],
	                [0.730688, 0.174989, 0.834800, "O"],
	                [0.201793, 0.934010, 0.835390, "O"],
	                [0.475415, 0.585704, 0.837642, "O"],
	                [0.490333, 0.136773, 0.757922, "H"],
	                [0.389613, 0.205628, 0.739703, "H"],
	                [0.742192, 0.196035, 0.776193, "H"],
	                [0.688831, 0.039443, 0.752018, "H"],
	                [0.220775, 0.584559, 0.749718, "H"],
	                [0.035773, 0.494984, 0.741960, "H"],
	                [0.221143, 0.961875, 0.751903, "H"],
	                [0.119960, 0.770566, 0.749472, "H"],
	                [0.523562, 0.690244, 0.749931, "H"],
	                [0.421235, 0.504708, 0.758125, "H"],
	                [0.741290, 0.873071, 0.788386, "H"],
	                [0.787272, 0.866253, 0.742869, "H"],
	                [0.373530, 0.133857, 0.837371, "H"],
	                [0.477063, 0.321318, 0.834104, "H"],
	                [0.162038, 0.601881, 0.810493, "H"],
	                [0.237896, 0.768249, 0.831684, "H"],
	                [0.635529, 0.780647, 0.839557, "H"],
	                [0.714354, 0.973671, 0.840412, "H"],
	                [0.805502, 0.232381, 0.857004, "H"],
	                [0.638573, 0.191286, 0.839749, "H"],
	                [0.150681, 0.909521, 0.862192, "H"],
	                [0.120598, 0.918610, 0.815761, "H"],
	                [0.464175, 0.560176, 0.866494, "H"],
	                [0.388747, 0.610459, 0.832578, "H"] ]

	outCell = uCellHelp.UnitCell.fromLattVects(lattVects, fractCoords=fractCoords)
	return outCell
