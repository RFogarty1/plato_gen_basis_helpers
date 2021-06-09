
import unittest
import unittest.mock as mock

import plato_pylib.shared.unit_convs as uConvHelp
import plato_pylib.shared.ucell_class as uCellHelp


import gen_basis_helpers.shared.plane_equations as planeEqnHelp
import gen_basis_helpers.analyse_md.get_indices_from_geom_impl as tCode


class TestGetSortedIndicesBasedOnDistanceFromAtomIdx(unittest.TestCase):

	def setUp(self):
		self.inpIdx = 1
		self.inpIndices = [0,3]

		self.lattParamsA, self.lattAnglesA = [10,10,10], [90,90,90]
		self.cartCoords = [ [1,1,1,"A"],
		                    [1,1,9,"B"],
		                    [2,2,7,"D"],
		                    [1,1,7,"C"]
		                  ]

		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParamsA, lattAngles=self.lattAnglesA)
		self.cellA.cartCoords = self.cartCoords

	def _runTestFunct(self):
		return tCode.getSortedIndicesBasedOnDistanceFromInpAtomIdx(self.cellA, self.inpIdx, self.inpIndices)

	def testExpectedCaseA(self):
		expIndices = self.inpIndices
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices,actIndices)

	def testExpectedCaseB(self):
		self.inpIndices = [0,2,3]
		expIndices = [0,3,2]
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices,actIndices)

#Also handles getIndicesGroupedBySurfLayerForFirstNLayers
class TestGetSurfaceIndicesStandard(unittest.TestCase):

	def setUp(self):
		self.coordsA = [ [0,1,2,"Y"], [1,2,3,"X"], [2,3,4,"X"], [3,4,5,"X"] ]
		self.surfEles = ["X"]
		self.top = True
		self.bottom = True
		self.distTol = 1e-1
		self.nLayers = 1
		self.exitCleanlyWhenOutOfLayers = False
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

	#Below are tests for getIndicesGroupedBySurfLayerForFirstNLayers
	#This is so i can reuse the setup code

	def _runGetGroupedTestFunct(self):
		args = [self.cellA, self.surfEles]
		kwargs = {"top":self.top, "distTol":self.distTol, "nLayers":self.nLayers,
		          "exitCleanlyWhenOutOfLayers":self.exitCleanlyWhenOutOfLayers}
		return tCode.getIndicesGroupedBySurfLayerForFirstNLayers(*args, **kwargs)

	def testThreeLayersTop_groupedFunction(self):
		self.nLayers = 3
		self.top = True
		expIndices = [ [3], [2], [1] ]
		actIndices = self._runGetGroupedTestFunct()
		self.assertEqual(expIndices, actIndices)

	def testTwoLayersBot_groupedFunction(self):
		self.nLayers = 2
		self.top = False
		expIndices = [ [1], [2] ]
		actIndices = self._runGetGroupedTestFunct()
		self.assertEqual(expIndices,actIndices)

	def testErrorExitWhenTooManyLayersRequested_groupedFunction(self):
		self.nLayers = 50
		with self.assertRaises(ValueError):
			self._runGetGroupedTestFunct()

	def testCleanExitWhenTooManyLayersRequestd_groupedFunction(self):
		self.exitCleanlyWhenOutOfLayers = True
		self.top, self.surfEles = True, ["X","Y"]
		self.nLayers = 50
		expIndices = [ [3], [2], [1], [0] ]
		actIndices = self._runGetGroupedTestFunct()
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


class TestGetHydroxylIndicesStandard(unittest.TestCase):

	def setUp(self):
		self.minOH, self.maxOH = 0.01, 1.1*uConvHelp.ANG_TO_BOHR
		self.maxNebDist = 1.1*uConvHelp.ANG_TO_BOHR
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = _loadGeomSingleMgLayerHydroxylatedWithWaterMonolayer()
		kwargDict = {"minOH":self.minOH, "maxOH":self.maxOH, "maxNebDist":self.maxNebDist}
		self.testObjA = tCode.GetHydroxylMoleculeIndicesFromGeomStandard(**kwargDict)

	def _runTestFunct(self):
		return self.testObjA.getIndicesFromInpGeom(self.cellA)

	def testExpectedIndicesA(self):
		expIndices = _loadExpectedHydroxylIndices_singleMgLayerHydroxylatedWithWaterMonolayer()
		actIndices = self._runTestFunct()
		self.assertEqual( sorted(expIndices), sorted(actIndices) )

	def testNoneDetectedWhenMaxOHTooShort(self):
		self.maxOH = 0.5*uConvHelp.ANG_TO_BOHR
		self.createTestObjs()
		expIndices = list()
		actIndices = self._runTestFunct()
		self.assertEqual(expIndices, actIndices)

	def testNoneDetectedForXOH_geom(self):
		""" e.g. dont want to detect OH if its part of methanol/water etc. """
		self.maxNebDist = 1.5*uConvHelp.ANG_TO_BOHR
		self.createTestObjs()

		cartCoords = [ [2,2,2,"X"], [2,2,2+(1.4*uConvHelp.ANG_TO_BOHR),"O"], [2,2,2+((1.4+1.0)*uConvHelp.ANG_TO_BOHR),"H"],
		               [7,7,7,"O"], [7,7,7+(1.0*uConvHelp.ANG_TO_BOHR),"H"] ]
		self.cellA.cartCoords = cartCoords

		expIndices = [ [3,4] ]
		actIndices = self._runTestFunct()

		self.assertEqual( sorted(expIndices), sorted(actIndices) )

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



class TestGetWaterIndicesWithinDistOfInpIndices(unittest.TestCase):

	def setUp(self):
		self.cellA = _loadGeomTestVaryTypesOfSurfAtomA()
		self.inpIndicesA = [5,9]
		self.cutoffDist = 5 #bonds about 2.2 Angstrom i guess
		self.oxyInCutoff = True
		self.hyInCutoff = False
		self.waterDetector = None

	def _runTestFunct(self):
		args = [self.cellA, self.inpIndicesA, self.cutoffDist]
		kwargs = {"oxyInCutoff":self.oxyInCutoff, "hyInCutoff":self.hyInCutoff, "waterDetector":self.waterDetector}
		return tCode.getIndicesOfWaterWithinCutoffDistOfInpIndices(*args,**kwargs)

	def testExpectedForCaseA(self):
		#expVals found by dumping to *.cell file and looking in vesta
		expVals = [ [19,32,33], [21,36,37] ]
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)  

	def testExpectedUsingHydrogen_noneFound(self):
		self.oxyInCutoff, self.hyInCutoff = False,True
		expVals = []
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)

	#Tried to switch O/H coords at first, but failed since wastn detected as a "water" molecule
	def testExpectedUsingHydrogen_twoFound(self):
		self.oxyInCutoff, self.hyInCutoff = False, True
		self.cutoffDist = 6
		expVals = [ [19,32,33], [21,36,37] ]
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)

class TestGetHydroxylIndicesWithinDistOfInpIndices(unittest.TestCase):

	def setUp(self):
		self.hydroxyDetector = None
		self.inpIndices = [5]
		self.cutoffDist = 2.2*uConvHelp.ANG_TO_BOHR
		self.oxyInCutoff, self.hyInCutoff = True, True
		self.minDist = None
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = _loadGeomSingleMgLayerHydroxylatedWithWaterMonolayer()

	def _runTestFunct(self):
		args = [self.cellA, self.inpIndices, self.cutoffDist]
		kwargs = {"oxyInCutoff":self.oxyInCutoff, "hyInCutoff":self.hyInCutoff, "minDist":self.minDist}
		return tCode.getIndicesOfHydroxylWithinCutoffOfInpIndices(*args, **kwargs)

	def testExpCaseA(self):
		expVals = [ [13,14], [21,22], [25,26] ]
		actVals = self._runTestFunct()
		self.assertEqual( sorted(expVals), sorted(actVals) )

	def testNoneReturnedIfOxyInCutoffFalse(self):
		self.oxyInCutoff = False
		expVals = list()
		actVals = self._runTestFunct()
		self.assertEqual( sorted(expVals), sorted(actVals) )

	def testExpectedHyInCutoffWithLargerCutoffDist(self):
		self.oxyInCutoff = False
		self.cutoffDist = 2.8*uConvHelp.ANG_TO_BOHR
		expVals = [ [13,14], [21,22], [25,26] ]
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)

	def testExpectedMinAndMaxDistanceSpecified(self):
		self.hyInCutoff = False
		self.cutoffDist = 3.9*uConvHelp.ANG_TO_BOHR
		self.minDist = 2.4*uConvHelp.ANG_TO_BOHR
		self.createTestObjs()
		expVals = [ [9,10], [11,12], [19,20] ]
		actVals = self._runTestFunct()
		self.assertEqual( sorted(expVals), sorted(actVals) )

class TestGetIndicesOfBilayerClosestToSurface(unittest.TestCase):

	def setUp(self):
		self.cellA = _loadGeomTestVaryTypesOfSurfAtomA()
		self.surfaceDetector = tCode.GetSurfaceIndicesFromGeomStandard(["Mg"], top=True, bottom=False, distTol=2, nLayers=1)
		self.maxBilayerThickness = 1.5
		self.waterDetector = None
		self.expNumberWater = None
		self.planeEqn = None

	def _runTestFunct(self):
		args = [self.cellA, self.surfaceDetector]
		kwargs = {"maxBilayerThickness":self.maxBilayerThickness, "waterDetector":self.waterDetector,
		          "expNumberWater":self.expNumberWater, "planeEqn":self.planeEqn}
		return tCode.getIndicesOfWaterBilayerClosestToSurface(*args, **kwargs)

	def testExpectedIndicesA(self):
		expVals = [ [18, 30, 31], [19, 32, 33], [20, 34, 35],
		            [21, 36, 37], [22, 38, 39], [23, 40, 41] ]
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)

	def testRaisesWhenUnexpectedNumberReturned(self):
		self.expNumberWater = 20
		with self.assertRaises(AssertionError):
			self._runTestFunct()

	def testExpectedForZeroWaterPresent(self):
		self.cellA.cartCoords = [x for x in self.cellA.cartCoords if x[-1].upper()=="MG"]
		expVals = []
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)

	def testRaisesIfBothSurfaceDetectorAndPlaneEqnSet(self):
		self.planeEqn = 6
		with self.assertRaises(ValueError):
			self._runTestFunct()

	def testRaisesIfNeitherSurfaceDetectorNorPlaneEqnSet(self):
		self.surfaceDetector=None
		with self.assertRaises(ValueError):
			self._runTestFunct()

	def testExpectedResultUsingPlaneEqnInput(self):
		self.surfaceDetector=None
		self.planeEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,63) #This lies at the top of the cell; hence the nearest layer is the layer furthest from the surface
		expSecondLayer = [ [24, 42, 43], [25, 44, 45], [26, 46, 47],
		                   [27, 48, 49], [28, 50, 51], [29, 52, 53] ]
		actSecondLayer = self._runTestFunct()
		self.assertEqual(expSecondLayer,actSecondLayer)		


class TestGetNBilayersClosestToSurface(unittest.TestCase):

	def setUp(self):
		self.cellA = _loadGeomTestVaryTypesOfSurfAtomA()
		self.surfaceDetector = tCode.GetSurfaceIndicesFromGeomStandard(["Mg"], top=True, bottom=False, distTol=2, nLayers=1)
		self.maxBilayerThickness = 1.5
		self.waterDetector = None
		self.expWaterPerLayer = None
		self.maxNLayers = None
		self.planeEqn = None

	def _runTestFunct(self):
		args = [self.cellA, self.surfaceDetector]
		kwargs = {"maxBilayerThickness":self.maxBilayerThickness, "waterDetector":self.waterDetector,
		          "expWaterPerLayer":self.expWaterPerLayer, "maxNLayers":self.maxNLayers, "planeEqn":self.planeEqn}
		return tCode.getIndicesOfWaterBilayersStartingClosestToSurface(*args,**kwargs)

	def _loadExpFirstAndSecondLayerDefaultSettings(self):
		expFirstLayer  = [ [18, 30, 31], [19, 32, 33], [20, 34, 35],
		                   [21, 36, 37], [22, 38, 39], [23, 40, 41] ]

		expSecondLayer = [ [24, 42, 43], [25, 44, 45], [26, 46, 47],
		                   [27, 48, 49], [28, 50, 51], [29, 52, 53] ]
		return [expFirstLayer, expSecondLayer]

	def testExpectedIndicesA(self):
		expVals = self._loadExpFirstAndSecondLayerDefaultSettings()
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)

	def testExpectedFirstLayerOnly(self):
		self.maxNLayers = 1
		expVals = [self._loadExpFirstAndSecondLayerDefaultSettings()[0]]
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)

	def testRaisesForUnexpectedNumberOfWaterInABilayer(self):
		self.expWaterPerLayer = 18
		with self.assertRaises(AssertionError):
			self._runTestFunct()

	def testExpectedUsingPlaneEqnInput(self):
		self.surfaceDetector=None
		self.planeEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,63) #This lies at the top of the cell; hence the nearest layer is the layer furthest from the surface
		expFirstLayer, expSecondLayer = self._loadExpFirstAndSecondLayerDefaultSettings()
		expVals = [expSecondLayer, expFirstLayer]
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)


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


def _loadGeomSingleMgLayerHydroxylatedWithWaterMonolayer():
	lattVectors = [ [x*uConvHelp.BOHR_TO_ANG for x in [18.19806 ,        0,        0]],
	                [x*uConvHelp.BOHR_TO_ANG for x in [-9.099031, 15.75998,        0]],
	                [x*uConvHelp.BOHR_TO_ANG for x in [-9.099031, 15.75998, 63.97856]] ]

	cartCoords = [ [  0.1214542359, 1.8196267123, 22.1955723719,  "Mg"],
	               [  3.3284047979, 1.8212532892, 22.1914370011,  "Mg"],
	               [  6.5383078083, 1.8262059327, 22.1816696526,  "Mg"],
	               [ -1.4864883933, 4.6011105484, 22.1915290289,  "Mg"],
	               [ -3.0916117901, 7.3860392044, 22.1818495864,  "Mg"],
	               [  1.7233630404, 4.6060280692, 22.1817849271,  "Mg"],
	               [  0.1214463163, 7.3795391972, 22.1956747834,  "Mg"],
	               [  4.9363577945, 4.5995439582, 22.1955293125,  "Mg"],
	               [  3.3283650735, 7.3810176883, 22.1914579493,  "Mg"],
	               [   4.9357408756, 6.4558465066, 23.1428645602, "O"],
	               [   4.9317129760, 6.4578090856, 24.1261352568, "H"],
	               [  -1.4875861102, 6.4563188250, 23.0955659872, "O"],
	               [  -1.5390639112, 6.4602746494, 24.0812622605, "H"],
	               [   1.7273417794, 6.4542934871, 23.1366283102, "O"],
	               [   1.7165806377, 6.4529634491, 24.1185714081, "H"],
	               [   8.1422929808, 0.8964480936, 23.0955843628, "O"],
	               [   8.0912046708, 0.9003698378, 24.0812871965, "H"],
	               [   6.5424373543, 3.6745332540, 23.1365246049, "O"],
	               [   6.5316883407, 3.6730287745, 24.1184459409, "H"],
	               [   1.7273456379, 0.8945823024, 23.1366219878, "O"],
	               [   1.7166759222, 0.8932007172, 24.1185678419, "H"],
	               [   0.1208818476, 3.6756884827, 23.1430080822, "O"],
	               [   0.1171114663, 3.6776744068, 24.1262808732, "H"],
	               [   4.9356624802, 0.8959088297, 23.1428943007, "O"],
	               [   4.9319239561, 0.8978710944, 24.1261676150, "H"],
	               [   3.3274315554, 3.6763986244, 23.0956350496, "O"],
	               [   3.2759433699, 3.6802622178, 24.0813194442, "H"],
	               [  -1.8874432309, 6.4767520522, 25.9624287592, "O"],
	               [  -2.4829221472, 5.6946685286, 26.0560040799, "H"],
	               [  -2.4775581249, 7.2632244931, 26.0560011187, "H"],
	               [   1.1522083310, 6.4580466268, 26.2169360954, "O"],
	               [   1.3331680081, 6.4732780923, 27.1807979360, "H"],
	               [   0.1621470753, 6.4589103759, 26.1563564863, "H"],
	               [   7.7424923019, 0.9169105569, 25.9626207374, "O"],
	               [   7.1470233194, 0.1348774287, 26.0563652099, "H"],
	               [   7.1524611077, 1.7034217508, 26.0559810826, "H"],
	               [   5.9671298549, 3.6781874135, 26.2169632216, "O"],
	               [   6.1480902956, 3.6928928451, 27.1808290255, "H"],
	               [   4.9770596881, 3.6791149719, 26.1564070752, "H"],
	               [   1.1521037718, 0.8981902086, 26.2169285841, "O"],
	               [   1.3330809663, 0.9133915044, 27.1807843755, "H"],
	               [   0.1620338414, 0.8992287287, 26.1563771918, "H"],
	               [   2.9275561518, 3.6967973390, 25.9627056626, "O"],
	               [   2.3321023090, 2.9147150091, 26.0564069628, "H"],
	               [   2.3374510526, 4.4832542552, 26.0562773955, "H"] ]

	outCell = uCellHelp.UnitCell.fromLattVects(lattVectors)
	outCell.cartCoords = cartCoords
	outCell.convAngToBohr()
	return outCell

def _loadExpectedHydroxylIndices_singleMgLayerHydroxylatedWithWaterMonolayer():
	outIndices = [ [9,10] , [11,12], [13,14], [15,16], [17,18],
	               [19,20], [21,22], [23,24], [25,26] ]
	return outIndices

