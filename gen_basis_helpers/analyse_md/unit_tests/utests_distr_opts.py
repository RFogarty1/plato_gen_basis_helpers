
import copy
import itertools as it
import unittest

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.shared.plane_equations as planeEqnHelp

import gen_basis_helpers.analyse_md.binned_res as binResHelp
import gen_basis_helpers.analyse_md.distr_opt_objs as tCode



class TestWaterOrientationOpts(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.cartCoords = [ [0,0,0,"O"], [1,1,1,"H"], [2,2,2,"H"],
		                    [3,3,3,"O"], [4,4,4,"H"], [5,5,5,"H"],
		                    [8,8,8,"X"], [9,9,9,"Y"] ]

		#Options for the actual opts obj
		self.checkEdges = True
		self.binEdges = [-90,0,90]
		self.oxyIndices = [0,3]
		self.hyIndices = [ [1,2], [4,5] ]
		self.angleType = "roll"
		self.primaryIdxType = "O"

		self.createTestObjs()

	def createTestObjs(self):
		#Create the cell
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

		#Create other stuff
		self.binObjA = binResHelp.BinnedResultsStandard.fromBinEdges(self.binEdges)
		args = [self.binObjA, self.oxyIndices, self.hyIndices]
		kwargs = {"angleType":self.angleType, "checkEdges":self.checkEdges}
		self.testObjA = tCode.WaterOrientationOptions(*args,**kwargs)

	def testExpectedDomainsUponChanging(self):
		angleTypes = ["roll", "pitch", "azimuth"]
		expDomains = [ [-90,90], [-90,90], [-180,180] ]
		for aType,expDomain in it.zip_longest(angleTypes,expDomains):
			self.testObjA.angleType = aType
			actDomain = self.testObjA.domain
			self.assertAlmostEqual(expDomain[0], actDomain[0])
			self.assertAlmostEqual(expDomain[1], actDomain[1])

	def testRaisesIfCheckEdgesAndWeSetToAngleTypeWithLowerDomain(self):
		#Check error on initiation
		self.binEdges = [0,50,92]
		with self.assertRaises(ValueError):
			self.createTestObjs()

		#Check error when changin from azimuth (large domain) to pitch (smaller domain)
		self.angleType="azimuth"
		self.createTestObjs()
		with self.assertRaises(ValueError):
			self.testObjA.angleType = "roll"

	def testAltConstructor(self):
		expObj = self.testObjA

		waterIndices = [ [1,2,0],[4,5,3] ]
		currArgs = [self.binObjA, waterIndices, self.cellA]
		currKwargs = {"primaryIdxType":self.primaryIdxType, "checkEdges":self.checkEdges, "angleType":self.angleType}
		actObj = tCode.WaterOrientationOptions.fromWaterIndicesAndGeom(*currArgs, **currKwargs)

		self.assertEqual(expObj,actObj)


#Testing alternative constructor mostly
class TestWaterMinDistOpts(unittest.TestCase):

	def setUp(self):
		#Geom for testing the alternative initializer
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.cartCoords = [ [0,0,0,"O"], [1,1,1,"H"], [2,2,2,"H"],
		                    [3,3,3,"O"], [4,4,4,"H"], [5,5,5,"H"],
		                    [8,8,8,"X"], [9,9,9,"Y"] ]

		#Options for the actual opts obj
		self.oxyIndices = [0,3]
		self.hyIndices = [[1,2],[4,5]]
		self.toIndices = [6,7]
		self.primaryIdxType = "O"
		self.minDistType = "all"
		self.binResObj = 5 #Int has functionally similar equality method to the binRes obj so...

		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams,lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

		#Create test obj
		currArgs = [self.binResObj, self.oxyIndices, self.hyIndices, self.toIndices]
		currKwargs = {"primaryIdxType":self.primaryIdxType, "minDistType":self.minDistType}
		self.testObj = tCode.WaterMinDistOptions(*currArgs, **currKwargs)

	def testAltConstructorCaseA(self):
		waterIndices = [ [0,1,2], [3,4,5] ]
		expObj = self.testObj
		currArgs = [self.binResObj, waterIndices, self.toIndices, self.cellA]
		currKwargs = {"primaryIdxType":self.primaryIdxType, "minDistType":self.minDistType}
		actObj = tCode.WaterMinDistOptions.fromWaterIndicesAndGeom(*currArgs, **currKwargs)
		self.assertEqual(expObj, actObj)

class TestWaterMinPlanarDistOpts(unittest.TestCase):

	def setUp(self):
		#Geom for testing the alternative initializer
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.cartCoords = [ [0,0,0,"O"], [1,1,1,"H"], [2,2,2,"H"],
		                    [3,3,3,"O"], [4,4,4,"H"], [5,5,5,"H"] ]

		#Options for the actual opts obj
		self.oxyIndices = [0,3]
		self.hyIndices = [[1,2],[4,5]]
		self.planeEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,5)
		self.primaryIdxType = "O"
		self.minDistType = "all"
		self.binResObj = 5 #Int has functionally similar equality method to the binRes obj so...

		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams,lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

		#Create test obj
		currArgs = [self.binResObj, self.oxyIndices, self.hyIndices]
		currKwargs = {"planeEqn":self.planeEqn, "primaryIdxType":self.primaryIdxType, "minDistType":self.minDistType}
		self.testObj = tCode.WaterMinPlanarDistOptions(*currArgs, **currKwargs)
		
	def testAltConstructorCaseA(self):
		waterIndices = [ [0,1,2], [3,4,5] ]
		expObj = self.testObj
		currKwargs = {"planeEqn":self.planeEqn,"primaryIdxType":self.primaryIdxType, "minDistType":self.minDistType}
		actObj = tCode.WaterMinPlanarDistOptions.fromWaterIndicesAndGeom(self.binResObj, waterIndices, self.cellA, **currKwargs)
		self.assertEqual(expObj,actObj)


#Mainly just to test the alternative constructor (using water indices + geom) works
class TestWaterPlanarDistOpts(unittest.TestCase):

	def setUp(self):
		#Geom for testing the alternative initializer
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.cartCoords = [ [0,0,0,"O"], [1,1,1,"H"], [2,2,2,"H"],
		                    [3,3,3,"O"], [4,4,4,"H"], [5,5,5,"H"] ]

		#Options for the actual opts obj
		self.oxyIndices = [3,0]
		self.hyIndices = [[4,5],[1,2]]
		self.planeEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,5)
		self.primaryIdxType = "O"
		self.binResObj = 5 #Int has functionally similar equality method to the binRes obj so...

		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams,lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

		currArgs = [self.binResObj, self.oxyIndices, self.hyIndices]
		currKwargs = {"planeEqn":self.planeEqn, "primaryIdxType":self.primaryIdxType}
		self.testObj = tCode.WaterPlanarDistOptions(*currArgs, **currKwargs)

	def testCmpEqual(self):
		objA = copy.deepcopy(self.testObj)
		self.createTestObjs()
		objB = self.testObj
		self.assertEqual(objA, objB)

	def testCmpUnequal_diffPlaneEqn(self):
		objA = copy.deepcopy(self.testObj)
		self.planeEqn = planeEqnHelp.ThreeDimPlaneEquation(1,0,1,3)
		self.createTestObjs()
		objB = self.testObj
		self.assertNotEqual(objA, objB)

	def testCmpUnequal_diffHyIndices(self):
		objA = copy.deepcopy(self.testObj)
		self.hyIndices = [[8,9],[1,6]]
		self.createTestObjs()
		objB = self.testObj
		self.assertNotEqual(objA,objB)

	def testInitializer_fromGeomAndWaterIndices(self):
		waterIndices = [ [5,3,4], [1,2,0] ] #Oxygen in different slot for each
		currKwargs = {"planeEqn":self.planeEqn, "primaryIdxType":self.primaryIdxType}
		expObj = tCode.WaterPlanarDistOptions(self.binResObj, [3,0], [[5,4],[1,2]], **currKwargs)
		actObj = tCode.WaterPlanarDistOptions.fromWaterIndicesAndGeom(self.binResObj, waterIndices, self.cellA, **currKwargs)

		self.assertEqual(expObj,actObj)



