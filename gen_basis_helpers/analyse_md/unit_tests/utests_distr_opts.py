
import copy
import unittest

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.shared.plane_equations as planeEqnHelp

import gen_basis_helpers.analyse_md.distr_opt_objs as tCode


#Testing alternative constructor mostly

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



