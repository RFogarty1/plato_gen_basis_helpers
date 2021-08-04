
import copy
import itertools as it
import unittest

import numpy as np

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.shared.plane_equations as planeEqnHelp

import gen_basis_helpers.analyse_md.atom_combo_core as tCode



class TestSparseMatrixCalculatorA(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.cartCoords = [ [2,2,4,"Mg"],
		                    [2,2,1,"Mg"],
		                    [2,2,9,"Mg"],
		                    [2,2,8,"Mg"] ]

		self.distIndicesFrom = [0,2]
		self.distIndicesTo = [1,3]

		self.planarEqnA = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,2)
		self.planarEqnB = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,2)

		self.planarIndicesA = [0,2]
		self.planarIndicesB = [3]

		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

		#Create individual populators
		self.distsPopulatorA = tCode._DistMatrixPopulator( self.distIndicesFrom, self.distIndicesTo )
		self.planarPopulatorA = tCode._PlanarDistMatrixPopulator( self.planarIndicesA, self.planarEqnA )
		self.planarPopulatorB = tCode._PlanarDistMatrixPopulator( self.planarIndicesB, self.planarEqnB )

		#Create the overall calculator
		currPopulators = [self.distsPopulatorA, self.planarPopulatorA, self.planarPopulatorB]
		self.testObj = tCode._SparseMatrixCalculatorStandard(currPopulators)

	def _runTestFunct(self):
		self.testObj.calcMatricesForGeom(self.cellA)

	def _loadExpectedDictA(self):
		#Load the distance matrix
		distMatrix = np.empty( (4,4) )
		distMatrix[:] = np.nan
		distMatrix[0][1], distMatrix[1][0] = 3,3
		distMatrix[0][3], distMatrix[3][0] = 4,4
		distMatrix[2][1], distMatrix[1][2] = 2,2
		distMatrix[2][3], distMatrix[3][2] = 1,1

		#Load the plane equations
		uniquePlaneEqns = [self.planarEqnA]
		planarMatrix = np.array( [2, np.nan, 3, 4] ) 

		outDict = dict()
		outDict["uniquePlaneEquations"] = uniquePlaneEqns
		outDict["planarDists"] = [ planarMatrix ] 
		outDict["distMatrix"] = distMatrix
		return outDict

	def _testExpectedAndActualDictsEqual(self, expDict, actDict):
		#Distance matrices
		expDistMatrix, actDistMatrix = expDict["distMatrix"], actDict["distMatrix"]
		self.assertTrue( np.allclose(expDistMatrix,actDistMatrix,equal_nan=True) )

		#plane equations
		for expPlaneEqn, actPlaneEqn in it.zip_longest(expDict["uniquePlaneEquations"],actDict["uniquePlaneEquations"]):
			self.assertEqual(expPlaneEqn,actPlaneEqn)

		#Planar distance matrices
		for expMatrix, actMatrix in it.zip_longest(expDict["planarDists"], actDict["planarDists"]):
			self.assertTrue( np.allclose(expMatrix,actMatrix,equal_nan=True) )


	#Should be pretty comprehensive; multiple of same type and also two different types simultaneously
	def testExpectedCaseA(self):
		expDict = self._loadExpectedDictA()
		self._runTestFunct()
		actDict = self.testObj.outDict
		self._testExpectedAndActualDictsEqual(expDict, actDict)


	def testExpectedCase_matricesAlreadyPopulated(self):
		""" This allows us to check that new matrices WILL overwrite the previous ones """
		self._runTestFunct()
		self.testObj.outDict["distMatrix"][0][1] += 2
		expDict = self._loadExpectedDictA()
		self._runTestFunct()
		actDict = self.testObj.outDict
		self._testExpectedAndActualDictsEqual(expDict, actDict)



class TestSparseMatrixCalculatorEquality(unittest.TestCase):

	def setUp(self):
		self.populators = [4,7,2] #Better than mocks here; since these have effectively similar criteria for comparison to actual populators
		self.createTestObjs()

	def createTestObjs(self):
		self.testObj = tCode._SparseMatrixCalculatorStandard(self.populators)

	def testEqualCmp(self):
		objA = copy.deepcopy(self.testObj)
		self.createTestObjs()
		objB = self.testObj
		self.assertEqual(objA,objB)

	def testUnequalCmp_diffValPopulators(self):
		objA = copy.deepcopy(self.testObj)
		self.populators[-1] += 2
		self.createTestObjs()
		objB = self.testObj
		self.assertNotEqual(objA, objB)

	def testUnequalCmp_diffLenPopulators(self):
		objA = copy.deepcopy(self.testObj)
		self.populators.append(4)
		self.createTestObjs()
		objB = self.testObj
		self.assertNotEqual(objA,objB)


class TestPlanarDistMatrixPopulatorEquality(unittest.TestCase):

	def setUp(self):
		self.indices = [0,3]
		self.planeEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,2)
		self.objLevel = 2
		self.createTestObjs()

	def createTestObjs(self):
		self.testObj = tCode._PlanarDistMatrixPopulator(self.indices, self.planeEqn, level=self.objLevel)

	def testEqualCmp(self):
		objA = copy.deepcopy(self.testObj)
		self.createTestObjs()
		objB = self.testObj
		self.assertEqual(objA, objB)

	def testUnequalCmp_diffIndices(self):
		objA = copy.deepcopy(self.testObj)
		self.indices[1] += 2
		self.createTestObjs()
		objB = self.testObj
		self.assertNotEqual(objA,objB)

	def testUnequalCmp_diffPlaneEqn(self):
		objA = copy.deepcopy(self.testObj)
		self.planeEqn = planeEqnHelp.ThreeDimPlaneEquation(1,0,1,3)
		self.createTestObjs()
		objB = self.testObj
		self.assertNotEqual(objA,objB)

	def testUnequalCmp_diffTypes(self):
		objA = self.testObj
		objB = 5
		self.assertNotEqual(objA, objB)
		self.assertNotEqual(objB, objA)


class TestDistMatrixPopulatorEquality(unittest.TestCase):

	def setUp(self):
		self.fromIndices = [0,3]
		self.toIndices = [1,2]
		self.level = 4
		self.createTestObjs()

	def createTestObjs(self):
		self.testObj = tCode._DistMatrixPopulator(self.fromIndices, self.toIndices, level=self.level)

	def testEqualCmp(self):
		objA = copy.deepcopy(self.testObj)
		self.createTestObjs()
		objB = self.testObj
		self.assertEqual(objA,objB)

	def testUnequalCmp_diffToIndices(self):
		objA = copy.deepcopy(self.testObj)
		self.toIndices[1] += 2
		self.createTestObjs()
		objB = self.testObj
		self.assertNotEqual(objA, objB)

	def testUnequalCmp_diffTypes(self):
		objA = self.testObj
		objB = 5
		self.assertNotEqual(objA, objB)
		self.assertNotEqual(objB, objA)


class TestPlanarDistMatrixPopulator(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.cartCoords = [ [2,2,4,"Mg"],
		                    [2,2,1,"Mg"],
		                    [2,2,9,"Mg"],
		                    [2,2,8,"Mg"] ]

		self.indices = [0,2]
		self.planeEqnA = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,2)
		self.objLevel = 0
		self.callLevel = 0

		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords
		self.testObj = tCode._PlanarDistMatrixPopulator(self.indices, self.planeEqnA, level=self.objLevel)
		self.outDict = dict()

	def _runTestFunct(self):
		self.testObj.populateMatrices(self.cellA, self.outDict, self.callLevel)

	def testExpected_emptyDict(self):
		expMatrix = np.empty( (4) )
		expMatrix[:] = np.nan
		expMatrix[0] = 2
		expMatrix[2] = 3
		expUniqueEquations = [self.planeEqnA]

		self._runTestFunct()
		actMatrix = self.outDict["planarDists"][0]
		actUniqueEquations = self.outDict["uniquePlaneEquations"]

		self.assertEqual(expUniqueEquations,actUniqueEquations)
		self.assertTrue( np.allclose(expMatrix,actMatrix,equal_nan=True) )

	def testExpected_firstPlaneEqnPopulated(self):
		prevPlaneEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,5)
		prevPlanarDistMatrix = np.empty( (4) )
		prevPlanarDistMatrix[:] = np.nan

		self.outDict["uniquePlaneEquations"] = [ prevPlaneEqn ]
		self.outDict["planarDists"] = [ prevPlanarDistMatrix ]  

		#
		expNewMatrix = np.empty( (4) )
		expNewMatrix[:] = np.nan
		expNewMatrix[0], expNewMatrix[2] = 2,3

		expMatrices = [prevPlanarDistMatrix, expNewMatrix]
		expUniquePlaneEqns = [prevPlaneEqn, self.planeEqnA]

		self._runTestFunct()

		self.assertEqual(expUniquePlaneEqns, self.outDict["uniquePlaneEquations"])
		for expMatrix,actMatrix in it.zip_longest(expMatrices, self.outDict["planarDists"]):
			self.assertTrue( np.allclose(expMatrix,actMatrix,equal_nan=True) )

	def testExpected_partOfPlaneEqnMatrix(self):
		#Load initial environment
		initMatrix = np.empty( (4) )
		initMatrix[0] = 2
		self.outDict["uniquePlaneEquations"] = [ self.planeEqnA ]
		self.outDict["planarDists"] = [initMatrix]

		#
		expNewMatrix = [2, np.nan, 3, np.nan]
		expMatrices = [ expNewMatrix ]
		expUniquePlaneEquations = [ self.planeEqnA ]
		
		self._runTestFunct()
		actUniquePlaneEquations = self.outDict["uniquePlaneEquations"]
		actMatrices = self.outDict["planarDists"]

		self.assertEqual(expUniquePlaneEquations, self.outDict["uniquePlaneEquations"])
		for expMatrix,actMatrix in it.zip_longest(expMatrices,actMatrices):
			self.assertTrue( np.allclose(expMatrix,actMatrix,equal_nan=True) )


class TestDistMatrixPopulator(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.cartCoords = [ [2,2,4,"Mg"],
		                    [2,2,1,"Mg"],
		                    [2,2,9,"Mg"],
		                    [2,2,8,"Mg"] ]

		self.fromIndices = [0,3]
		self.toIndices = [1,2]
		self.objLevel = 0
		self.callLevel = 0

		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords
		self.testObj = tCode._DistMatrixPopulator(self.fromIndices, self.toIndices, level=self.objLevel)
		self.outDict = dict()

	def _runTestFunct(self):
		self.testObj.populateMatrices(self.cellA, self.outDict, self.callLevel)

	def testExpectedMatrixMissing(self):
		#setup expected matrix
		expMatrix = np.empty( (4,4) )
		expMatrix[:] = np.nan
		expMatrix[0][1], expMatrix[1][0] = 3, 3
		expMatrix[0][2], expMatrix[2][0] = 5, 5

		expMatrix[3][1], expMatrix[1][3] = 3, 3
		expMatrix[3][2], expMatrix[2][3] = 1, 1

		#Test actual matrix matches
		self._runTestFunct()
		actMatrix = self.outDict["distMatrix"]

		self.assertTrue( np.allclose(expMatrix,actMatrix,equal_nan=True) )

	def testExpectedMatrixPartiallyPresent(self):
		#Partially populate starting matrix
		self.outDict["distMatrix"] = np.empty( (4,4) )
		self.outDict["distMatrix"][:] = np.nan
		self.outDict["distMatrix"][3][2], self.outDict["distMatrix"][2][3] = 1,1
		self.outDict["distMatrix"][0][2], self.outDict["distMatrix"][2][0] = 5,5

		#Setup expected matrix
		expMatrix = np.empty( (4,4) )
		expMatrix[:] = np.nan
		expMatrix[0][1], expMatrix[1][0] = 3, 3
		expMatrix[0][2], expMatrix[2][0] = 5, 5

		expMatrix[3][1], expMatrix[1][3] = 3, 3
		expMatrix[3][2], expMatrix[2][3] = 1, 1

		#Test actual matrix matches
		self._runTestFunct()
		actMatrix = self.outDict["distMatrix"]

		self.assertTrue( np.allclose(expMatrix, actMatrix, equal_nan=True) )

	def testNothingHappensWhenLevelWrong(self):
		self.callLevel += 1
		self._runTestFunct()
		with self.assertRaises(KeyError):
			unused = self.outDict["distMatrix"]





