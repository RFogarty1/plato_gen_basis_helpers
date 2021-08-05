
import copy
import itertools as it
import unittest

import numpy as np

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.shared.plane_equations as planeEqnHelp

import gen_basis_helpers.analyse_md.atom_combo_core as atomComboCoreHelp

import gen_basis_helpers.analyse_md.atom_combo_populators as tCode


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

@unittest.skip("")
class TestDiscHBondCounterWithDistFilterPopulator(unittest.TestCase):

	def setUp(self):
		#1) All geometric parameters for testing
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		#roll 90,pitch=45, OH len~1, HOH angle 104.5. Then just added translation vectors
		self.waterACoords = [ [0,0,0,"O"], [-0.13,0,0.99,"H"], [0.99,0,-0.13,"H"] ]
		#translation = [2,0,0]; no rotations
		self.waterBCoords = [ [2,0,0,"O"], [2.61, 0.79, 0, "H"], [2.61,-0.79,0,"H"] ]
		#translation = [2+(2*0.61), 2*0.79, 0]
		self.waterCCoords = [ [3.22, 1.58, 0, "O"],  [3.83, 2.37, 0, "H"], [3.83, 0.79, 0, "H"] ]
		#translation = [2+(2*0.61), -2*0.79,0]
		self.waterDCoords = [ [3.22, -1.58, 0, "O"], [3.83, -0.79, 0, "H"], [3.83, -2.37, 0, "H"] ]

		self.xCoord = [[0,0,0,"X"]]
		self.coords = self.waterACoords + self.waterBCoords + self.waterCCoords + self.waterDCoords + self.xCoord

		#
		self.oxyIndices = [0,3,6,9]
		self.hyIndices = [ [1,2], [4,5], [7,8], [10,11] ]
		self.distFilterIndices = [12]
		self.distFilterVals = [0,3], [3,5] #AB should be one group, with CD as the other. + easy to flip this
		self.maxOO = 3 #AC h-bond would be possible iff this was set high enough i suspect
		self.acceptor = True
		self.donor = True

		self.createTestObjs()

	def createTestObjs(self):
		#Geom
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams,lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Get the test obj
		currArgs = [self.oxyIndices, self.hyIndices, self.distFilterIndices, self.distFilterVals]
		currKwargs = {"acceptor":self.acceptor, "donor":self.donor, "maxOO":self.maxOO}
		self.testObj = tCode._DiscHBondCounterBetweenGroupsWithOxyDistFilterPopulator(*currArgs, **currKwargs)

		#Get an outDict[covers one interface]
		self.outDict = dict()

		#Get a sparse matrix calculator
		self.sparseMatrixCalculator = atomComboCoreHelp._SparseMatrixCalculatorStandard([self.testObj])

	def _runTestFunct(self):
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

	def _runTestFunct_outDictInterface(self):
		levelA, levelB = 0, 1
		self.testObj.populateMatrices(self.cellA, self.outDict, levelA)
		self.testObj.populateMatrices(self.cellA, self.outDict, levelB)


	#For now I'm assuming all O-O dists calculated. I could optimise a bit and calculate only some of these based on maxDist but...
	def _loadExpDistMatrixA(self):
		expDistMatrix = np.empty( (13,13) )
		expDistMatrix[:] = np.nan

		#O-O distances
		expDistMatrix[0][0], expDistMatrix[3][3] = 0,0
		expDistMatrix[6][6], expDistMatrix[9][9] = 0,0

		expDistMatrix[0][3], expDistMatrix[3][0] = 2,2
		expDistMatrix[0][6], expDistMatrix[6][0] = 3.586753406633916, 3.586753406633916
		expDistMatrix[0][9], expDistMatrix[9][0] = 3.586753406633916, 3.586753406633916

		expDistMatrix[3][6], expDistMatrix[6][3] = 1.9961963831246665, 1.9961963831246665
		expDistMatrix[3][9], expDistMatrix[9][3] = 1.9961963831246665, 1.9961963831246665

		expDistMatrix[6][9], expDistMatrix[9][6] = 2*1.58, 2*1.58

		#O-filter atom cases
		expDistMatrix[0][12], expDistMatrix[12][0] = 0, 0
		expDistMatrix[3][12], expDistMatrix[12][3] = 2, 2
		expDistMatrix[6][12], expDistMatrix[12][6] = 3.586753406633916, 3.586753406633916
		expDistMatrix[9][12], expDistMatrix[12][9] = 3.586753406633916, 3.586753406633916

		return expDistMatrix



	def testExpected_noMatrixPresent_donorAndAcceptor(self):
		expDistMatrix = self._loadExpDistMatrixA()
		anglesMatrix = np.empty( (13,13,13) )
		anglesMatrix[:] = np.nan

		anglesMatrix[3][6][7]  = 179.99999914296953
		anglesMatrix[3][6][8]  = 75.34718105754148
		anglesMatrix[6][3][4]  = 5.416432150147496e-06 #B->C donor
		anglesMatrix[6][3][5]  = 104.65281894245852
		anglesMatrix[3][9][10] = 75.34718105754148
		anglesMatrix[3][9][11] = 179.99999914296953
		anglesMatrix[9][3][4]  = 104.65281894245852
		anglesMatrix[9][3][5]  = 5.416432150147496e-06 #B->D donor

		self._runTestFunct()
		actDistMatrix = self.sparseMatrixCalculator.outDict["distMatrix"]
		actAnglesMatrix = self.sparseMatrixCalculator.outDict["angleMatrix"]

#		import pdb
#		pdb.set_trace()

		self.assertTrue( np.allclose(expDistMatrix, actDistMatrix ,equal_nan=True) )
		self.assertTrue( np.allclose(anglesMatrix, actAnglesMatrix, equal_nan=True) )

#Angles needed by default (BC and BD should be all that are needed):
# [3,6,7] B->C acceptor possibility
# [3,6,8] B->C acceptor possibility
# [6,3,4] B->C donor possibility
# [6,3,5] B->C donor possibility
# [3,9,10] B->D acceptor possibility
# [3,9,11] B->D acceptor possibility
# [9,3,4] B->D donor possibility
# [9,3,5] B->D donor possibility

#Code used to calculate angles
#		import gen_basis_helpers.analyse_md.calc_dists as calcDistsHelp
#		angleIndicesNeeded = [ [3,6,7], [3,6,8], [6,3,4], [6,3,5], [3,9,10], [3,9,11], [9,3,4], [9,3,5] ]
#		angles = calcDistsHelp.getInterAtomicAnglesForInpGeom(self.cellA, angleIndicesNeeded)
#		import pdb
#		pdb.set_trace()


	#Should be super EZ to pass
	def testExpected_noMatrixPresent_acceptorOnly(self):
		self.donor = False
		self.createTestObjs()

		expDistMatrix = self._loadExpDistMatrixA()
		anglesMatrix = np.empty( (13,13,13) )
		anglesMatrix[:] = np.nan

		anglesMatrix[3][6][7]  = 179.99999914296953
		anglesMatrix[3][6][8]  = 75.34718105754148
		anglesMatrix[3][9][10] = 75.34718105754148
		anglesMatrix[3][9][11] = 179.99999914296953

		self._runTestFunct()
		actDistMatrix = self.sparseMatrixCalculator.outDict["distMatrix"]
		actAnglesMatrix = self.sparseMatrixCalculator.outDict["angleMatrix"]

		self.assertTrue( np.allclose(expDistMatrix, actDistMatrix ,equal_nan=True) )
		self.assertTrue( np.allclose(anglesMatrix, actAnglesMatrix, equal_nan=True) )


	def testExpected_noMatrixPresent_donorOnly(self):
		self.acceptor = False
		self.createTestObjs()

		expDistMatrix = self._loadExpDistMatrixA()
		anglesMatrix = np.empty( (13,13,13) )
		anglesMatrix[:] = np.nan

		anglesMatrix[6][3][4]  = 5.416432150147496e-06 #B->C donor
		anglesMatrix[6][3][5]  = 104.65281894245852
		anglesMatrix[9][3][4]  = 104.65281894245852
		anglesMatrix[9][3][5]  = 5.416432150147496e-06 #B->D donor

		self._runTestFunct()
		actDistMatrix = self.sparseMatrixCalculator.outDict["distMatrix"]
		actAnglesMatrix = self.sparseMatrixCalculator.outDict["angleMatrix"]

		self.assertTrue( np.allclose(expDistMatrix, actDistMatrix ,equal_nan=True) )
		self.assertTrue( np.allclose(anglesMatrix, actAnglesMatrix, equal_nan=True) )


	def testExpected_matrixPresent_donorAndAcceptor(self):
		#Partially populate matrices
		self.outDict["distMatrix"] = np.empty( (13,13) )
		self.outDict["distMatrix"][:] = np.nan

		self.outDict["distMatrix"][0][3],  self.outDict["distMatrix"][3][0] = 2, 2

		self.outDict["angleMatrix"] = np.empty( (13,13,13) )
		self.outDict["angleMatrix"][:] = np.nan
		self.outDict["angleMatrix"][6][3][5] = 104.65281894245852

		#Get expected
		expDistMatrix = self._loadExpDistMatrixA()
		expAnglesMatrix = np.empty( (13,13,13) )
		expAnglesMatrix[:] = np.nan

		expAnglesMatrix[3][6][7]  = 179.99999914296953
		expAnglesMatrix[3][6][8]  = 75.34718105754148
		expAnglesMatrix[6][3][4]  = 5.416432150147496e-06 #B->C donor
		expAnglesMatrix[6][3][5]  = 104.65281894245852
		expAnglesMatrix[3][9][10] = 75.34718105754148
		expAnglesMatrix[3][9][11] = 179.99999914296953
		expAnglesMatrix[9][3][4]  = 104.65281894245852
		expAnglesMatrix[9][3][5]  = 5.416432150147496e-06 #B->D donor

		#Run with the testObj interface
		self._runTestFunct_outDictInterface()

		#Test actual vs expected
		actAngleMatrix, actDistMatrix = self.outDict["angleMatrix"], self.outDict["distMatrix"]

		self.assertTrue( np.allclose(expDistMatrix, actDistMatrix ,equal_nan=True) )
		self.assertTrue( np.allclose(expAnglesMatrix, actAngleMatrix, equal_nan=True) )




















