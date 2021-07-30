
import copy
import itertools as it
import unittest
import unittest.mock as mock

import numpy as np

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.binned_res as binResHelp
import gen_basis_helpers.analyse_md.calc_distrib_core as calcDistribCoreHelp
import gen_basis_helpers.analyse_md.calc_radial_distrib_impl as calcRadImpl
import gen_basis_helpers.analyse_md.traj_core as trajCoreHelp

import gen_basis_helpers.shared.plane_equations as planeEqnHelp


import gen_basis_helpers.analyse_md.atom_combo_distr as tCode



class TestGetMultipleCombinedDistrBins(unittest.TestCase):

	def setUp(self):
		#Some coords
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.coords = [  [2,2,2,"Mg"],
		                 [2,4,2,"O" ],
		                 [2,2,4,"O" ],
		                 [9,2,4,"Mg"] ]

		#Setup bins
		self.distBinEdgesA = [0, 1.5, 3.1, 4.5] #Mg-O Dists are 2,2,3, and a bit >3
		self.planarBinEdgesA = [0,3,5] #Dists for Mg are obv 2 and 4 (using surface plane)

		self.distBinEdgesB = [0,4,8] #Dists are 3.61, 4.1, 3
		self.planarBinEdgesB = [0,5,8]

		self.indicesA_optsA = [0,3] 
		self.indicesB_optsA = [1,2]

		self.indicesA_optsB = [0,1,2]
		self.indicesB_optsB = [3]

		self.createTestObjs()

	def createTestObjs(self):
		#Create co-ordinates
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Create a simple one-step trajectory
		self.stepA = trajCoreHelp.TrajStepFlexible(unitCell=self.cellA)
		self.trajA = trajCoreHelp.TrajectoryInMemory([self.stepA])

		#Create bins then options objects
		self.distBinsA = binResHelp.BinnedResultsStandard.fromBinEdges(self.distBinEdgesA)
		self.planarBinsA = binResHelp.BinnedResultsStandard.fromBinEdges(self.planarBinEdgesA)

		self.distBinsB = binResHelp.BinnedResultsStandard.fromBinEdges(self.distBinEdgesB)
		self.planarBinsB = binResHelp.BinnedResultsStandard.fromBinEdges(self.planarBinEdgesB)

#binResObj, indicesA, indicesB, volume=None, minDistAToB=False):
		self.minDistsOptsA = calcDistribCoreHelp.CalcRdfOptions(self.distBinsA, self.indicesA_optsA, self.indicesB_optsA, minDistAToB=True)
		self.minDistsOptsB = calcDistribCoreHelp.CalcRdfOptions(self.distBinsB, self.indicesA_optsB, self.indicesB_optsB, minDistAToB=True)

#binResObj, indices, planeEqn=None, volume=None):
		self.planarDistOptsA = calcRadImpl.CalcPlanarRdfOptions(self.planarBinsA, self.indicesA_optsA)
		self.planarDistOptsB = calcRadImpl.CalcPlanarRdfOptions(self.planarBinsB, self.indicesA_optsB)

		self.optsObjs = [ [self.minDistsOptsA, self.planarDistOptsA], [self.minDistsOptsB, self.planarDistOptsB] ]

	def _getExpectedCaseA(self):
		edgesA, edgesB = [self.distBinEdgesA, self.planarBinEdgesA], [self.distBinEdgesB, self.planarBinEdgesB]
		expBinObjA, expBinObjB = binResHelp.NDimensionalBinnedResults(edgesA), binResHelp.NDimensionalBinnedResults(edgesB)
		expBinObjA.initialiseCountsMatrix()
		expBinObjB.initialiseCountsMatrix()

		#Sort out the unnormalised counts
		expBinObjA.binVals["counts"][1][0] = 1
		expBinObjA.binVals["counts"][1][1] = 1

		expBinObjB.binVals["counts"][0][0] = 2
		expBinObjB.binVals["counts"][1][0] = 1

		#Sort the normalised counts (same as unnormalised for the simple case with 1 step)
		expBinObjA.initialiseCountsMatrix(countKey="normalised_counts")
		expBinObjB.initialiseCountsMatrix(countKey="normalised_counts")
		np.copyto( expBinObjA.binVals["normalised_counts"], expBinObjA.binVals["counts"] )
		np.copyto( expBinObjB.binVals["normalised_counts"], expBinObjB.binVals["counts"] )

		return [expBinObjA, expBinObjB]

	def _runTestFunct(self):
		return tCode.getAtomicComboDistrBinsFromOptsObjs( self.trajA, self.optsObjs )

	def testCaseA(self):
		expBins = self._getExpectedCaseA()
		actBins = self._runTestFunct()


class TestCheckIndicesConsistent(unittest.TestCase):

	def setUp(self):
		self.primaryIndicesA = [1,4,2]
		self.primaryIndicesB = [1,4,2]
		self.createTestObjs()

	def createTestObjs(self):
		dudBinObj, dudSecondaryIndices = None, [5]

		currArgs = [dudBinObj, self.primaryIndicesA, dudSecondaryIndices]
		self.objA = calcDistribCoreHelp.CalcRdfOptions(*currArgs, minDistAToB=True)
		self.objB = calcRadImpl.CalcPlanarRdfOptions(dudBinObj, self.primaryIndicesB)

	def _runTestFunct(self):
		return tCode._checkIndicesConsistentInOptsObjGroup([self.objA,self.objB])

	def testDoesntRaiseWhenConsistent(self):
		self._runTestFunct()

	def testRaisesForInconsistentA(self):
		self.primaryIndicesA.append(4)
		self.createTestObjs()
		with self.assertRaises(ValueError):
			self._runTestFunct()

	def testRaisesForInconsistentB(self):
		self.primaryIndicesA[-1] += 1
		self.createTestObjs()
		with self.assertRaises(ValueError):
			self._runTestFunct()

	def testRaisesForInconsistent_diffOrdered(self):
		self.primaryIndicesA = sorted(self.primaryIndicesA, reverse=True)
		self.primaryIndicesB = sorted(self.primaryIndicesB)
		self.createTestObjs()
		with self.assertRaises(ValueError):
			self._runTestFunct()


#Likely will remove this at some point in the future
class TestGetPrimaryIndicesFromOptsObj(unittest.TestCase):

	def setUp(self):
		self.primaryIndices = [3,4,5]
		self.createTestObjs()

	def createTestObjs(self):
		dudBinObj = None

		currArgs = [dudBinObj, self.primaryIndices, [2,8,9]]
		self.distOptsObj = calcDistribCoreHelp.CalcRdfOptions(*currArgs, minDistAToB=True)

		currArgs = [dudBinObj, self.primaryIndices]
		self.planarDistsOptsObj = calcRadImpl.CalcPlanarRdfOptions(*currArgs)


	def testExpected_radDist(self):
		expIndices = self.primaryIndices
		actIndices = tCode._getPrimaryIndicesFromOptObj(self.distOptsObj)
		self.assertEqual(expIndices, actIndices)

	def testExpected_planarDistObj(self):
		expIndices = self.primaryIndices
		actIndices = tCode._getPrimaryIndicesFromOptObj(self.planarDistsOptsObj)
		self.assertEqual(expIndices,actIndices)


class TestGetMatrixCalculator(unittest.TestCase):

	def setUp(self):
		#These should generally be the same regardless
		self.distIndicesFromA = [1,2]
		self.distIndicesFromB = [1,2]

		#
		self.distIndicesToA = [5,6,7]
		self.distIndicesToB = [5,7,12,14] 

		#Planar options stuff
		self.planarIndicesA = [2,4]
		self.planarIndicesB = [5,7,8]

		self.planarEqnA = planeEqnHelp.ThreeDimPlaneEquation(0,1,1,4)
		self.planarEqnB = None

		self.createTestObjs()

	def createTestObjs(self):
		dudBinsObj = None

		#Sort out dists objs
		currArgsA = [dudBinsObj, self.distIndicesFromA, self.distIndicesToA]
		currArgsB = [dudBinsObj, self.distIndicesFromB, self.distIndicesToB]
		self.distsOptsA = calcDistribCoreHelp.CalcRdfOptions(*currArgsA, minDistAToB=True)
		self.distsOptsB = calcDistribCoreHelp.CalcRdfOptions(*currArgsB, minDistAToB=True)

		#Sort out planar dists help
		self.planarOptsA = calcRadImpl.CalcPlanarRdfOptions(dudBinsObj, self.planarIndicesA, planeEqn=self.planarEqnA)
		self.planarOptsB = calcRadImpl.CalcPlanarRdfOptions(dudBinsObj, self.planarIndicesB, planeEqn=self.planarEqnB)

		self.useOptsObjs = [ [self.distsOptsA], [self.distsOptsB] ]

	def _runTestFunct(self):
		return tCode._getMatrixCalculatorForOptsObjsGroups( self.useOptsObjs )

	def testExpCase_rdfOnly(self):
		self.useOptsObjs = [ [self.distsOptsA], [self.distsOptsB] ]

		currKwargs = {"distIndicesFrom":[self.distIndicesFromA, self.distIndicesFromB],
		              "distIndicesTo":[self.distIndicesToA, self.distIndicesToB]}
		expCalculator = tCode._SparseMatrixCalculator(**currKwargs)
		actCalculator = self._runTestFunct()

		self.assertEqual(expCalculator, actCalculator)

	def testRaises_rdf_minDistNotSet(self):
		self.distsOptsA.minDistAToB = False
		self.useOptsObjs = [ [self.distsOptsA], [self.distsOptsB] ]
		with self.assertRaises(ValueError):
			self._runTestFunct()

	def testExpected_planeDistsOnly(self):
		""" Note this includes setting a default value for a plane equation """
		self.useOptsObjs = [ [self.planarOptsA], [self.planarOptsB] ]
		currKwargs = {"planeDistEquations":[self.planarEqnA, planeEqnHelp.ThreeDimPlaneEquation(0,0,1,0)],
		              "planeDistIndices":[self.planarIndicesA, self.planarIndicesB]}
		expCalculator = tCode._SparseMatrixCalculator(**currKwargs)
		actCalculator = self._runTestFunct()

		self.assertEqual(expCalculator, actCalculator)


class TestSparseMatrixCalculator_calcMatrices(unittest.TestCase):

	def setUp(self):
		#Some coords
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.coords = [  [2,2,2,"Mg"],
		                 [2,2,6, "O"],
		                 [2,2,4, "O"],
		                 [2,2,5,"Mg"],
		                 [2,2,8, "O"],
		                 [2,2,9,"Mg"] ]


		#Plane equations part
		self.planeEquations = [ planeEqnHelp.ThreeDimPlaneEquation(0,0,1,0),
		                        planeEqnHelp.ThreeDimPlaneEquation(0,0,1,3),
		                        planeEqnHelp.ThreeDimPlaneEquation(0,0,1,0) ]

		self.planeIndices = [  [0,1], [1,4], [0,5] ]

		#Distances part
		self.distIndicesFrom = [ [0], [0,1], [0] ]
		self.distIndicesTo   = [ [3,4], [4,5], [2] ]

		self.createTestObjs()

	def createTestObjs(self):
		#Create the geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Create the test obj
		currKwargs = {"planeDistEquations":self.planeEquations, "planeDistIndices":self.planeIndices,
		              "distIndicesFrom":self.distIndicesFrom, "distIndicesTo":self.distIndicesTo}
		self.testObj = tCode._SparseMatrixCalculator(**currKwargs)

	def _runTestFunct(self):
		self.testObj.calcMatricesForGeom(self.cellA)

	def testUniquePlaneEquations(self):
		self.planeEquations.append( planeEqnHelp.ThreeDimPlaneEquation(0,0,1,0) )
		expUniqueVals = [ self.planeEquations[0], self.planeEquations[1] ]
		actUniqueVals = self.testObj.uniquePlaneEquations
		self.assertEqual(expUniqueVals, actUniqueVals)

	def testExpectedPlaneDistArray(self):
		expArray = np.empty( (2,6) )
		expArray[:] = np.nan

		#The first UNIQUE plane equation
		expArray[0][0] = 2
		expArray[0][1] = 4 
		expArray[0][5] = 1 

		#The second UNIQUE plane equation
		expArray[1][1] = 3
		expArray[1][4] = 5

		self._runTestFunct()
		actArray = self.testObj.planarDists

		self.assertTrue( np.allclose(expArray,actArray, equal_nan=True) )

	#Note: We could optimise further; and some of these entries maybe wouldnt appear
	def testExpectedDistMatrixA(self):
		expArray = np.empty( (6,6) )
		expArray[:] = np.nan

		#
		expArray[0][2], expArray[2][0] = 2, 2
		expArray[0][3], expArray[3][0] = 3, 3
		expArray[0][4], expArray[4][0] = 4, 4
		expArray[0][5], expArray[5][0] = 3, 3

		#
		expArray[1][2], expArray[2][1] = 2, 2
		expArray[1][3], expArray[3][1] = 1, 1
		expArray[1][4], expArray[4][1] = 2, 2
		expArray[1][5], expArray[5][1] = 3, 3

		self._runTestFunct()
		actArray = self.testObj.distMatrix

		self.assertTrue( np.allclose(expArray,actArray, equal_nan=True) )


class TestSparseMatrixCalculator_simpleEquality(unittest.TestCase):

	def setUp(self):
		self.distIndicesFrom = [ [2,4,5], [1] ]
		self.distIndicesTo = [ [6,7,8], [7,8] ]
		self.planeDistEquations = [ planeEqnHelp.ThreeDimPlaneEquation(0,0,1,2) ]
		self.planeDistIndices = [ [5,6] ]

		self.createTestObjs()

	def createTestObjs(self):
		currKwargs = {"distIndicesFrom": self.distIndicesFrom, "distIndicesTo": self.distIndicesTo,
		              "planeDistEquations": self.planeDistEquations, "planeDistIndices":self.planeDistIndices}
		self.testObj = tCode._SparseMatrixCalculator(**currKwargs)

	def testEqualObjsCompareEqual(self):
		objA = copy.deepcopy(self.testObj)
		self.createTestObjs()
		objB = self.testObj
		self.assertEqual(objA, objB)

	def testUnequalObjsCompareUnequal_distIndicesNoneInOne(self):
		objA = copy.deepcopy(self.testObj)
		self.distIndicesFrom = None
		self.distIndicesTo = None
		self.createTestObjs()
		objB = self.testObj
		self.assertNotEqual(objA, objB)

	def testUnequalObjsCompareUnequal_diffLengthDistIndices(self):
		objA = copy.deepcopy(self.testObj)
		self.distIndicesFrom.append([4])
		self.createTestObjs()
		objB = self.testObj
		self.assertNotEqual(objA, objB)

	def testDiffLengthPlaneEquations(self):
		objA = copy.deepcopy(self.testObj)
		self.planeDistEquations.append( planeEqnHelp.ThreeDimPlaneEquation(1,0,0,2) ) 
		self.createTestObjs()
		objB = self.testObj
		self.assertNotEqual(objA, objB)

	def testDiffPlaneEquationVals(self):
		objA = copy.deepcopy(self.testObj)
		self.planeDistEquations[0] = planeEqnHelp.ThreeDimPlaneEquation(0,1,1,4)
		self.createTestObjs()
		objB = self.testObj
		self.assertNotEqual(objA, objB)


class TestGetMultiDimValsToBin(unittest.TestCase):

	def setUp(self):
		#Some coords
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.coords = [  [2,2,2,"Mg"],
		                 [2,2,6, "O"],
		                 [2,2,4, "O"],
		                 [2,2,5,"Mg"] ]

		#Options for planar distance
		self.planeEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,5)
		self.planarIndices = [0,1]
		self.fromIndices = [0,1]
		self.toIndices = [2,3]

		self.createTestObjs()

	def createTestObjs(self):
		#Create the geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Create sparse matrix object + populate arrays
		currKwargs = {"planeDistEquations":[self.planeEqn], "planeDistIndices":[self.planarIndices],
		              "distIndicesFrom":[self.fromIndices], "distIndicesTo":[self.toIndices]}
		self.sparseMatrixObj = tCode._SparseMatrixCalculator(**currKwargs)
		self.sparseMatrixObj.calcMatricesForGeom(self.cellA)

		#Create opts for both planar and min dists
		dudBinsObj = None
		self.planarOpts = calcRadImpl.CalcPlanarRdfOptions(dudBinsObj, self.planarIndices, planeEqn=self.planeEqn)
		self.optsObj = calcDistribCoreHelp.CalcRdfOptions(dudBinsObj, self.fromIndices, self.toIndices, minDistAToB=True)

		#Create the binner object
		self.testObj = tCode._GetMultiDimValsToBinFromSparseMatrices.fromOptsObjs([self.planarOpts, self.optsObj])

	def testExpectedCaseA(self):
		#[planarA,radA], [planarB,radB]
		expVals = [ [3,2], [1,1] ]
		actVals = self.testObj.getValsToBin(self.sparseMatrixObj)
		self.assertTrue( np.allclose( np.array(expVals), np.array(actVals) ) )


class TestGetOneDimValsToBin_planarDists(unittest.TestCase):

	def setUp(self):
		#Some coords
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.coords = [  [2,2,2,"Mg"],
		                 [2,2,6, "O"],
		                 [2,2,4, "O"],
		                 [2,2,5,"Mg"] ]

		#Options for the planar distance
		self.planeEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,2)
		self.planarIndices = [2,0,3]

		self.createTestObjs()

	def createTestObjs(self):
		#Create the geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Create an options object
		dudBinsObj = None
		self.planarOpts = calcRadImpl.CalcPlanarRdfOptions(dudBinsObj, self.planarIndices, planeEqn=self.planeEqn)

		#Create the sparse matrix object + populate arrays
		currKwargs = {"planeDistEquations":[self.planeEqn], "planeDistIndices":[self.planarIndices]}
		self.sparseMatrixObj = tCode._SparseMatrixCalculator(**currKwargs)
		self.sparseMatrixObj.calcMatricesForGeom(self.cellA)

		#Create the binner object
		self.testObj = tCode._PlanarDistsGetOneDimValsToBin.fromOptsObjs(self.planarOpts)

	def _runTestFunct(self):
		return self.testObj.getValsToBin(self.sparseMatrixObj)

	def testExpectedValsA(self):
		expVals = [2, 0, 3]
		actVals = self._runTestFunct()
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expVals,actVals)]


@unittest.skip("")
class TestGetOneDimValsToBin_minDists(unittest.TestCase):

	def setUp(self):
		#Some coords
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.coords = [  [2,2,4,"Mg"],
		                 [2,2,6, "O"],
		                 [2,2,5, "O"],
		                 [2,2,9,"Mg"] ]

		#Options for the distances to calculate
		self.fromIndices = [1,2]
		self.toIndices = [0,3]

		self.createTestObjs()

	def createTestObjs(self):
		#Create the geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Create an options object
		dudBinsObj = None
		currArgs = [dudBinsObj, self.fromIndices, self.toIndices]
		self.optsObj = calcDistribCoreHelp.CalcRdfOptions(*currArgs, minDistAToB=True)

		#Create the sparse matrix object + populate arrays
		currKwargs = {"distIndicesFrom":[self.fromIndices], "distIndicesTo":[self.toIndices]}
		self.sparseMatrixObj = tCode._SparseMatrixCalculator(**currKwargs)
		self.sparseMatrixObj.calcMatricesForGeom(self.cellA)

		#Create the binner object
		self.testObj = tCode._MinDistsGetOneDimValsToBin(self.fromIndices, self.toIndices)

	def _runTestFunct(self):
		return self.testObj.getValsToBin(self.sparseMatrixObj)
	
	def testExpectedValsA(self):
		expVals = [2, 1]
		actVals = self._runTestFunct()
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expVals,actVals)]








