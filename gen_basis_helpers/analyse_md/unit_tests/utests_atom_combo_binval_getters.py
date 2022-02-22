
import copy
import itertools as it
import math
import unittest

import numpy as np

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.atom_combo_core as atomComboCoreHelp
import gen_basis_helpers.analyse_md.atom_combo_opts_obj_maps as atomComboObjsMapHelp
import gen_basis_helpers.analyse_md.atom_combo_populators as atomComboPopulators
import gen_basis_helpers.analyse_md.binned_res as binResHelp
import gen_basis_helpers.analyse_md.calc_distrib_core as calcDistribCoreHelp
import gen_basis_helpers.analyse_md.calc_radial_distrib_impl as calcRadImpl
import gen_basis_helpers.analyse_md.distr_opt_objs as distrOptObjHelp
import gen_basis_helpers.shared.plane_equations as planeEqnHelp

import gen_basis_helpers.analyse_md.atom_combo_binval_getters as tCode



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

		#Create opts for both planar and min dists
		dudBinsObj = None
		self.planarOpts = calcRadImpl.CalcPlanarRdfOptions(dudBinsObj, self.planarIndices, planeEqn=self.planeEqn)
		self.optsObj = calcDistribCoreHelp.CalcRdfOptions(dudBinsObj, self.fromIndices, self.toIndices, minDistAToB=True)


		#Create sparse matrix object + populate arrays
		self.sparseMatrixObj = atomComboObjsMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.planarOpts, self.optsObj])
		self.sparseMatrixObj.calcMatricesForGeom(self.cellA)

		#Create the individual binner objects
		self.planarBinner = tCode._PlanarDistsGetOneDimValsToBin(self.planeEqn, self.planarIndices)
		self.minDistBinner = tCode._MinDistsGetOneDimValsToBin(self.fromIndices, self.toIndices)

		#Create the binner object
		self.testObj = atomComboCoreHelp._GetMultiDimValsToBinFromSparseMatrices([self.planarBinner, self.minDistBinner])

	def testExpectedCaseA(self):
		#[planarA,radA], [planarB,radB]
		expVals = [ [3,2], [1,1] ]
		actVals = self.testObj.getValsToBin(self.sparseMatrixObj)
		self.assertTrue( np.allclose( np.array(expVals), np.array(actVals) ) )


class TestHozMinDistsGetter(unittest.TestCase):

	def setUp(self):
		#Geometry options
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.coords = [ [2,0,2,"A"],
		                [4,0,3,"B"],
		                [5,0,3,"C"],
		                [6,0,4,"D"] ]

		#Distribution options
		self.binResObj = binResHelp.BinnedResultsStandard.fromBinEdges([-0.1,10,20,30]) #Should be unneeded but...
		self.indicesFrom = [0,2]
		self.indicesTo = [1,3]
		self.minDistAToB = True
		self.minDistVal = 0.01

		self.createTestObjs()

	def createTestObjs(self):
		#Geom
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Create opts object
		currArgs = [self.binResObj, self.indicesFrom, self.indicesTo]
		self.optsObj = distrOptObjHelp.CalcHozDistOptions(*currArgs, minDistAToB=self.minDistAToB, minDistVal=self.minDistVal)

		#Get matrix calculator + populate
		self.sparseMatrixObj = atomComboObjsMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.optsObj])
		self.sparseMatrixObj.calcMatricesForGeom(self.cellA)

		#get binval getter
		self.testObj = atomComboObjsMapHelp.getMultiDimBinValGetterFromOptsObjs([self.optsObj])

	def _runTestFunct(self):
		return self.testObj.getValsToBin(self.sparseMatrixObj)

	def testExpectedCaseA(self):
		expVals = [ (2,), (1,) ] #2,1
		actVals = self._runTestFunct()

		for expIter, actIter in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expIter,actIter)]

	def testExpectedNoMinDistBelowThreshold(self):
		""" Common case where our group is exactly atom and we want a min-dist without including self. That needs to return zero (which can be binned if needed) instead of np.nan (which will rightfully throw an error)"""
		self.indicesFrom = [0]
		self.indicesTo = [0]
		self.createTestObjs()
		expVals = [(0,)]
		actVals = self._runTestFunct()

		for expIter, actIter in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expIter,actIter)]


class TestRadialDistribWithPlanarDists(unittest.TestCase):
	""" Testing radial distrib binner alone misses too much, hence want a planar dist val too """ 

	def setUp(self):
		#Create geometry
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.coords =  [ [0,0,7,"A"],
		                 [0,0,8,"B"],
		                 [0,0,9,"C"],
		                 [0,0,1,"D"],
		                 [0,0,2,"E"] ]

		#Simple options
		self.indicesA = [0,1]
		self.indicesB = [2,3,4]
		self.minDistAToB = False
		self.dudBinResObj = None
		self.rdfBinObj = binResHelp.BinnedResultsStandard.fromBinEdges([-0.1,10,20,30]) #Needed for filter functions
		self.rdfFilterBasedOnBins = True

		self.planeEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,0)

		self.createTestObjs()

	def createTestObjs(self):
		#Geom
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Options objects
		currArgs = [self.rdfBinObj, self.indicesA, self.indicesB]
		self.rdfOptObj = distrOptObjHelp.CalcRdfOptions(*currArgs, minDistAToB=self.minDistAToB, filterBasedOnBins=self.rdfFilterBasedOnBins)
		self.planarOptObj = calcRadImpl.CalcPlanarRdfOptions(self.dudBinResObj, self.indicesA, planeEqn=self.planeEqn)

		#Create sparse matrix populator and populate it
		self.sparseMatrixObj = atomComboObjsMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.rdfOptObj, self.planarOptObj])
		self.sparseMatrixObj.calcMatricesForGeom(self.cellA)

		#Create the binner object
		self.testObj = atomComboObjsMapHelp.getMultiDimBinValGetterFromOptsObjs([self.rdfOptObj,self.planarOptObj])

	def testExpectedCaseA(self):
		distAC, distAD, distAE = 2, 4, 5
		distBC, distBD, distBE = 1, 3, 4
		planarDistA, planarDistB = 3, 2
		expBinVals = [ [distAC,planarDistA], [distAD,planarDistA], [distAE,planarDistA],
		               [distBC,planarDistB], [distBD,planarDistB], [distBE,planarDistB] ]
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixObj)

		self.assertTrue(np.allclose(np.array(expBinVals), np.array(actBinVals)))

	def testExpectedCaseB_rdfFilterApplied(self):
		self.rdfBinObj = binResHelp.BinnedResultsStandard.fromBinEdges([-0.1,4.1])
		self.createTestObjs()

		distAC, distAD, distAE = 2, 4, 5 #distAE should be filtered out
		distBC, distBD, distBE = 1, 3, 4
		planarDistA, planarDistB = 3, 2

		expBinVals = [ [distAC,planarDistA], [distAD,planarDistA],
		               [distBC,planarDistB], [distBD,planarDistB], [distBE,planarDistB] ]

		actBinVals = self.testObj.getValsToBin(self.sparseMatrixObj)

		self.assertTrue(np.allclose(np.array(expBinVals), np.array(actBinVals)))


#Create equality tests for a couple so we can later test functions to generate from opts objs (for a couple)

class TestPlanarDistsEquality(unittest.TestCase):

	def setUp(self):
		self.planeEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,5)
		self.indices = [3,4,2]
		self.createTestObjs()

	def createTestObjs(self):
		self.testObj = tCode._PlanarDistsGetOneDimValsToBin(self.planeEqn, self.indices)

	def testCmpEqual(self):
		objA = copy.deepcopy(self.testObj)
		self.createTestObjs()
		objB = self.testObj
		self.assertEqual(objA, objB)

	def testCmpUnequal_diffLenIndices(self):
		objA = copy.deepcopy(self.testObj)
		self.indices.append(4)
		self.createTestObjs()
		objB = self.testObj
		self.assertNotEqual(objA, objB)

	def testCmpUnequal_diffPlaneEqns(self):
		objA = copy.deepcopy(self.testObj)
		self.planeEqn = planeEqnHelp.ThreeDimPlaneEquation(1,2,1,3)
		self.createTestObjs()
		objB = self.testObj
		self.assertNotEqual(objA, objB)

	def testCmpUnequal_diffObjTypes(self):
		objA = self.testObj
		objB = 4
		self.assertNotEqual(objA, objB)
		self.assertNotEqual(objB, objA)

class TestMinDistsEquality(unittest.TestCase):

	def setUp(self):
		self.fromIndices = [4,1,5]
		self.toIndices = [8,1]
		self.createTestObjs()

	def createTestObjs(self):
		self.testObj = tCode._MinDistsGetOneDimValsToBin(self.fromIndices, self.toIndices)

	def testCmpEqual(self):
		objA = copy.deepcopy(self.testObj)
		self.createTestObjs()
		objB = self.testObj
		self.assertEqual(objA, objB)

	def testCmpUnequal_diffToIndices(self):
		objA = copy.deepcopy(self.testObj)
		self.toIndices[-1] += 2
		self.createTestObjs()
		objB = self.testObj
		self.assertNotEqual(objA, objB)

	def testCmpUnequal_diffLenFromIndices(self):
		objA = copy.deepcopy(self.testObj)
		self.fromIndices.append(7)
		self.createTestObjs()
		objB = self.testObj
		self.assertNotEqual(objA, objB)

	def testCmpUnequal_diffTypes(self):
		objA = self.testObj
		objB = 6
		self.assertNotEqual(objA, objB)


class TestWaterOrientationBinValsGetter(unittest.TestCase):

	def setUp(self):
		#1) Geometric parameters
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		water_azi90   = [ [0,0,0,"O"], [-1,1,0,"H"], [1, 1, 0,"H"] ]
		water_azi120  = [ [0,0,0,"O"], [-1.37, 0.37, 0,"H"], [0.37,1.37,0,"H"] ] #Actually -60.11? So maybe 119.89 in MD land [119.88652694042403]
		water_roll70  = [ [0,0,0,"O"], [1,0.34,0.94,"H"] , [1,-0.34,-0.94,"H"] ] #Actually 70.11483488614456
		water_pitch_80_azi_20 = [ [0.0, 0.0, 0.0, 'O'], [-0.17, 0.78, 0.60, 'H'], [0.37, -0.71, 0.60, 'H'] ]#79.98644481907608,19.686775804346176


		self.cartCoords = water_azi90 + water_azi120 + water_roll70 + water_pitch_80_azi_20

		#Options object parameters
		self.binResObj = None #Irrelevant
		self.oxyIndices = [0,3,6,9]
		self.hyIndices = [ [1,2], [4,5], [7,8], [10,11] ]

		self.createTestObjs()

	def createTestObjs(self):
		#Sort out geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

		#Sort out options objects
		currArgs = [self.binResObj, self.oxyIndices, self.hyIndices]
		currKwargs = {"checkEdges":False}
		self.optsObjRoll = distrOptObjHelp.WaterOrientationOptions(*currArgs, **currKwargs, angleType="roll")
		self.optsObjPitch = distrOptObjHelp.WaterOrientationOptions(*currArgs, **currKwargs, angleType="pitch")
		self.optsObjAzi = distrOptObjHelp.WaterOrientationOptions(*currArgs, **currKwargs, angleType="azimuth")

		#Get the sparse matrix calculator and populate it
		self.sparseCalculator = atomComboObjsMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.optsObjRoll, self.optsObjPitch, self.optsObjAzi])
		self.sparseCalculator.calcMatricesForGeom(self.cellA)

		#Get the bin val getter object
		self.testObj = atomComboObjsMapHelp.getMultiDimBinValGetterFromOptsObjs([self.optsObjRoll, self.optsObjPitch, self.optsObjAzi])

	def testExpectedCaseA(self):
		#Roll, pitch, azimuth; Note that roll is maybe reversed in sign compared to the adsorbate code. Though sign is always going to be arbitrary
		#Adsorbate code is a bit weird w.r.t oh-distances and roll i think...
		expVals = [ [0                  , 0                ,90],
		            [0                  , 0                ,119.88652694042403],
		            [-70.11483488614456 , 0                ,0],
		            [-0.23900590074245542, 79.98644481907608,19.686775804346176] ]
		actVals = self.testObj.getValsToBin(self.sparseCalculator)

		self.assertTrue( np.allclose(np.array(expVals), np.array(actVals)) )

	def testExpectedCase_absRoll(self):
		#Swap from roll to the abs roll
		self.optsObjRoll.angleType = "abs_roll"
		self.testObj = atomComboObjsMapHelp.getMultiDimBinValGetterFromOptsObjs([self.optsObjRoll, self.optsObjPitch, self.optsObjAzi])

		#
		expVals = [ [0                  , 0                ,90],
		            [0                  , 0                ,119.88652694042403],
		            [ 70.11483488614456 , 0                ,0],
		            [ 0.23900590074245542, 79.98644481907608,19.686775804346176] ]
		actVals = self.testObj.getValsToBin(self.sparseCalculator)

		self.assertTrue( np.allclose(np.array(expVals), np.array(actVals)) )

	def testExpected_diffIndicesEach(self):
		""" This tests the populators a bit better; since they have to actually deal with partially-populated matrices"""
		#Create options objects
		rollOxyIndices, rollHyIndices = [0,3], [ [1,2], [4,5] ]
		pitchOxyIndices, pitchHyIndices = [3,6,9], [ [4,5], [7,8], [10,11] ]

		currArgs = [self.binResObj, rollOxyIndices, rollHyIndices]
		currKwargs = {"checkEdges":False}
		optsObjRoll = distrOptObjHelp.WaterOrientationOptions(*currArgs, **currKwargs, angleType="roll")

		currArgs = [self.binResObj, pitchOxyIndices, pitchHyIndices]
		optsObjPitch = distrOptObjHelp.WaterOrientationOptions(*currArgs, **currKwargs, angleType="pitch")

		#Create sparse calculator + populate
		sparseCalculator = atomComboObjsMapHelp.getSparseMatrixCalculatorFromOptsObjIter([optsObjRoll, optsObjPitch])
		sparseCalculator.calcMatricesForGeom(self.cellA)

		#Create bin val getter
		testObjRoll = atomComboObjsMapHelp.getOneDimBinValGetterFromOptsObj(optsObjRoll)
		testObjPitch = atomComboObjsMapHelp.getOneDimBinValGetterFromOptsObj(optsObjPitch)

		#Figure out expected values
		expValsRoll = [ 0, 0]
		expValsPitch = [ 0, 0, 79.98644481907608 ]

		#run funct + compare expected and actual values
		actValsRoll = testObjRoll.getValsToBin( sparseCalculator )
		actValsPitch = testObjPitch.getValsToBin( sparseCalculator )

		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expValsRoll, actValsRoll)]
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expValsPitch, actValsPitch)]


class TestWaterMinDistsPlusMinDistFilterBinValsGetter(unittest.TestCase):

	def setUp(self):
		#1) All geometric parameters for testing
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		#pitch=90, OH len~1, HOH angle 104.5. Then just added translation vectors (also rounded coords)
		#Then collapsed along x/y so only the z-values matter 
		self.waterACoords = [ [0.0, 0.0, 0.0, 'O'], [0, 0, 0.59, 'H'], [0, 0, 0.63, 'H'] ]
		self.waterBCoords = [ [0.0, 0.0, 5.0, 'O'], [0, 0, 5.59, 'H'], [0, 0, 5.63, 'H']]
		self.toIdxCoords = [ [0,0,3,"X"], [0,0,4,"Y"], [0,0,5,"Z"] ]
		self.coords = self.waterACoords + self.waterBCoords + self.toIdxCoords

		#2) Defining the option object args
		self.binResObj = None #Shouldnt matter for getting values to bin
		self.oxyIndices = [0,3]
		self.hyIndices = [ [1,2], [4,5] ]
		self.toIndices = [6,7]
		self.primaryIdxType = "O"

		self.filterToIndices = [8]
		self.filterDistWindow = [0,20] #Basically will encapsulate all of them here
		self.minDistType = "all"

		self.createTestObjs()

	def createTestObjs(self):
		#Geom
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Options object
		currArgs = [self.binResObj, self.oxyIndices, self.hyIndices, self.toIndices, self.filterToIndices,
		            self.filterDistWindow]
		currKwargs = {"primaryIdxType":self.primaryIdxType, "minDistType":self.minDistType}
		self.optsObj = distrOptObjHelp.WaterMinDistPlusMinDistFilterOptions(*currArgs, **currKwargs)

		#Create sparse matrix calculator + populate it
		self.sparseCalculator = atomComboObjsMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.optsObj])
		self.sparseCalculator.calcMatricesForGeom(self.cellA)

		#Create the test obj
		self.testObj = atomComboObjsMapHelp.getOneDimBinValGetterFromOptsObj(self.optsObj)

	def _runTestFunct(self):
		return self.testObj.getValsToBin(self.sparseCalculator)

	def testExpectedCase_largeFilterDistWindow(self):
		expVals = [2.37, 1.0]
		actVals = self._runTestFunct()
		[self.assertAlmostEqual(exp,act, places=6) for exp,act in it.zip_longest(expVals,actVals)]

	def testExpectedCase_smallFilterWindow(self):
		self.filterDistWindow = [0,1.5] #Should mean only "Y" is in play
		self.createTestObjs()
		expVals = [3.37, 1.0]
		actVals = self._runTestFunct()
		[self.assertAlmostEqual(exp,act, places=6) for exp,act in it.zip_longest(expVals, actVals)]

	def testRaisesWhenToIndicesAllFilteredOut(self):
		""" If toIndices() gets filtered to an empty list i want code to crash for now; this situation wont happen in current planned use case + I dont want to figure out how to best handle it yet"""
		self.filterDistWindow = [0,0.1]
		self.createTestObjs()
		with self.assertRaises(NotImplementedError):
			self._runTestFunct()

#This is also sort of a test of the populator + opt objs getter; which i havnt bothered testing separately
class TestWaterMinDistsGetter(unittest.TestCase):

	def setUp(self):
		#1) All geometric parameters for testing
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		#pitch=90, OH len~1, HOH angle 104.5. Then just added translation vectors (also rounded coords)
		#Then collapsed along x/y so only the z-values matter 
		self.waterACoords = [ [0.0, 0.0, 0.0, 'O'], [0, 0, 0.59, 'H'], [0, 0, 0.63, 'H'] ]
		self.waterBCoords = [ [0.0, 0.0, 5.0, 'O'], [0, 0, 5.59, 'H'], [0, 0, 5.63, 'H']]
		self.toIdxCoords = [ [0,0,3,"X"], [0,0,4,"Y"] ]
		self.coords = self.waterACoords + self.waterBCoords + self.toIdxCoords

		#2) For the options object
		self.binResObj = None #Shouldnt matter for getting values to bin
		self.oxyIndices = [0,3]
		self.hyIndices = [ [1,2], [4,5] ]
		self.toIndices = [6]
		self.minDistType = "all"

		self.createTestObjs()

	def createTestObjs(self):
		#Geom
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Create an options object
		currArgs = [self.binResObj, self.oxyIndices, self.hyIndices, self.toIndices]
		currKwargs = {"primaryIdxType":"O", "minDistType":self.minDistType}
		self.optsObj = distrOptObjHelp.WaterMinDistOptions(*currArgs, **currKwargs)

		#Sparse matrix calculator + populate it ()
		self.sparseCalculator = atomComboObjsMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.optsObj])
		self.sparseCalculator.calcMatricesForGeom(self.cellA)

		#Create the test object
		self.testObj = atomComboObjsMapHelp.getOneDimBinValGetterFromOptsObj(self.optsObj)

	def _runTestFunct(self):
		return self.testObj.getValsToBin(self.sparseCalculator)

	def testExpected_all(self):
		expVals = [2.37,2]
		actVals = self._runTestFunct()
		[self.assertAlmostEqual(exp,act,places=6) for exp,act in it.zip_longest(expVals, actVals)]

	def testExpected_hydrogenOnly(self):
		self.minDistType = "H"
		self.createTestObjs()
		expVals = [2.37,2.59]
		actVals = self._runTestFunct()
		[self.assertAlmostEqual(exp,act,places=6) for exp,act in it.zip_longest(expVals, actVals)]

class TestWaterPlanarMinDistGetter(unittest.TestCase):

	def setUp(self):
		#1) All geometric parameters for testing
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		#pitch=90, OH len~1, HOH angle 104.5. Then just added translation vectors (also rounded coords)
		self.waterACoords = [ [0.0, 0.0, 0.0, 'O'], [0, 0.79, 0.59, 'H'], [0, -0.79, 0.63, 'H'] ]
		self.waterBCoords = [ [0.0, 0.0, 5.0, 'O'], [0,  0.79, 5.59, 'H'], [0, -0.79, 5.63, 'H']]
		self.coords = self.waterACoords + self.waterBCoords
		self.outDict = dict()

		#Options for the populator
		self.minDistType = "all"
		self.planeEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,4)
		self.oxyIndices = [0,3]
		self.hyIndices = [ [1,2], [4,5] ]

		self.createTestObjs()

	def createTestObjs(self):
		#Geom
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams,lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Sparse matrix calculator + populate it
		currArgs = [self.oxyIndices, self.hyIndices, self.planeEqn]
		currKwargs = {"primaryOnly":False}
		self.populator = atomComboPopulators._WaterPlanarDistPopulator(*currArgs, **currKwargs)

		self.sparseMatrixObj = atomComboCoreHelp._SparseMatrixCalculatorStandard([self.populator])
		self.sparseMatrixObj.calcMatricesForGeom(self.cellA)

		#Create the test object
		currKwargs = {"minDistType":self.minDistType}
		self.testObj = tCode._WaterPlanarMinDistBinValGetter(*currArgs, **currKwargs)

	def _runTestFunct(self):
		return self.testObj.getValsToBin(self.sparseMatrixObj)

	def testExpected_minDistAll(self):
		expVals = [3.37, 1]
		actVals = self._runTestFunct()
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expVals,actVals)]

	def testExpected_minDistOxy(self):
		self.minDistType = "O"
		self.createTestObjs()
		expVals = [4,1]
		actVals = self._runTestFunct()
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expVals,actVals)]

	def testExpected_minDistHy(self):
		self.minDistType = "H"
		self.createTestObjs()
		expVals = [3.37,1.59]
		actVals = self._runTestFunct()
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expVals,actVals)]

class TestWaterPlanarDistBinValGetter(unittest.TestCase):

	#Note: Largely duplicating the populators
	def setUp(self):
		#1) All geometric parameters for testing
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		#pitch=90, OH len~1, HOH angle 104.5. Then just added translation vectors (also rounded coords)
		self.waterACoords = [ [0.0, 0.0, 0.0, 'O'], [0, 0.79, 0.61, 'H'], [0, -0.79, 0.61, 'H'] ]
		self.waterBCoords = [ [0.0, 0.0, 5.0, 'O'], [0,  0.79, 5.61, 'H'], [0, -0.79, 5.61, 'H']]
		self.coords = self.waterACoords + self.waterBCoords
		self.outDict = dict()

		#Options for the populator
		self.primaryIdxType = "O"
		self.planeEqn = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,4)
		self.oxyIndices = [0,3]
		self.hyIndices = [ [1,2], [4,5] ]

		self.createTestObjs()

	def createTestObjs(self):
		#Geom
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams,lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Sparse matrix calculator + populate it
		currArgs = [self.oxyIndices, self.hyIndices, self.planeEqn]
		currKwargs = {"primaryIdxType":self.primaryIdxType}
		self.populator = atomComboPopulators._WaterPlanarDistPopulator(*currArgs, **currKwargs)

		self.sparseMatrixObj = atomComboCoreHelp._SparseMatrixCalculatorStandard([self.populator])
		self.sparseMatrixObj.calcMatricesForGeom(self.cellA)

		#Create the test object
		self.testObj = tCode._WaterPlanarDistBinValGetter(*currArgs, **currKwargs)

	def _runTestFunct(self):
		return self.testObj.getValsToBin(self.sparseMatrixObj)

	def testExpected_oxyIndices(self):
		expVals = [4,1]
		actVals = self._runTestFunct()
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expVals,actVals)]

	def testExpected_haIndices(self):
		self.primaryIdxType = "HA"
		self.createTestObjs()
		expVals = [3.39, 1.61]
		actVals = self._runTestFunct()
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expVals, actVals)]


#Testing classes that try to bin distances/angles for hydrogen bonds
class TestGetHBondParametersBetweenGenericGroups(unittest.TestCase):

	#Stolen from TestCountHBondsBetweenGenericGroups
	def setUp(self):
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		#Coords; 
		co2Coords =   [ [0,0,0,"O"], [1,0,0,"C"], [2,0,0,"O"] ]
		waterCoords = [ [3,1,0,"H"], [4,2,0,"O"], [5,1,0,"H"] ]
		hfCoords =    [ [6,0,0,"F"], [7,0,0,"H"] ]

		self.cartCoords = co2Coords + waterCoords + hfCoords

		#from water; to co2 AND hf
		self.fromNonHyIndices = [ [4] ]
		self.fromHyIndices = [ [3,5] ]
		self.toNonHyIndices = [ [0,1,2], [6] ] #CO2 indices, then HF indices  
		self.toHyIndices = [ [], [7] ] 

		#Options for defining a h-bond; setting to slightly unphysical values to make the test simpler
		self.maxOO = 3 #O is actually sqrt(8)~2.8 away from the CO2 oxygen and HF Fluorine
		self.maxAngle = 50
		self.acceptor, self.donor = True, True

		#Misc
		self.binResObj = None
		self.primaryIndices = None
		self.optsObjInitializer = distrOptObjHelp.GetOODistsForHBondsBetweenGenericGroups

		self.createTestObjs()

	def createTestObjs(self):
		#Create the geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

		#Options object
		currArgs = [self.binResObj, self.fromNonHyIndices, self.fromHyIndices, self.toNonHyIndices, self.toHyIndices]
		currKwargs = {"acceptor":self.acceptor, "donor":self.donor, "maxOO":self.maxOO, "maxAngle":self.maxAngle, "primaryIndices":self.primaryIndices}
		self.optsObj = self.optsObjInitializer(*currArgs,**currKwargs)

		#Get a sparse matrix populator + populate it
		self.sparseCalculator = atomComboObjsMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.optsObj])
		self.sparseCalculator.calcMatricesForGeom(self.cellA)

		#Create the test object
		self.testObj = atomComboObjsMapHelp.getMultiDimBinValGetterFromOptsObjs([self.optsObj])

	def _runTestFunct(self):
		return self.testObj.getValsToBin(self.sparseCalculator)

	def testExpectedOOLengths(self):
		expVals = [(math.sqrt(8),), (math.sqrt(8),)]
		actVals = self._runTestFunct()
		for expIter,actIter in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expIter,actIter)]

	def testExpectedAngles(self):
		self.optsObjInitializer = distrOptObjHelp.GetOOHAnglesForHBondsBetweenGenericGroups
		self.createTestObjs()
		expVals = [ (0,), (0,) ]
		actVals = self._runTestFunct()
		for expIter,actIter in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expIter,actIter)]


class TestCountHBondsBetweenGenericGroups(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		#Coords: Do CO_2, water and HF combined geometry
		co2Coords = [  [0,0,0,"O"], [1,0,0,"C"], [2,0,0,"O"] ]
		waterCoords = [ [3,1,0,"H"],[4,2,0,"O"],[5,1,0,"H"] ] #H-O-H probably 90 degrees here
		hfCoords = [ [6,0,0,"F"], [7,0,0,"H"] ]

		self.cartCoords = co2Coords + waterCoords + hfCoords

		#From water; to co2 AND hf
		self.fromNonHyIndices = [ [4] ]
		self.fromHyIndices = [ [3,5] ]
		self.toNonHyIndices = [ [0,1,2], [6] ] #CO2 indices, then HF indices
		self.toHyIndices = [ [] ,[7] ]

		#Options for defining a h-bond; setting to slightly unphysical values to make the test simpler
		self.maxOO = 3 #O is actually sqrt(8)~2.8 away from the CO2 oxygen and HF Fluorine
		self.maxAngle = 50
		self.acceptor, self.donor = True, True

		#Misc
		self.binResObj = None
		self.primaryIndices = None

		self.createTestObjs()

	def createTestObjs(self):
		#Create the geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams,lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

		#Options object
		currArgs = [self.binResObj, self.fromNonHyIndices, self.fromHyIndices, self.toNonHyIndices, self.toHyIndices]
		currKwargs = {"acceptor":self.acceptor, "donor":self.donor, "maxOO":self.maxOO, "maxAngle":self.maxAngle, "primaryIndices":self.primaryIndices}
		self.optsObj = distrOptObjHelp.CountHBondsBetweenGenericGroupsOptions(*currArgs, **currKwargs)

#		#Get a sparse matrix populator + populate it
		self.sparseCalculator = atomComboObjsMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.optsObj])
		self.sparseCalculator.calcMatricesForGeom(self.cellA)

		#Create the test object
		self.testObj = atomComboObjsMapHelp.getMultiDimBinValGetterFromOptsObjs([self.optsObj])

	def _runTestFunct(self):
		return self.testObj.getValsToBin(self.sparseCalculator)

	def testExpectedTotalHBonds_fromWater(self):
		expVals = [(2,)]
		actVals = self._runTestFunct()
		for expIter,actIter in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expIter,actIter)]

	def testExpectedDonorHBonds_fromWater(self):
		self.acceptor = False
		self.createTestObjs()

		expVals = [(2,)]
		actVals = self._runTestFunct()
		for expIter,actIter in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expIter,actIter)]

	def testExpectedAcceptorHBonds_fromWater(self):
		self.donor = False
		self.createTestObjs()

		expVals = [(0,)]
		actVals = self._runTestFunct()
		for expIter,actIter in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expIter,actIter)]


	def testExpectedTotalHBonds_fromHF(self):
		#0) Setup
		self.fromNonHyIndices = [ [6] ]
		self.fromHyIndices = [ [7] ]
		self.toNonHyIndices = [ [0,1,2], [4] ] #CO2/water
		self.toHyIndices = [ [], [3,5] ]
		self.createTestObjs()

		#1) Run + compare
		expVals = [(1,)] #It only accepts one h-bond
		actVals = self._runTestFunct()
		for expIter,actIter in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expIter,actIter)]

	def testRaises_twoNonHyPlusHyGroups(self):
		""" Not suitable for this case, since code doesnt know which H is attached to which nonHyIdx """
		self.toHyIndices[0] = [4]
		with self.assertRaises(ValueError):
			self.createTestObjs()


class TestCountHBondsBetweenSpecifiedWaterGroups(unittest.TestCase):

	#Take from the TestDiscHBondCounterBetweenGroupsOxyDistFilter unit tests
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
		self.fromOxyIndices = [0,3]
		self.fromHyIndices = [ [1,2], [4,5] ]

		self.toOxyIndices = [6,9]
		self.toHyIndices  = [ [7,8], [10,11] ]

		self.maxOO = 3 #AC h-bond would be possible iff this was set high enough i suspect
		self.maxAngle = 35
		self.acceptor = True
		self.donor = True
		self.binResObj = None

		self.createTestObjs()

	def createTestObjs(self):
		#Create the geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams,lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Options object
		currArgs = [self.binResObj, self.fromOxyIndices, self.fromHyIndices, self.toOxyIndices, self.toHyIndices]
		currKwargs = {"acceptor":self.acceptor,"donor":self.donor, "maxOO":self.maxOO, "maxAngle":35}
		self.optsObj = distrOptObjHelp.CountHBondsBetweenWaterGroupsOptions(*currArgs, **currKwargs)

		#Get a sparse matrix populator + populate it
		self.sparseCalculator = atomComboObjsMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.optsObj])
		self.sparseCalculator.calcMatricesForGeom(self.cellA)

		#Create the test object
		self.testObj = atomComboObjsMapHelp.getMultiDimBinValGetterFromOptsObjs([self.optsObj])

	def _runTestFunct(self):
		return self.testObj.getValsToBin(self.sparseCalculator)

	def testExpectedTotalHBonds(self):
		expVals = [ (0,), (2,)]
		actVals = self._runTestFunct()
		for expIter,actIter in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expIter,actIter)]

	def testExpectedWithDonorAndAcceptorCounter(self):
		#Remake the objects
		currArgs = [self.binResObj, self.fromOxyIndices, self.fromHyIndices, self.toOxyIndices, self.toHyIndices]
		currKwargs = {"acceptor":False, "donor":True, "maxOO":self.maxOO, "maxAngle":35}
		self.donorOptsObj = distrOptObjHelp.CountHBondsBetweenWaterGroupsOptions(*currArgs, **currKwargs)

		currKwargs = {"acceptor":True, "donor":False, "maxOO":self.maxOO, "maxAngle":35}
		self.acceptorOptsObj = distrOptObjHelp.CountHBondsBetweenWaterGroupsOptions(*currArgs, **currKwargs)

		optsObjs = [self.acceptorOptsObj,self.donorOptsObj]

		#Get a sparse matrix populator + populate it
		self.sparseCalculator = atomComboObjsMapHelp.getSparseMatrixCalculatorFromOptsObjIter(optsObjs)
		self.sparseCalculator.calcMatricesForGeom(self.cellA)

		#Create the test object
		self.testObj = atomComboObjsMapHelp.getMultiDimBinValGetterFromOptsObjs(optsObjs)

		#Figure out expected + test		
		expVals = [ (0,0), (0,2) ]
		actVals = self._runTestFunct()

		for expIter,actIter in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expIter,actIter)]

	def testExpectedWithinGroup(self):
		self.toOxyIndices = [0,3]
		self.toHyIndices = [ [1,2], [4,5] ]
		self.acceptor = False
		self.createTestObjs()

		#Figure out expected + test		
		expVals = [ (1,), (0,) ]
		actVals = self._runTestFunct()

		for expIter,actIter in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expIter,actIter)]
	

class TestDiscHBondCounterBetweenGroupsOxyDistFilter(unittest.TestCase):

	#NOTE: Taken from the populators mostly
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
		self.maxAngle = 35
		self.acceptor = True
		self.donor = True

		self.createTestObjs()

	def createTestObjs(self):
		#Create the geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams,lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Get a sparse matrix populator + populate it
		currArgs = [self.oxyIndices, self.hyIndices, self.distFilterIndices, self.distFilterVals, self.acceptor, self.donor]
		currKwargs = {"maxOO":self.maxOO}		
		self.populator = atomComboPopulators._DiscHBondCounterBetweenGroupsWithOxyDistFilterPopulator(*currArgs, **currKwargs)

		self.sparseMatrixObj = atomComboCoreHelp._SparseMatrixCalculatorStandard([self.populator])
		self.sparseMatrixObj.calcMatricesForGeom(self.cellA)

		#Get the test obj
		currKwargs = {"maxOO":self.maxOO,"maxAngle":self.maxAngle}		
		self.testObj = tCode._DiscHBondCounterBetweenGroupsWithOxyDistFilterOneDimValGetter(*currArgs, **currKwargs)


	def _runTestFunct(self):
		return self.testObj.getValsToBin(self.sparseMatrixObj)

	def testExpectedCaseA_donorAndAcceptor(self):
		#Note that waterB donates to waterC and waterD; and thats all for default
		expVals = [0,2,0,0]
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)

	def testExpectedCaseA_donorOnly(self):
		self.donor, self.acceptor = True, False
		self.createTestObjs()
		expVals = [0,2,0,0]
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)

	def testExpectedCaseA_acceptorOnly(self):
		self.donor, self.acceptor = False, True
		self.createTestObjs()
		expVals = [0,0,0,0]
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)

	def testExpectedCase_flippedFilterValues(self):
		self.distFilterVals = [x for x in reversed(self.distFilterVals)]
		self.createTestObjs()
		expVals = [0,0,1,1]
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)

	#1 donor and 1 acceptor each
	def testExpectedCase_betweenSame(self):
		self.distFilterVals = [ self.distFilterVals[0], self.distFilterVals[0] ]
		self.createTestObjs()
		expVals = [1,1,0,0]
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)

	def testDistFilterIndicesNone(self):
		self.distFilterIndices = None
		self.createTestObjs()
		expVals = [1,3,1,1]
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)



class TestGetDistsForDiatomOpts(unittest.TestCase):

	def setUp(self):
		#All geometric params
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.hydroxylA = [ [0,0,0,"O"], [0,0,1,"H"] ] #1 length; 0 degree angle
		self.hydroxylB = [ [3,0,0,"O"], [3,1,1,"H"] ] #sqrt(2) length; 45 degree angle
		self.hydroxylC = [ [5,0,0,"O"], [7,0,0,"H"] ] #2 length;90 degrees angle

		self.coords = self.hydroxylA + self.hydroxylB + self.hydroxylC

		#
		self.binResObj = None
		self.diatomicIndices = [ [0,1], [2,3], [4,5] ] 
		self.createTestObjs()

	def createTestObjs(self):
		#Create the geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Create an options object
		self.optsObj = distrOptObjHelp.GetDistsForDiatomOpts(self.binResObj, self.diatomicIndices)

		#Get a sparse matrix populator + populate it
		self.sparseCalculator = atomComboObjsMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.optsObj])
		self.sparseCalculator.calcMatricesForGeom(self.cellA)

		#Create the test object
		self.testObj = atomComboObjsMapHelp.getMultiDimBinValGetterFromOptsObjs([self.optsObj])

	def _runTestFunct(self):
		return self.testObj.getValsToBin(self.sparseCalculator)

	def testExpectedValsA(self):

		expVals = [ (1,), (math.sqrt(2),), (2,) ]
		actVals = self._runTestFunct()

		for exp,act in it.zip_longest(expVals, actVals):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp,act)]


class TestGetAngleWithGenericVectorForDiatomOpts(unittest.TestCase):

	def setUp(self):
		#All geometric params
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.hydroxylA = [ [0,0,0,"O"], [0,0,1,"H"] ] #1 length; 0 degree angle
		self.hydroxylB = [ [3,0,0,"O"], [3,1,1,"H"] ] #sqrt(2) length; 45 degree angle
		self.hydroxylC = [ [5,0,0,"O"], [7,0,0,"H"] ] #2 length;90 degrees angle

		self.coords = self.hydroxylA + self.hydroxylB + self.hydroxylC

		#
		self.binResObj = None
		self.diatomicIndices = [ [0,1], [2,3], [4,5] ] 
		self.leftToRight = True
		self.inpVector = [ 0,0,1 ] 
		self.createTestObjs()


	def createTestObjs(self):
		#Create the geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Create an options object
		currArgs = [self.binResObj,self.diatomicIndices, self.inpVector]
		currKwargs = {"leftToRight":self.leftToRight}
		self.optsObj = distrOptObjHelp.GetAngleWithGenericVectorForDiatomOpts(*currArgs, **currKwargs)

		#Get a sparse matrix populator + populate it
		self.sparseCalculator = atomComboObjsMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.optsObj])
		self.sparseCalculator.calcMatricesForGeom(self.cellA)

		#Create the test object
		self.testObj = atomComboObjsMapHelp.getMultiDimBinValGetterFromOptsObjs([self.optsObj])

	def _runTestFunct(self):
		return self.testObj.getValsToBin(self.sparseCalculator)

	def testExpectedValsA(self):
		expVals = [(0,), (45,), (90,)]
		actVals = self._runTestFunct()
		for exp,act in it.zip_longest(expVals, actVals):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp,act)]

	def testExpectedReversedDirection(self):
		self.leftToRight=False
		self.createTestObjs()
		expVals = [(180,), (135,), (90,)]
		actVals = self._runTestFunct()
		for exp,act in it.zip_longest(expVals, actVals):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp,act)]

	def testExpectedFromOverlappingOpts_sameVector(self):
		""" Testing we get the correct results when using two options objects with same inpVector; mainly test for the populator for filling in a partially populated matrix """
		#Create the first populator
		self.diatomicIndices = [ [0,1], [2,3] ]
		self.createTestObjs()
		matrixCalcA = copy.deepcopy(self.sparseCalculator)

		#Create the second populator
		self.diatomicIndices = [ [4,5] ]
		self.createTestObjs()
		matrixCalcB = copy.deepcopy(self.sparseCalculator)

		#Create the total populator and a binval getter that encompasses all
		self.diatomicIndices = [ [0,1], [2,3], [4,5] ]
		self.createTestObjs()
		self.sparseCalculator = atomComboCoreHelp._SparseMatrixCalculatorStandard(matrixCalcA.populators + matrixCalcB.populators) 
		self.sparseCalculator.calcMatricesForGeom(self.cellA)

		expVals = [(0,), (45,), (90,)]
		actVals = self._runTestFunct()
		for exp,act in it.zip_longest(expVals, actVals):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp,act)]




