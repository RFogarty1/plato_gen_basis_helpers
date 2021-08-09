
import copy
import itertools as it
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


