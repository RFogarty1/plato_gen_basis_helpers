
import copy
import unittest

import numpy as np

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.atom_combo_core as atomComboCoreHelp
import gen_basis_helpers.analyse_md.atom_combo_opts_obj_maps as atomComboObjsMapHelp
import gen_basis_helpers.analyse_md.binned_res as binResHelp
import gen_basis_helpers.analyse_md.calc_distrib_core as calcDistribCoreHelp
import gen_basis_helpers.analyse_md.calc_radial_distrib_impl as calcRadImpl
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







