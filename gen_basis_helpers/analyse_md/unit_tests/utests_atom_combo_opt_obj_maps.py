

import unittest

import gen_basis_helpers.analyse_md.atom_combo_populators as atomComboPopulatorHelp
import gen_basis_helpers.analyse_md.atom_combo_binval_getters as binValGettersHelp
import gen_basis_helpers.analyse_md.atom_combo_core as atomComboCoreHelp
import gen_basis_helpers.analyse_md.calc_distrib_core as calcDistribCoreHelp
import gen_basis_helpers.analyse_md.calc_radial_distrib_impl as calcDistrImplHelp

import gen_basis_helpers.shared.plane_equations as planeEqnHelp

import gen_basis_helpers.analyse_md.atom_combo_opts_obj_maps as tCode



class TestGetMatrixPopulators(unittest.TestCase):


	def setUp(self):
		#Distance options obj
		self.fromIndicesA = [2,4]
		self.toIndicesA = [20,14]

		#Plane equation distance options obj
		self.planarIndicesA = [5,6,2,10]
		self.planeEqnA = planeEqnHelp.ThreeDimPlaneEquation(2,4,1,2)

		self.createTestObjs()

	def createTestObjs(self):
		dudBinObj = None

		#Distance options
		self.distOptsA = calcDistribCoreHelp.CalcRdfOptions(dudBinObj, self.fromIndicesA, self.toIndicesA)

		#Planar options
		self.planarOptsA = calcDistrImplHelp.CalcPlanarRdfOptions(dudBinObj, self.planarIndicesA, planeEqn=self.planeEqnA)

		#Expected objects
		self.expDistPopulator = atomComboPopulatorHelp._DistMatrixPopulator(self.fromIndicesA, self.toIndicesA)
		self.expPlanarPopulator = atomComboPopulatorHelp._PlanarDistMatrixPopulator(self.planarIndicesA, self.planeEqnA)

	def testExpected_distMatrixOnly(self):
		expObj = self.expDistPopulator
		actObj = tCode.getMatrixPopulatorFromOptsObj(self.distOptsA)
		self.assertEqual(expObj, actObj)

	def testExpected_planarDistMatrixOnly(self):
		expObj = self.expPlanarPopulator
		actObj = tCode.getMatrixPopulatorFromOptsObj(self.planarOptsA)
		self.assertEqual(expObj, actObj)

	def testExpectedFullMatrixCalculator(self):
		expObj = atomComboCoreHelp._SparseMatrixCalculatorStandard([self.expDistPopulator, self.expPlanarPopulator])
		actObj = tCode.getSparseMatrixCalculatorFromOptsObjIter([self.distOptsA, self.planarOptsA])
		self.assertEqual(expObj, actObj)



class TestGetBinValGetters(unittest.TestCase):

	def setUp(self):
		#Distance opts
		self.fromIndicesA = [2,4]
		self.toIndicesA = [20,14]
		self.minDistAToB = True

		#plane equation opts
		self.planarIndicesA = [5,6,2,10]
		self.planeEqnA = planeEqnHelp.ThreeDimPlaneEquation(2,4,1,2)

		self.createTestObjs()

	def createTestObjs(self):
		dudBinObj = None

		#Distance options
		self.distOptsA = calcDistribCoreHelp.CalcRdfOptions(dudBinObj, self.fromIndicesA, self.toIndicesA, minDistAToB=self.minDistAToB)

		#Planar options
		self.planarOptsA = calcDistrImplHelp.CalcPlanarRdfOptions(dudBinObj, self.planarIndicesA, planeEqn=self.planeEqnA)

		#Expected objects
		self.expDistBinner = binValGettersHelp._MinDistsGetOneDimValsToBin(self.fromIndicesA, self.toIndicesA)
		self.expPlanarBinner = binValGettersHelp._PlanarDistsGetOneDimValsToBin(self.planeEqnA, self.planarIndicesA)

	def testExpected_minDistGetterOnly(self):
		actObj = tCode.getOneDimBinValGetterFromOptsObj(self.distOptsA)
		expObj = self.expDistBinner
		self.assertEqual(expObj, actObj)

	def testExpected_planarDistGetterOnly(self):
		actObj = tCode.getOneDimBinValGetterFromOptsObj(self.planarOptsA)
		expObj = self.expPlanarBinner
		self.assertEqual(expObj, actObj)

	def testExpectedMultiDimGetter(self):
		actObj = tCode.getMultiDimBinValGetterFromOptsObjs( [self.distOptsA, self.planarOptsA] ) 
		expObj = atomComboCoreHelp._GetMultiDimValsToBinFromSparseMatrices([self.expDistBinner,self.expPlanarBinner])
		self.assertEqual(expObj, actObj)













