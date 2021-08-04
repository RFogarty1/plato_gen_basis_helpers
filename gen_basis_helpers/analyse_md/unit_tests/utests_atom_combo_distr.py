
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
		self.assertEqual(expBins, actBins)


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



