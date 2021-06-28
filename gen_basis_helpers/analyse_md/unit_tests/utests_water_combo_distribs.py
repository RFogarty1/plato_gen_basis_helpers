

import unittest
import unittest.mock as mock

import numpy as np

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.binned_res as binResHelp
import gen_basis_helpers.analyse_md.traj_core as trajCoreHelp
import gen_basis_helpers.analyse_md.water_rotations as waterRotHelp

import gen_basis_helpers.analyse_md.water_combo_distrs as tCode


class TestGetCombinedWaterRotationDistribs_iterOfOpts(unittest.TestCase):

	def setUp(self):
		#Some coords
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		#Bins
		self.aziEdgesA   = [-20, 10, 100]
		self.pitchEdgesA = [-30, 0, 30] 
		self.planarEdgesA = [0,2,7.1] #7.1 means centre of mass likely would exclude waterAzi90_pitch20 from bin; but using oxy means it wont

		self.aziEdgesB = [-50, 5, 100]
		self.pitchEdgesB = [-30, 10, 50]
		self.planarEdgesB = [0, 10]

		#Other options
		self.waterIndicesA = [ [0,1,2], [6,7,8], [9,10,11] ]
		self.waterIndicesB = [ [3,4,5], [6,7,8], [11,9,10] ] #Swap the order of elements in last; should break some possible implementations

		self.createTestObjs()

	def createTestObjs(self):
		#Define some rotated water molecules
		#setting up water with an angle of 90 degrees and sqrt(2) bondlengths makes this simpler
		waterAzi90_pitch20 = [ [0,0,0+7,"O"], [ -1, 0.94, 0.34+7, "H"] , [1   ,0.94, 0.34+7,"H" ] ] #Translated by +7 for planar dists stuff 
		waterAzi0_pitch20  = [ [0,0,0,"O"], [0.94, 1.0, 0.34, "H"] , [0.94, -1 , 0.34,"H"] ] 
		waterAzi0_pitchm20 = [ [0,0,0,"O"], [0.94, 1.0, -0.34, "H"], [0.94, -1 ,-0.34,"H"] ] 

		#Sort the trajectory; default is to use a one step trajectory
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = waterAzi90_pitch20 + waterAzi0_pitch20 + waterAzi0_pitchm20 + waterAzi0_pitchm20
		self.stepA = trajCoreHelp.TrajStepFlexible(unitCell=self.cellA)
		self.trajA = trajCoreHelp.TrajectoryInMemory([self.stepA])

		#Create the options objects
		self.aziBinsA = binResHelp.BinnedResultsStandard.fromBinEdges(self.aziEdgesA)
		self.pitchBinsA = binResHelp.BinnedResultsStandard.fromBinEdges(self.pitchEdgesA)
		self.planarBinsA = binResHelp.BinnedResultsStandard.fromBinEdges(self.planarEdgesA)

		self.aziBinsB = binResHelp.BinnedResultsStandard.fromBinEdges(self.aziEdgesB)
		self.pitchBinsB = binResHelp.BinnedResultsStandard.fromBinEdges(self.pitchEdgesB)
		self.planarBinsB = binResHelp.BinnedResultsStandard.fromBinEdges(self.planarEdgesB)

		#Create options objs
		self.aziOptsA = waterRotHelp.CalcStandardWaterOrientationDistribOptions(self.aziBinsA , self.waterIndicesA, angleType="azimuth")
		self.pitchOptsA = waterRotHelp.CalcStandardWaterOrientationDistribOptions(self.pitchBinsA, self.waterIndicesA, angleType="pitch")
		self.planarOptsA = tCode.CalcWaterPlanarDistribOptions_fromOxy(self.planarBinsA, self.waterIndicesA)

		self.aziOptsB = waterRotHelp.CalcStandardWaterOrientationDistribOptions(self.aziBinsB , self.waterIndicesB, angleType="azimuth")
		self.pitchOptsB = waterRotHelp.CalcStandardWaterOrientationDistribOptions(self.pitchBinsB, self.waterIndicesB, angleType="pitch")
		self.planarOptsB = tCode.CalcWaterPlanarDistribOptions_fromOxy(self.planarBinsB, self.waterIndicesB)

		self.optsObjs = [ [self.aziOptsA, self.pitchOptsA], [self.aziOptsB, self.pitchOptsB] ]

	def _runTestFunct(self):
#		return tCode.getMultipleCombinedWaterRotationDistribBinsFromOptObjs(self.trajA, self.optsObjs)
		return tCode.getMultipleWaterComboDistribBinsFromOptObjs(self.trajA, self.optsObjs)

	def _getExpectedCaseA(self):
		edgesA, edgesB = [self.aziEdgesA, self.pitchEdgesA], [self.aziEdgesB, self.pitchEdgesB]
		expBinObjA, expBinObjB = binResHelp.NDimensionalBinnedResults(edgesA), binResHelp.NDimensionalBinnedResults(edgesB)
		expBinObjA.initialiseCountsMatrix()
		expBinObjB.initialiseCountsMatrix()

		#Sort the unnormalised counts
		expBinObjA.binVals["counts"][0][0] = 2
		expBinObjA.binVals["counts"][1][1] = 1

		expBinObjB.binVals["counts"][0][1] = 1
		expBinObjB.binVals["counts"][0][0] = 2

		#Sort the normalised counts (same as unnormalised for the simple case)
		expBinObjA.initialiseCountsMatrix(countKey="normalised_counts")
		expBinObjB.initialiseCountsMatrix(countKey="normalised_counts")
		np.copyto( expBinObjA.binVals["normalised_counts"], expBinObjA.binVals["counts"] )
		np.copyto( expBinObjB.binVals["normalised_counts"], expBinObjB.binVals["counts"] )

		return [expBinObjA, expBinObjB]

	def _getExpectedCaseWithPlanarDists(self):
		edgesA, edgesB = [self.aziEdgesA, self.pitchEdgesA, self.planarEdgesA], [self.aziEdgesB, self.pitchEdgesB, self.planarEdgesB]
		expBinObjA, expBinObjB = binResHelp.NDimensionalBinnedResults(edgesA), binResHelp.NDimensionalBinnedResults(edgesB)
		expBinObjA.initialiseCountsMatrix()
		expBinObjB.initialiseCountsMatrix()

		#Sort unnormalised counts
		expBinObjA.binVals["counts"][0][0][0] = 2
		expBinObjA.binVals["counts"][1][1][1] = 1

		expBinObjB.binVals["counts"][0][1][0] = 1
		expBinObjB.binVals["counts"][0][0][0] = 2

		#Sort the normalised counts (same as unnormalised for the simple case)
		expBinObjA.initialiseCountsMatrix(countKey="normalised_counts")
		expBinObjB.initialiseCountsMatrix(countKey="normalised_counts")
		np.copyto( expBinObjA.binVals["normalised_counts"], expBinObjA.binVals["counts"] )
		np.copyto( expBinObjB.binVals["normalised_counts"], expBinObjB.binVals["counts"] )

		return [expBinObjA, expBinObjB]

	def testExpectedCaseA(self):
		expBinObjs = self._getExpectedCaseA()
		actBinObjs = self._runTestFunct()
		self.assertEqual(expBinObjs, actBinObjs)


	def testExpectedCaseWithPlanarDistsA(self):
		self.optsObjs = [ [self.aziOptsA, self.pitchOptsA, self.planarOptsA], [self.aziOptsB, self.pitchOptsB, self.planarOptsB] ]
		expBinObjs = self._getExpectedCaseWithPlanarDists()
		actBinObjs = self._runTestFunct()
		self.assertEqual(expBinObjs, actBinObjs)


class TestGetCombinedWaterRotationDistribs(unittest.TestCase):

	#Restrict to 2-dimensions for reasons of my sanity
	def setUp(self):
		#Some coords
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		#Bins
		self.aziEdges   = [-20, 10, 100]
		self.pitchEdges = [-30, 0, 30] 

		#Other options
		self.waterIndicesAll = [ [0,1,2], [3,4,5], [6,7,8], [9,10,11] ]

		self.createTestObjs()

	def createTestObjs(self):
		#Define some rotated water molecules
		#setting up water with an angle of 90 degrees and sqrt(2) bondlengths makes this simpler
		waterAzi90_pitch20 = [ [0,0,0,"O"], [ -1, 0.94, 0.34, "H"] , [1   ,0.94, 0.34,"H" ] ] #bin [1][1]
		waterAzi0_pitch20  = [ [0,0,0,"O"], [0.94, 1.0, 0.34, "H"] , [0.94, -1 , 0.34,"H"] ] #bin [0][1]
		waterAzi0_pitchm20 = [ [0,0,0,"O"], [0.94, 1.0, -0.34, "H"], [0.94, -1 ,-0.34,"H"] ] #bin [0][0]

		#Sort the trajectory; default is to use a one step trajectory
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = waterAzi90_pitch20 + waterAzi0_pitch20 + waterAzi0_pitchm20 + waterAzi0_pitchm20
		self.stepA = trajCoreHelp.TrajStepFlexible(unitCell=self.cellA)
		self.trajA = trajCoreHelp.TrajectoryInMemory([self.stepA])

		#Create the options objects
		self.aziBins = binResHelp.BinnedResultsStandard.fromBinEdges(self.aziEdges)
		self.pitchBins = binResHelp.BinnedResultsStandard.fromBinEdges(self.pitchEdges)

		self.aziOpts   = waterRotHelp.CalcStandardWaterOrientationDistribOptions(self.aziBins  , self.waterIndicesAll, angleType="azimuth")
		self.pitchOpts = waterRotHelp.CalcStandardWaterOrientationDistribOptions(self.pitchBins, self.waterIndicesAll, angleType="pitch")


	def _runTestFunct(self):
		optsObjs = [self.aziOpts, self.pitchOpts]
		return tCode.getCombinedWaterRotationDistribBinsFromOptsObjs(self.trajA, optsObjs)

	#No equality for this object yet (since its really hard to implement)
	def _getExpectedBinObjCaseA(self):
		edges = [self.aziEdges, self.pitchEdges]
		expBinObj = binResHelp.NDimensionalBinnedResults(edges)
		expBinObj.initialiseCountsMatrix()

		#Sort the unnormalised counts
		expBinObj.binVals["counts"][0][0] = 2
		expBinObj.binVals["counts"][0][1] = 1
		expBinObj.binVals["counts"][1][1] = 1

		#Sort the normalised counts (same as unnormalised for the simple case)
		expBinObj.initialiseCountsMatrix(countKey="normalised_counts")
		np.copyto( expBinObj.binVals["normalised_counts"], expBinObj.binVals["counts"] )

		return expBinObj

	def testExpectedCaseA(self):
		#Get the expected counts
		expBinObj = self._getExpectedBinObjCaseA()
		actBinObj = self._runTestFunct()
		self.assertEqual( expBinObj, actBinObj )

	def testWithNormalisedCounts(self):
		""" Test that we have normalised the counts by number of steps """
		self.trajA = trajCoreHelp.TrajectoryInMemory([self.stepA, self.stepA])
		expBinObj = self._getExpectedBinObjCaseA()
		#Note we only increase counts; NOT normalised counts
		expBinObj.binVals["counts"][0][0] *= 2
		expBinObj.binVals["counts"][0][1] *= 2
		expBinObj.binVals["counts"][1][1] *= 2
		actBinObj = self._runTestFunct()

		self.assertEqual( expBinObj, actBinObj )

	def testRaisesIfIndicesDifferentBetweenOptObjs(self):
		self.aziOpts.waterIndices = [ [1,2,3], [4,5,6] ]
		self.pitchOpts.waterIndices = [ [4,5,6] ]

		with self.assertRaises(ValueError):
			self._runTestFunct()


