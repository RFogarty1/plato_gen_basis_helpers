

import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp
import plato_pylib.shared.unit_convs as uConvHelp

import gen_basis_helpers.analyse_md.traj_core as trajHelp
import gen_basis_helpers.analyse_md.binned_res as binResHelp

import gen_basis_helpers.analyse_md.surf_norm_z.get_linear_density as tCode

class TestLinearDensity(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.coordsA = [ [2,2,3,"X"], [2,2,4,"X"], [2,2,6,"Y"], [2,2,8,"X"] ]
		self.coordsB = [ [2,2,2,"X"], [2,2,6,"X"], [2,2,7,"Y"], [3,3,6,"X"] ] 
		self.massDict = {"X":10, "Y":100}
		self.massUnitConv = 1
		self.lengthUnitConv = 1 #Ang to cm generally the sensible one i guess
		self.groupIndices = None
		self.binWidth = 5
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellB = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords, self.cellB.cartCoords = self.coordsA, self.coordsB

		self.trajStepA = trajHelp.TrajStepBase(unitCell=self.cellA)
		self.trajStepB = trajHelp.TrajStepBase(unitCell=self.cellB)
		self.fullTrajA = trajHelp.TrajectoryInMemory([self.trajStepA, self.trajStepB])

	def _runTestFunct(self):
		return tCode.getLinearDensityZ(self.fullTrajA, self.binWidth, groupIndices=self.groupIndices, massDict=self.massDict)

	def testExpectedOutputAllAtomsA(self):
		expObj = self._loadExpectedAllAtoms()
		actObj = self._runTestFunct()
		self.assertEqual(expObj, actObj)

	def testExpectedOutput_onlyX(self):
		self.groupIndices = [0,1,3]
		expObj = self._loadExpectedXAtomsOnly()
		actObj = self._runTestFunct()
		self.assertEqual(expObj, actObj)

	#Value should be the same as cubic here
	def testNonCubicCell(self):
		self.lattAngles = [60,60,90]
		self.createTestObjs()
		expObj = self._loadExpectedAllAtoms()
		actObj = self._runTestFunct()
		self.assertEqual(expObj, actObj)

	def _loadExpectedAllAtoms(self):
		expBins = [2.5, 7.5]
		volumeSlice = 10*10*self.binWidth
		expValsA = [(10+10)/volumeSlice, (100+10)/volumeSlice]
		expValsB = [10/volumeSlice     , (100+10+10)/volumeSlice]
		expVals = [(a+b)/2 for a,b in zip(expValsA, expValsB)]
		expVals = [x/uConvHelp.AVOGADRO_NUMBER for x in expVals] #Since we use molar mass in the dictionary
		binVals = {"lin_den":expVals}
		expObj = binResHelp.BinnedResultsStandard.fromConstantBinWidth(expBins, self.binWidth, binVals)
		return expObj

	def _loadExpectedXAtomsOnly(self):
		expBins = [2.5, 7.5]
		volumeSlice = 10*10*self.binWidth
		expValsA = [(10+10)/volumeSlice, (10)/volumeSlice]
		expValsB = [10/volumeSlice     , (10+10)/volumeSlice]
		expVals = [(a+b)/2 for a,b in zip(expValsA, expValsB)]
		expVals = [x/uConvHelp.AVOGADRO_NUMBER for x in expVals] #Since we use molar mass in the dictionary
		binVals = {"lin_den":expVals}
		expObj = binResHelp.BinnedResultsStandard.fromConstantBinWidth(expBins, self.binWidth, binVals)
		return expObj



