

import copy
import unittest
import unittest.mock as mock

import gen_basis_helpers.analyse_md.traj_core as trajHelp
import gen_basis_helpers.analyse_md.manip_traj as tCode

class TestGetTrajSampledEveryNSteps(unittest.TestCase):

	def setUp(self):
		self.steps = [0,1,2,3,4,5,6]
		self.sampleEveryN = 2
		self.inPlace = False
		self.createView = True
		self.createTestObjs()

	def createTestObjs(self):
		self.trajSteps = [trajHelp.TrajStepBase(step=step) for step in self.steps]
		self.inpTrajA = trajHelp.TrajectoryInMemory(self.trajSteps)

	def _runTestFunct(self):
		currKwargs = {"inPlace":self.inPlace, "createView":self.createView}
		return tCode.getTrajSampledEveryNSteps(self.inpTrajA, self.sampleEveryN, **currKwargs)

	def testExpValsA(self):
		expSteps = [trajHelp.TrajStepBase(step=s) for s in [0,2,4,6]]
		expTraj = trajHelp.TrajectoryInMemory(expSteps)
		actTraj = self._runTestFunct()
		self.assertEqual(expTraj,actTraj)

	def testNewTrajIndependentOfOldWhenCreateViewFalse(self):
		self.createView = False
		self.sampleEveryN = 1

		#Check their equal
		expTraj = self.inpTrajA
		actTraj = self._runTestFunct()
		self.assertEqual(expTraj, actTraj)	
	
		#Check their independent
		expTraj.trajSteps[0].step += 1
		self.assertNotEqual(expTraj, actTraj)

	def testInPlaceExpected(self):
		self.inPlace = True
		self.createView = False
		expSteps = [trajHelp.TrajStepBase(step=s) for s in [0,2,4,6]]
		expTraj = trajHelp.TrajectoryInMemory(expSteps)
		self._runTestFunct()
		actTraj = self.inpTrajA
		self.assertEqual(expTraj, actTraj)

	def testRaisesIfInPlaceAndCreateViewBothTrue(self):
		self.inPlace, self.createView = True, True
		with self.assertRaises(ValueError):
			self._runTestFunct()

class TestGetTrajSplitIntoPieces(unittest.TestCase):

	def setUp(self):
		self.steps = [0,1,2,3,4,5,6]
		self.nStepsEach = 2
		self.createView = True
		self.createTestObjs()

	def createTestObjs(self):
		self.trajSteps = [trajHelp.TrajStepBase(step=step) for step in self.steps]
		self.inpTrajA = trajHelp.TrajectoryInMemory(self.trajSteps)

	def _runTestFunct(self):
		return tCode.getTrajSplitIntoEqualSections(self.inpTrajA, self.nStepsEach, createView=self.createView)

	def _getExpectedValsA(self):
		expTrajA = trajHelp.TrajectoryInMemory( [trajHelp.TrajStepBase(step=step) for step in [0,1]] )
		expTrajB = trajHelp.TrajectoryInMemory( [trajHelp.TrajStepBase(step=step) for step in [2,3]] )
		expTrajC = trajHelp.TrajectoryInMemory( [trajHelp.TrajStepBase(step=step) for step in [4,5]] )
		return [expTrajA, expTrajB, expTrajC]

	def testExpectedSplitsA(self):
		""" Note we ignore the last step here since the only other option would be to return unequal slices """
		expOutput = self._getExpectedValsA()
		actOutput = self._runTestFunct()
		self.assertEqual(expOutput, actOutput)

	def testOutputTrajsIndependentOfInpWhenCreateViewFalse(self):
		self.createView = False
		actOutput = self._runTestFunct()
		self.assertEqual( actOutput[0].trajSteps[0], self.inpTrajA.trajSteps[0] )
		actOutput[0].trajSteps[0].step += 10
		self.assertNotEqual( actOutput[0].trajSteps[0], self.inpTrajA.trajSteps[0] )




