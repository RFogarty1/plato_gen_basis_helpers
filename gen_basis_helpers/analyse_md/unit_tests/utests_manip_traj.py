

import copy
import itertools as it
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

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


class TestAddVelocitiesToTrajInMem(unittest.TestCase):

	def setUp(self):
		self.posConvFactor = 1
		self.timeConvFactor = 1
		self.everyN = 1

		#Global cell
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		#Other
		self.times = [2,4,8]
		self.coordsA = [ [1,1,1,"X"], [2,2,2,"Y"] ]  
		self.coordsB = [ [2,2,1,"X"], [4,1,4,"Y"] ]
		self.coordsC = [ [4,1,4,"X"], [8,8,8,"Y"] ]

		self.velKey = "velocities_from_pos"
		self.createTestObjs()

	def createTestObjs(self):
		#Create the cells
		cells = [uCellHelp.UnitCell(lattParams=[x for x in self.lattParams], lattAngles=[x for x in self.lattAngles]) for a in range(3)]
		for idx,coords in enumerate([self.coordsA, self.coordsB, self.coordsC]):
			cells[idx].cartCoords = coords

		#Create the trajectory
		self.stepA = trajHelp.TrajStepFlexible(unitCell=cells[0], step=0, time=self.times[0])
		self.stepB = trajHelp.TrajStepFlexible(unitCell=cells[1], step=1, time=self.times[1])
		self.stepC = trajHelp.TrajStepFlexible(unitCell=cells[2], step=2, time=self.times[2])
		self.inpTraj = trajHelp.TrajectoryInMemory( [self.stepA, self.stepB, self.stepC] )

	def _runTestFunct(self):
		args = [self.inpTraj]
		kwargs = {"posConvFactor":self.posConvFactor, "timeConvFactor":self.timeConvFactor,
		          "velKey":self.velKey, "everyN":self.everyN}
		return tCode.addVelocitiesToTrajInMemNVT(*args, **kwargs)

	def _loadVelocitiesCaseA(self):
		deltaTimeA, deltaTimeB = 2, 4

		#So coordsB-coordsA
		velA = [ [(2-1)/deltaTimeA, (2-1)/deltaTimeA, (1-1)/deltaTimeA],
		         [(4-2)/deltaTimeA, (1-2)/deltaTimeA, (4-2)/deltaTimeA] ]

		#CoordsC - coordsB here
		velB = [ [(4-2)/deltaTimeB, (1-2)/deltaTimeB, (4-1)/deltaTimeB],
		         [(8-4)/deltaTimeB, (8-1)/deltaTimeB, (8-4)/deltaTimeB] ]

		return [velA, velB]

	def testExpectedCaseA(self):
		expTraj = copy.deepcopy(self.inpTraj)
		expVelocities = self._loadVelocitiesCaseA()

		expTraj.trajSteps[0].addExtraAttrDict( {self.velKey: {"value":expVelocities[0], "cmpType":"numericalArray"} })
		expTraj.trajSteps[1].addExtraAttrDict( {self.velKey: {"value":expVelocities[1], "cmpType":"numericalArray"} })
		self._runTestFunct()

		self.assertEqual(expTraj, self.inpTraj)
		self.assertEqual(self.inpTraj, expTraj)

	def testExpectedEveryTwo(self):
		self.everyN = 2

		expTraj = copy.deepcopy(self.inpTraj)
		expVelocities = self._loadVelocitiesCaseA()
		expTraj.trajSteps[0].addExtraAttrDict( {self.velKey: {"value":expVelocities[0], "cmpType":"numericalArray"} } )
		self._runTestFunct()
		self.assertEqual(expTraj, self.inpTraj)
		self.assertEqual(self.inpTraj, expTraj)

	def testExpectedCaseA_plusConversionFactors(self):
		self.posConvFactor = 2
		self.timeConvFactor = 10
		self.createTestObjs()

		#Figure out expected
		expTraj = copy.deepcopy(self.inpTraj)
		expVelocities = self._loadVelocitiesCaseA()
		convFactor = self.posConvFactor/self.timeConvFactor

		for idx,expVelOneAtom in enumerate(expVelocities[0]):
			expVelocities[0][idx] = [ x*convFactor for x in expVelocities[0][idx] ]
		for idx,expVelOneAtom in enumerate(expVelocities[1]):
			expVelocities[1][idx] = [ x*convFactor for x in expVelocities[1][idx] ]

		expTraj.trajSteps[0].addExtraAttrDict( {self.velKey: {"value":expVelocities[0], "cmpType":"numericalArray"} })
		expTraj.trajSteps[1].addExtraAttrDict( {self.velKey: {"value":expVelocities[1], "cmpType":"numericalArray"} })

		#Run test and compare
		self._runTestFunct()
		self.assertEqual(expTraj, self.inpTraj)
		self.assertEqual(self.inpTraj, expTraj)



class TestAddAtomicTempsToTraj(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.velsA = [ [0,0,100], [0,0,-100], [0,200,0] ]
		self.velsB = [ [0,0,200], [0,0,-300], [0,500,0] ]

		self.coordsA = [ [0,0,0,"X"], [1,1,1,"Y"], [2,2,2,"Y"] ]
		self.coordsB =  [ [0,0,0,"X"], [1,1,1,"Y"], [2,2,2,"Y"] ]
		self.massDict = {"X":2, "Y":4}
		self.velKey = "vel_key"
		self.atomTempsKey = "atom_temps"
		self.inpIndices = [0,1,2]
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellB = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coordsA
		self.cellB.cartCoords = self.coordsB

		extraAttrDictA = {self.velKey: {"value":self.velsA, "cmpType":"numericalArray"}}
		extraAttrDictB = {self.velKey: {"value":self.velsB, "cmpType":"numericalArray"}}
		self.trajStepA = trajHelp.TrajStepFlexible(unitCell=self.cellA, extraAttrDict=extraAttrDictA)
		self.trajStepB = trajHelp.TrajStepFlexible(unitCell=self.cellB, extraAttrDict=extraAttrDictB)

		self.inpTrajA = trajHelp.TrajectoryInMemory([self.trajStepA,self.trajStepB])

	def _runTrajStepFunct(self):
		args = [self.trajStepA]
		kwargs = {"velKey":self.velKey, "massDict":self.massDict}
		return tCode._getAtomicTempArrayForTrajStep(*args, **kwargs)

	def _runModTrajFunct(self):
		args = [self.inpTrajA]
		kwargs = {"atomTempKey":self.atomTempsKey, "massDict":self.massDict, "velKey":self.velKey}
		return tCode.addAtomicTempsToTraj(*args, **kwargs)

	def testExpectedCaseForTrajStepA(self):
		#Expected velocitiy simply calculated in excel quickly
		expTemps = [0.80181570000775, 1.6036314000155, 6.414525600062]
		actTemps = self._runTrajStepFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expTemps, actTemps)]

	def testExpectedForTwoStepTraj(self):
		expTempsA = [0.80181570000775, 1.6036314000155, 6.414525600062]
		expTempsB = [3.207262800031, 14.4326826001395, 40.0907850003875]

		expTraj = copy.deepcopy(self.inpTrajA)
		expTraj.trajSteps[0].addExtraAttrDict( {self.atomTempsKey: {"value":expTempsA, "cmpType":"numericalArray"} } )
		expTraj.trajSteps[1].addExtraAttrDict( {self.atomTempsKey: {"value":expTempsB, "cmpType":"numericalArray"} } )

		self._runModTrajFunct()
		actTraj = self.inpTrajA

		self.assertEqual(expTraj, actTraj)

	def testExpectedWhenOneStepHasNoVelocities(self):
		expTempsB = [3.207262800031, 14.4326826001395, 40.0907850003875]
		self.inpTrajA.trajSteps[0].__dict__.pop(self.velKey)
		self.inpTrajA.trajSteps[0].numericalArrayCmpAttrs.pop()
		expTraj = copy.deepcopy(self.inpTrajA)
		expTraj.trajSteps[1].addExtraAttrDict( {self.atomTempsKey: {"value":expTempsB, "cmpType":"numericalArray"} } )

		self._runModTrajFunct()
		actTraj = self.inpTrajA

		self.assertEqual(expTraj, actTraj)

class TestGetTimeVsTempForTraj(unittest.TestCase):

	def setUp(self):
		self.inpIndices = None
		self.tempsA = [2.1, 4.2, 1]
		self.tempsB = [3.5, 5.3, 1]
		self.times = [10, 20]
		self.atomTempKey = "atom_temps"

		self.createTestObjs()

	def createTestObjs(self):
		self.stepA = trajHelp.TrajStepFlexible(time=self.times[0], extraAttrDict={self.atomTempKey:{"value":self.tempsA, "cmpType":"numericalArray"}})
		self.stepB = trajHelp.TrajStepFlexible(time=self.times[1],extraAttrDict={self.atomTempKey:{"value":self.tempsB, "cmpType":"numericalArray"}})
		self.trajA = trajHelp.TrajectoryInMemory( [self.stepA, self.stepB] )

	def _runTestFunct(self):
		args = [self.trajA]
		kwargs = {"inpIndices": self.inpIndices, "atomTempKey": self.atomTempKey}
		return tCode.getTimeVsTempForTraj(*args, **kwargs)

	def testExpectedNoIndicesPassed(self):
		expVals = [ [self.times[0], (2.1+4.2+1)/3], [self.times[1], (3.5+5.3+1)/3] ]
		actVals = self._runTestFunct()

		for exp,act in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp,act)]

	def testExpected_indicesPassed(self):
		self.inpIndices = [0,2]
		expVals = [ [self.times[0], (2.1+1)/2], [self.times[1], (3.5+1)/2] ]
		actVals = self._runTestFunct()

		for exp,act in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp,act)]

	def testExpectedIfOneStepHasNoAtomTemps(self):
		self.trajA.trajSteps[0].__dict__.pop(self.atomTempKey)
		expVals = [[self.times[1], (3.5+5.3+1)/3]]
		actVals = self._runTestFunct()

		for exp,act in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp,act)]


class TestGetTrajBetweenTimeSteps(unittest.TestCase):

	def setUp(self):
		self.timeTol = 1e-5
		self.minTime = 0
		self.maxTime = None
		self.timeTol = 1e-5

		self.steps = [1,2,3,4]
		self.times = [20, 10,40,30]

		self.createTestFunct()

	def createTestFunct(self):
		self.stepA = trajHelp.TrajStepFlexible(step=self.steps[0], time=self.times[0])
		self.stepB = trajHelp.TrajStepFlexible(step=self.steps[1], time=self.times[1])
		self.stepC = trajHelp.TrajStepFlexible(step=self.steps[2], time=self.times[2])
		self.stepD = trajHelp.TrajStepFlexible(step=self.steps[3], time=self.times[3])

		self.allSteps = [self.stepA, self.stepB, self.stepC, self.stepD]
		self.inpTraj = trajHelp.TrajectoryInMemory(self.allSteps)

	def _runTestFunct(self):
		currKwargs = {"minTime":self.minTime, "maxTime":self.maxTime, "timeTol":self.timeTol}
		return tCode.getTrajBetweenTimes(self.inpTraj, **currKwargs)

	def _getTrajFromStepIndices(self, stepIndices):
		outSteps = list()
		for idx in stepIndices:
			currStepIdx, currTime = self.steps[idx], self.times[idx]
			currTrajStep = trajHelp.TrajStepFlexible(step=currStepIdx, time=currTime)
			outSteps.append(currTrajStep)
		return trajHelp.TrajectoryInMemory(outSteps)

	def testExpectedVals_NothingSet(self):
		expStepIndices = [0,1,2,3]
		expTraj = self._getTrajFromStepIndices(expStepIndices)
		actTraj = self._runTestFunct()
		self.assertEqual(expTraj, actTraj)

	def testExpectedVals_minSet(self):
		self.minTime = 19
		expStepIndices = [0,2,3]
		expTraj = self._getTrajFromStepIndices(expStepIndices)
		actTraj = self._runTestFunct()
		self.assertEqual(expTraj,actTraj)

	def testExpectedVals_maxSet(self):
		self.maxTime = 31
		expStepIndices = [0,1,3]
		expTraj = self._getTrajFromStepIndices(expStepIndices)
		actTraj = self._runTestFunct()
		self.assertEqual(expTraj,actTraj)

	def testExpectedvals_minAndMaxSet(self):
		self.minTime, self.maxTime = 11,31
		expStepIndices = [0,3]
		expTraj = self._getTrajFromStepIndices(expStepIndices)
		actTraj = self._runTestFunct()
		self.assertEqual(expTraj, actTraj) 

	def testExpectedVals_timeTolRelevant(self):
		self.timeTol = 1
		self.minTime = min(self.times) - 0.9*self.timeTol
		expStepIndices = [0,1,2,3]
		expTraj = self._getTrajFromStepIndices(expStepIndices)
		actTraj = self._runTestFunct()
		self.assertEqual(expTraj, actTraj)

	def testRaisesIfMinTimeGreaterThanMaxTime(self):
		self.minTime = 10
		self.maxTime = self.minTime - 1
		with self.assertRaises(ValueError):
			self._runTestFunct()

