
import copy
import os
import itertools as it
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp
import gen_basis_helpers.analyse_md.traj_core as tCode



class TestTrajInMemoryReadWrite(unittest.TestCase):

	def setUp(self):
		self.unitCellA = uCellHelp.UnitCell(lattParams=[10,10,10],lattAngles=[90,90,90])
		self.timeA, self.timeB = 2, 3
		self.stepA, self.stepB = 0, 50
		self.fileNameA = "temp_file_a.traj"
		self.createTestObjs()

	def tearDown(self):
		os.remove(self.fileNameA)

	def createTestObjs(self):
		self.trajStepA = tCode.TrajStepBase(time=self.timeA, step=self.stepA, unitCell=self.unitCellA)
		self.trajStepB = tCode.TrajStepBase(time=self.timeB, step=self.stepB, unitCell=self.unitCellA)
		self.testObjA = tCode.TrajectoryInMemory([self.trajStepA, self.trajStepB])

	def testReadAndWriteConsistentA(self):
		expObj = self.testObjA
		tCode.dumpTrajObjToFile(self.testObjA, self.fileNameA)
		actObj = tCode.readTrajObjFromFileToTrajectoryInMemory(self.fileNameA)
		self.assertEqual(expObj, actObj)

	def testReadAndWriteConsistent_readLastStepOnly(self):
		expObj = self.trajStepB
		tCode.dumpTrajObjToFile(self.testObjA, self.fileNameA)
		actObj = tCode.readLastTrajStepFromFile(self.fileNameA)
		self.assertEqual(expObj,actObj)


class TestTrajInMemory(unittest.TestCase):
	
	def setUp(self):
		self.trajStepA = 4
		self.trajStepB = 5
		self.createTestObjs()

	def createTestObjs(self):
		self.fullTrajA = [self.trajStepA, self.trajStepB]
		self.testObjA = tCode.TrajectoryInMemory(self.fullTrajA)

	def testExpectedIterableReturned(self):
		for exp,act in it.zip_longest(self.fullTrajA, self.testObjA):
			self.assertEqual(exp,act)

	def testEqualObjsCompareEqual(self):
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = self.testObjA
		self.assertEqual(objA,objB)

	def testUnequalObjsCompareUnequal(self):
		objA = copy.deepcopy(self.testObjA)
		self.trajStepB += 1
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	@mock.patch("gen_basis_helpers.analyse_md.traj_core.TrajStepBase")
	def testToDictAndFromDictConsistent(self, mockedTrajStepCls):
		#Define expected
		expDictA, expDictB = mock.Mock(), mock.Mock()
		self.trajStepA, self.trajStepB = mock.Mock(), mock.Mock()

		#Set mock functions appropriately
		def _fromDictFunct(inpDict):
			if inpDict==expDictA:
				return self.trajStepA
			elif inpDict==expDictB:
				return self.trajStepB	
			else:
				raise ValueError("")

		self.trajStepA.toDict.side_effect = lambda:expDictA
		self.trajStepB.toDict.side_effect = lambda:expDictB
		mockedTrajStepCls.fromDict.side_effect = _fromDictFunct

		#Run/compare
		self.createTestObjs()
		objA = self.testObjA
		tempDict = self.testObjA.toDict()
		objB = tCode.TrajectoryInMemory.fromDict(tempDict)
		self.assertEqual(objA,objB)

	def testApplyFunctionToAllTrajSteps(self):
		mockFunct = mock.Mock()
		self.testObjA.applyFunctToEachTrajStep(mockFunct)
		for tStep in self.fullTrajA:
			mockFunct.assert_any_call(tStep)



class TestTrajStepBase(unittest.TestCase):
	
	def setUp(self):
		self.step = 4
		self.unitCell = 7 #Exactly comparable; mock is anoying since deepcopy doesnt work well for this
		self.time = 2.5
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"step":self.step, "unitCell":self.unitCell, "time":self.time}
		self.testObjA = tCode.TrajStepBase(**kwargDict)

	def testEqualObjsCompareEqual(self):
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = self.testObjA
		self.assertEqual(objA,objB)

	def testEqualObjsCompareEqual_bothTimesNone(self):
		self.time = None
		self.createTestObjs()
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = self.testObjA
		self.assertEqual(objA,objB)
		self.assertEqual(objB,objA)

	def testUnequalObjsCompareUnequal_diffStep(self):
		objA = copy.deepcopy(self.testObjA)
		self.step += 1
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA,objB)
		self.assertNotEqual(objB,objA)

	def testUnequalObjsCompareUnequal_diffTime(self):
		objA = copy.deepcopy(self.testObjA)
		self.time += 0.1
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)
		self.assertNotEqual(objB, objA)

	def testUnequalObjsCompareUnequal_oneTimeNone(self):
		objA = copy.deepcopy(self.testObjA)
		self.time = None
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)
		self.assertNotEqual(objB, objA)

	def testToDictAndFromDictConsistent(self):
		self.unitCell = uCellHelp.UnitCell(lattParams=[1,2,3], lattAngles=[90,90,90])
		self.createTestObjs()
		expObj = self.testObjA
		outDict = self.testObjA.toDict()
		actObj = tCode.TrajStepBase.fromDict(outDict)
		self.assertEqual(expObj,actObj)

class TestTrajStepFlexible(unittest.TestCase):

	def setUp(self):
		self.step = 4
		self.unitCell = 7 #Exactly comparable; mock is anoying since deepcopy doesnt work well for this
		self.time = 2.5

		self.extraAttrA = {"maxVelocity": {"value":2.5, "cmpType":"numerical"}}
		self.extraAttrB = {"velocities":  {"value":[ [2.5, 3.4, 5.4], [3.2, 5.2, 6.4] ],
		                                   "cmpType":"numericalArray"}}
		self.createTestObjs()

	def createTestObjs(self):
		self.extraAttrDict = copy.deepcopy( self.extraAttrA )
		self.extraAttrDict.update( copy.deepcopy(self.extraAttrB) )
		currKwargs = {"unitCell":self.unitCell, "step":self.step, "time":self.time,
		              "extraAttrDict":self.extraAttrDict}
		self.testObjA = tCode.TrajStepFlexible(**currKwargs)

	def testExtraAttrDictCantOverwrite(self):
		self.extraAttrA = {"numericalCmpAttrs":{}}
		with self.assertRaises(ValueError):
			self.createTestObjs()

	def testTwoUnequalObjsCompareUnequal_newNumericalAttrDiffers(self):
		objA = copy.deepcopy(self.testObjA)
		self.extraAttrA["maxVelocity"]["value"] += 1.0
		self.createTestObjs()
		self.assertNotEqual(objA, self.testObjA)
		self.assertNotEqual(self.testObjA, objA)

	def testTwoUnequalObjsCompareUnequal_newNumericalArrayDiffers(self):
		objA = copy.deepcopy(self.testObjA)
		self.extraAttrB["velocities"]["value"][1][1] += 1.0
		self.createTestObjs()
		self.assertNotEqual(objA, self.testObjA)
		self.assertNotEqual(self.testObjA, objA)

	def testTwoUnequalObjsCompareUnequal_newNumericalArrayDiffLengths2ndDim(self):
		objA = copy.deepcopy(self.testObjA)
		self.extraAttrB["velocities"]["value"][1].append(2.2)
		self.createTestObjs()
		self.assertNotEqual(objA, self.testObjA)
		self.assertNotEqual(self.testObjA, objA)

	def testTwoUnequalObjsCompareUnequal_newNumericalOneAsNone(self):
		objA = copy.deepcopy(self.testObjA)
		self.extraAttrB["velocities"]["value"] = None
		self.createTestObjs()
		self.assertNotEqual(objA, self.testObjA)
		self.assertNotEqual(self.testObjA, objA)

	#Note this test should fail if cmpAttrs arent set properly with to/from dict
	def testToDictAndFromDictConsistent(self):
		#We expect more keys than this; but want to test at least these are present
		self.unitCell = uCellHelp.UnitCell(lattParams=[1,2,3], lattAngles=[90,90,90])
		self.createTestObjs()
		expKeysPartial = {"maxVelocity", "velocities"} 
		outDict = self.testObjA.toDict()
		objB = tCode.TrajStepFlexible.fromDict(outDict)

		for key in expKeysPartial:
			self.assertTrue( key in outDict.keys() )

		self.assertEqual(self.testObjA, objB)
		self.assertEqual(objB, self.testObjA)

	def testAddAttrsAfterInit(self):
		extraAttrDict = copy.deepcopy(self.extraAttrB)
		expObj = copy.deepcopy(self.testObjA)
		self.extraAttrB = dict()
		self.createTestObjs()
		self.testObjA.addExtraAttrDict(extraAttrDict)
		self.assertEqual(expObj, self.testObjA)



#Tough case:
# [  [1,3], [4,8], [5,7] ] -> [ [1,3], [4,8] ] BUT very tough to pull off
# Alternatively, I could just raise in that case
class TestMergeTrajInMemory(unittest.TestCase):

	def setUp(self):
		self.stepsA = [1,2,3]
		self.stepsB = [4,5,6]
		self.stepsC = [7,8,9]
		self.createTestObjs()

	def createTestObjs(self):
		self.trajA = tCode.TrajectoryInMemory( [tCode.TrajStepBase(step=x) for x in self.stepsA] )
		self.trajB = tCode.TrajectoryInMemory( [tCode.TrajStepBase(step=x) for x in self.stepsB] )
		self.trajC = tCode.TrajectoryInMemory( [tCode.TrajStepBase(step=x) for x in self.stepsC] )
		self.trajListA = [self.trajA, self.trajB, self.trajC]

		allSteps = self.stepsA + self.stepsB + self.stepsC
		self.expTrajStd = tCode.TrajectoryInMemory( [tCode.TrajStepBase(step=x) for x in allSteps] )

	def testExpectedFromMergingWhenAlreadyOrdered(self):
		expTraj = self.expTrajStd
		actTraj = tCode.getMergedTrajInMemory(self.trajListA)
		self.assertEqual(expTraj,actTraj)

	def testExpectedWhenMergingOutOfOrderCases(self):
		self.trajListA = [self.trajB, self.trajA, self.trajC]
		expTraj = self.expTrajStd
		actTraj = tCode.getMergedTrajInMemory(self.trajListA)
		self.assertEqual(expTraj, actTraj)

	def testRaisesIfTrajectoriesOverlap_overlapStratNone_trimStratNone(self):
		self.stepsA.append(self.stepsB[-1])
		self.createTestObjs()
		with self.assertRaises(ValueError):
			tCode.getMergedTrajInMemory(self.trajListA, overlapStrat=None, trimStrat=None)

	def testExpectedStepsWithUnitOverlap_overlapStratSimple(self):
		self.stepsA.append(4)
		self.createTestObjs()
		expSteps = [1,2,3,4,5,6,7,8,9]
		expTraj = tCode.TrajectoryInMemory( [tCode.TrajStepBase(step=x) for x in expSteps] )
		actTraj = tCode.getMergedTrajInMemory(self.trajListA, overlapStrat="simple")
		self.assertEqual(expTraj,actTraj)

	def testExpectedStepsWithOverhang_trimStratSimple(self):
		self.stepsA = [1,2,3,4,5]
		self.stepsB = [4,6] #Deleting step 5 makes it more likely we really are taking the step 4 from this set of trajs
		self.stepsC = [7,8,9]
		self.createTestObjs()
		expSteps = [1,2,3,4,6,7,8,9]
		expTraj = tCode.TrajectoryInMemory( [tCode.TrajStepBase(step=x) for x in expSteps] )
		actTraj = tCode.getMergedTrajInMemory(self.trajListA, overlapStrat="simple", trimStrat="simple")

		self.assertEqual(expTraj,actTraj)


class TestGetGeomClosestToTime(unittest.TestCase):

	def setUp(self):
		self.equiDistTol = 1e-2
		self.prioritiseLate = True
		self.quenchTime = 5
		self.createTestObjs()

	def createTestObjs(self):
		self.testTrajA = _loadTrajInMemoryA()
		
	def _runTestFunct(self):
		kwargs = {"equiDist":self.equiDistTol, "prioritiseLate":self.prioritiseLate}
		return tCode.getTimeAndGeomClosestToInpTimeForInpTraj(self.quenchTime, self.testTrajA, **kwargs)


	@mock.patch("gen_basis_helpers.analyse_md.traj_core.copy.deepcopy")
	def testPrioritiseLate(self, mockCopy):
		mockCopy.side_effect = lambda x: x
		expTime = 6
		expGeom = self.testTrajA.trajSteps[2].unitCell
		actTime,actGeom = self._runTestFunct()
		self.assertAlmostEqual(expTime, actTime)
		self.assertAlmostEqual(expGeom, actGeom)


	@mock.patch("gen_basis_helpers.analyse_md.traj_core.copy.deepcopy")
	def testSimplePrioritiseEarlyA(self, mockCopy):
		self.prioritiseLate = False
		self.createTestObjs()
		mockCopy.side_effect = lambda x: x
		expTime = 4
		expGeom = self.testTrajA.trajSteps[1].unitCell
		actTime,actGeom = self._runTestFunct()
		self.assertAlmostEqual(expTime, actTime)
		self.assertAlmostEqual(expGeom, actGeom)

	@mock.patch("gen_basis_helpers.analyse_md.traj_core.copy.deepcopy")
	def testEquidistHasExpectedEffect(self, mockCopy):
		self.equiDistTol = 0.3
		self.prioritiseLate = False
		self.quenchTime = 5.1 #Closer 6, but eqidist makes t=4 and t=6 seem equally distanced from it (abs(1.1-0.9)<0.3)
		mockCopy.side_effect = lambda x: x
		self.createTestObjs()
		expTime = 4
		expGeom = self.testTrajA.trajSteps[1].unitCell
		actTime, actGeom = self._runTestFunct()
		self.assertAlmostEqual(expTime, actTime)
		self.assertAlmostEqual(expGeom, actGeom)

	@mock.patch("gen_basis_helpers.analyse_md.traj_core.copy.deepcopy")
	def testTimeZero(self, mockCopy):
		self.quenchTime = 0
		mockCopy.side_effect = lambda x: x
		self.createTestObjs()
		expTime = 2
		expGeom = self.testTrajA.trajSteps[0].unitCell
		actTime, actGeom = self._runTestFunct()
		self.assertAlmostEqual(expTime, actTime)
		self.assertAlmostEqual(expGeom, actGeom)



def _loadTrajInMemoryA():
	outStepVals = [0,1,2,3]
	outTimeVals = [2,4,6,8]
	outGeomVals = [mock.Mock() for x in outStepVals] 

	outIter = [ [a,b,c] for a,b,c in it.zip_longest(outGeomVals, outStepVals, outTimeVals) ] 
	outTrajSteps = [tCode.TrajStepBase(unitCell=cell, step=step, time=time) for cell,step,time in outIter]

	return tCode.TrajectoryInMemory(outTrajSteps)

