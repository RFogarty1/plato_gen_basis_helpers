
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

