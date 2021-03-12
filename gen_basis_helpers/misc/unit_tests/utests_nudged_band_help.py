


import copy
import unittest
import unittest.mock as mock

import gen_basis_helpers.misc.nudged_band_paths as tCode


class TestNudgedBandIO(unittest.TestCase):

	def setUp(self):
		self.outPathA = "fake_path_a"
		self.objA = mock.Mock()

	@mock.patch("gen_basis_helpers.misc.nudged_band_paths.sharedIoHelp.readObjWithFromDictFromJsonFile")
	def testExpectedCallMadeForReadingPath(self, mockReadFunct):
		tCode.readNudgedBandPathwayFromJsonFileStandard(self.outPathA)
		mockReadFunct.assert_called_with(tCode.NudgedBandPathStandard, self.outPathA)

	@mock.patch("gen_basis_helpers.misc.nudged_band_paths.sharedIoHelp.dumpObjWithToDictToJson")
	def testExpecetedCallMadeForWritingToPath(self, mockDumpFunct):
		tCode.dumpNudgedBandPathwayToFile(self.objA, self.outPathA)
		mockDumpFunct.assert_called_with(self.objA, self.outPathA)


class TestNudgedBandPath(unittest.TestCase):

	def setUp(self):
		self.distA = 4
		self.distB = 5
		self.distC = 6
		self.createTestObjs()

	def createTestObjs(self):
		self.stepA = tCode.NudgedBandStepStandard(dist=self.distA)
		self.stepB = tCode.NudgedBandStepStandard(dist=self.distB)
		self.stepC = tCode.NudgedBandStepStandard(dist=self.distC)

		self.pathAB  = tCode.NudgedBandPathStandard([self.stepA, self.stepB])
		self.pathABC = tCode.NudgedBandPathStandard([self.stepA, self.stepB, self.stepC])
		self.pathAC  = tCode.NudgedBandPathStandard([self.stepA, self.stepC])

	def testEqualObjsCompareEqual(self):
		objA = copy.deepcopy(self.pathAB)
		self.createTestObjs()
		objB = self.pathAB
		self.assertEqual(objA,objB)

	def testUnequalObjsCompareUnequal_extraStep(self):
		self.assertNotEqual(self.pathAB, self.pathABC)

	def testUnequalObjsCompareUnequal_secondStepDifferent(self):
		self.assertNotEqual(self.pathAB, self.pathAC)

	def testToAndFromDictConsistent(self):
		inpDict = self.pathAB.toDict()
		objB = tCode.NudgedBandPathStandard.fromDict(inpDict)
		self.assertEqual(self.pathAB,objB)

class TestNudgedBandStep(unittest.TestCase):

	def setUp(self):
		self.geom = 5 #Would be geom object; equality works basically same for integer though
		self.energies = 7 #Would be energies object; equality works basically same for integer though
		self.dist = 1.5
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"geom":self.geom,"energies":self.energies,"dist":self.dist}
		self.testObjA = tCode.NudgedBandStepStandard(**kwargDict)

	def testEqualObjsCompareEqualA(self):
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = self.testObjA
		self.assertEqual(objA, objB)

	def testUnequlObjsCompareUnequal_diffDist(self):
		objA = copy.deepcopy(self.testObjA)
		self.dist += 0.6
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA,objB)

	def testUnequalObjsCompareUnequal_diffGeom(self):
		objA = copy.deepcopy(self.testObjA)
		self.energies += 1
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testUnequalObjsCompareUnequal_geomNone(self):
		objA = copy.deepcopy(self.testObjA)
		self.geom = None
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testUnequalObjsCompareUnequal_distNone(self):
		objA = copy.deepcopy(self.testObjA)
		self.dist = None
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)

	def testToAndFromDictConsistent(self):
		self.geom, self.energies = mock.Mock(), mock.Mock()
		self.geom.toDict.side_effect = {"key_geom":"val_geom"}
		self.energies.toDict.side_effect = {"key_energies":"val_energies"}
		self.createTestObjs()
		inpDict = self.testObjA.toDict()
		objB = tCode.NudgedBandStepStandard.fromDict(inpDict)

		self.geom.toDict.assert_called_with()
		self.energies.toDict.assert_called_with()

		self.assertEqual(self.testObjA,objB)



