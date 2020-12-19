
import collections
import unittest
import unittest.mock as mock

import gen_basis_helpers.lammps_interface.misc_objs as tCode

class TestNVTEnsemble(unittest.TestCase):

	def setUp(self):
		self.startTemp = 300
		self.finalTemp = 500
		self.dampTime = 200
		self.thermostat = "Nose-Hoover"
		self.numbFmt = "{:.1f}"
		self.createTestObjs()

	def createTestObjs(self):
		currKwargs = {"thermostat":self.thermostat, "endTemp":self.finalTemp, "dampTime":self.dampTime,
		              "numbFmt":self.numbFmt}
		self.testObjA = tCode.NVTEnsemble(self.startTemp, **currKwargs)

	def testExpStrFromSimpleOptsA(self):
		expStr = "all nvt temp 300.0 500.0 200.0"
		actStr = self.testObjA.fixStr
		self.assertEqual(expStr, actStr)

	def testRaisesIfDampTimeNotSet(self):
		self.dampTime = None
		self.createTestObjs()
		with self.assertRaises(ValueError):
			self.testObjA.fixStr


class TestNPTEnsemble(unittest.TestCase):

	def setUp(self):
		self.startTemp = 300
		self.finalTemp = 500
		self.dampTimeTemp = 200
		self.pressureDims = "z"
		self.startPressure  = 10
		self.endPressure = 20
		self.dampTimePressure = 2000
		self.numbFmtTemp = "{:.2f}"
		self.numbFmtPressure = "{:.2f}"
		self.numbFmtTime = "{:.2f}"
		self.numbFmtAll = "{:.2f}"
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"pressureDims":self.pressureDims, "endTemp":self.finalTemp,
		             "endPressure":self.endPressure, "dampTimeTemp":self.dampTimeTemp,
		             "dampTimePressure":self.dampTimePressure, "numbFmtTime":self.numbFmtTime,
		             "numbFmtPressure":self.numbFmtPressure, "numbFmtTemp":self.numbFmtTemp}
		self.testObjA = tCode.NPTEnsembleNoseHooverStandard(self.startTemp, self.startPressure, **kwargDict)

	def testExpStrFromSimpleOptsA(self):
		expStr = "all npt temp 300.00 500.00 200.00 z 10.00 20.00 2000.00"
		actStr = self.testObjA.fixStr
		self.assertEqual(expStr, actStr)

	def testRaiseIfPressureDampTimeNotSet(self):
		self.dampTimePressure = None
		self.createTestObjs()
		with self.assertRaises(ValueError):
			self.testObjA.fixStr

	def testRaisesIfTempDampTimeNotSet(self):
		self.dampTimeTemp = None
		self.createTestObjs()
		with self.assertRaises(ValueError):
			self.testObjA.fixStr


class TestCreateVelocityObj(unittest.TestCase):

	def setUp(self):
		self.temp = 400
		self.seed = 300
		self.dist = "fake_dist"
		self.group = "fake_group"
		self.createTestObjs()

	def createTestObjs(self):
		currKwargs = {"seed":self.seed, "group":self.group, "dist":self.dist}
		self.testObjA = tCode.VelocityCreateCommObj(self.temp, **currKwargs)

	def testExpectedStrFromSimpleOptsA(self):
		expStr = "fake_group create 400.0 300 dist fake_dist"
		actStr = self.testObjA.commandStr
		self.assertEqual(expStr,actStr)


class TestDumpObjStandard(unittest.TestCase):

	def setUp(self):
		self.everyNSteps = 30
		self.groupId = "fake_group"
		self.dumpType = "atom"
		self.fileExt = "lammpstrj"
		self.scale = True
		self.createTestObjs()

	def createTestObjs(self):
		currKwargs = {"groupId":self.groupId, "dumpType":self.dumpType, "fileExt":self.fileExt,
		              "scale":self.scale}
		self.testObjA = tCode.DumpCommObjStandard(self.everyNSteps, **currKwargs)

	def testExpectedDictFromSimpleOptsA(self):
		currArgs = [ ["dump","myDump fake_group atom 30 dump.lammpstrj"],
		             ["dump_modify", "myDump scale yes"] ]
		expDict = collections.OrderedDict(currArgs)
		actDict = self.testObjA.commandDict
		self.assertEqual(expDict,actDict)


class TestReflectiveWallFaceObj(unittest.TestCase):

	def setUp(self):
		self.face = "xlo"
		self.groupId = "fake_group"
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.ReflectiveWallFace(self.face, groupId=self.groupId)

	def testExpectedFixCommandFromDictA(self):
		expStr = "fake_group wall/reflect xlo EDGE"
		actStr = self.testObjA.fixStr
		self.assertEqual(expStr,actStr)



