
import os
import unittest
import unittest.mock as mock

import gen_basis_helpers.ref_data.targ_dir_from_path as tCode

class TestGetTargDirFromPath(unittest.TestCase):

	def setUp(self):
		self.baseDirA = os.path.join("basedir","extraBit")
		self.targSubDir = "fake_sub_dirA"
		self.testInpPathExtraBit = "test_path_a"
		self.testInpPathA = os.path.join(self.baseDirA,self.testInpPathExtraBit)
		self.createTestObjs()

	@mock.patch("gen_basis_helpers.ref_data.targ_dir_from_path.os.path.abspath")
	def createTestObjs(self, mockedAbsPath):
		mockedAbsPath.side_effect = lambda inpPath: inpPath
		self.testObjA = tCode.GetSpecialDirFromSubDir(self.baseDirA, self.targSubDir)

	#Note i cant mock the os.path.join part, because other os.path bits depend on it
	def testExpCallMade_noKwargs(self):
		expArgs = [self.baseDirA, self.targSubDir, self.testInpPathExtraBit]
		expPath = os.path.join(self.baseDirA, self.targSubDir, self.testInpPathExtraBit)
		actPath = self.testObjA(self.testInpPathA)
		self.assertEqual(expPath, actPath)

	def testRaisesIfCurrDirNotSubdir(self):
		self.testInpPathA = "fake_dir"
		with self.assertRaises(AssertionError):
			self.testObjA(self.testInpPathA)

	def testExtensionAdded(self):
		extension="extension"
		expPath = os.path.join(self.baseDirA, self.targSubDir, self.testInpPathExtraBit, extension)
		actPath = self.testObjA(self.testInpPathA, extension=extension)
		self.assertEqual(expPath,actPath)

