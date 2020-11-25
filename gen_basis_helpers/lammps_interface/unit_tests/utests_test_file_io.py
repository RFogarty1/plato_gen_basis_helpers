
import collections

import unittest
import unittest.mock as mock

import gen_basis_helpers.lammps_interface.file_io as tCode

class TestTokenizeDataFile(unittest.TestCase):

	def setUp(self):
		self.createTestObjs()

	def createTestObjs(self):
		self.fileStrA = _loadDataFileStrA()

	@mock.patch("gen_basis_helpers.lammps_interface.file_io._getStrFromInpPath")
	def testExpectedTokensFromDudFileA(self, mockedGetStrFromInpPath):
		mockedGetStrFromInpPath.side_effect = lambda x: _loadDataFileStrA()
		testInpPath = "fake_path"
		expDict = _getExpectedDictFromFileStrA()
		actDict = tCode.tokenizeDataFile(testInpPath)
		mockedGetStrFromInpPath.assert_called_with(testInpPath)
		self.assertEqual(expDict,actDict)

	@mock.patch("gen_basis_helpers.lammps_interface.file_io._writeFileFromStr")
	def testExpectedFileStrFromTokens_simpleInput(self, mockedWriteToFile):
		expStr = _loadDataFileStrA()
		expPath, tokens = "fake_path", _getExpectedDictFromFileStrA()
		tCode.writeDataFileFromTokens(expPath, tokens)
		mockedWriteToFile.assert_called_with(expPath, expStr)


class TestWriteScriptFile(unittest.TestCase):

	def setUp(self):
		self.createTestObjs()

	def createTestObjs(self):
		self.fileStrA = _loadScriptFileStrA()

	@mock.patch("gen_basis_helpers.lammps_interface.file_io._writeFileFromStr")
	def testWriteScriptFileA(self, mockedWriteToFile):
		expPath, expStr = "fake_path_a", _loadScriptFileStrA()
		inpList = [ ["units","real"], ["atom_style","atomic\n"], ["boundary","p p p"] ]
		tokenDict = collections.OrderedDict(inpList)
		tCode.writeScriptFileFromTokens(expPath, tokenDict)
		mockedWriteToFile.assert_called_with(expPath, expStr)


def _loadDataFileStrA():
	outStr = """LAMMPS Atom File

         4500 atoms
         3000 bonds

Masses

	1	15.9994

Atoms

      1    1  1 -0.8472   12.12456   28.09298   22.27452  0  1  0

"""
	return outStr


def _getExpectedDictFromFileStrA():
	headerStr = "         4500 atoms\n         3000 bonds"
	massStr = "	1	15.9994"
	atomStr = "      1    1  1 -0.8472   12.12456   28.09298   22.27452  0  1  0"
	outDict = collections.OrderedDict([ ["LAMMPS Atom File",headerStr], ["Masses",massStr], ["Atoms",atomStr] ])
	return outDict


def _loadScriptFileStrA():
	outStr = """units real
atom_style atomic
boundary p p p
"""
	return outStr

