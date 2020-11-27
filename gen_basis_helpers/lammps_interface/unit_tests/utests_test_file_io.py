
import collections

import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.lammps_interface.lammps_geom as geomHelp
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


class TestGetUnitCellObjFromDataFile(unittest.TestCase):


	def setUp(self):
		self.lattParamsA = [2,3,4]
		self.lattAnglesA = [90,90,90]
		self.fractCoordsA = [ [0.5,0.5,0.5,"X"],
		                      [0.7,0.7,0.7,"Y"],
		                      [0.7,0.7,0.7,"Y"] ]
		self.atomStyle = "full"
		self.massDict = {"X":4.5,"Y":2.8}
		self.createTestObjs()

	def createTestObjs(self):
		kwargs = {"lattParams":self.lattParamsA, "lattAngles":self.lattAnglesA}
		self.testCellA = uCellHelp.UnitCell(**kwargs)
		self.testCellA.fractCoords = self.fractCoordsA
		self.lammpsObjA = self._createLammpsGeomFromUCell(self.testCellA)
		self.geomToDataDict = geomHelp.GetDataDictFromLammpsGeomAtomStyleFull()

	def _createLammpsGeomFromUCell(self, inpCell):
		kwargs = {"eleToTypeIdx":{"X":1,"Y":2}, "eleToCharge":{"X":0,"Y":0}, "eleToMass":self.massDict,
		          "geomToBondInfo": lambda *args:list(), "geomToAngleInfo": lambda *args:list(),
		          "geomToMoleculeIDs":lambda *args:list()} 
		outObj = geomHelp.LammpsGeom(inpCell, **kwargs)
		return outObj

	@mock.patch("gen_basis_helpers.lammps_interface.file_io.tokenizeDataFile")
	def testExpectedGeomFromFileA(self, mockTokenizeDataFile):
		fakeFilePath = "fake_path"
		dataDict = self.geomToDataDict(self.lammpsObjA)
		mockTokenizeDataFile.side_effect = lambda arg: dataDict
		expCell = self.testCellA
		actCell = tCode.getUCellObjectFromDataFile(fakeFilePath, atomStyle=self.atomStyle, massDict=self.massDict)
		mockTokenizeDataFile.assert_called_with(fakeFilePath)
		self.assertEqual(expCell,actCell)


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

