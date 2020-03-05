
import os

import unittest
import unittest.mock as mock

import gen_basis_helpers.castep.castep_creator as tCode

class TestCastepCalcObj(unittest.TestCase):

	def setUp(self):
		self.basePath = os.path.join("fake","folder","fake_file")
		self.cellFileDict = {"fake_kwarg":"fake_val"}
		self.paramFileDict = {"fake_param":"fake_param_kwarg"}
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.CastepCalcObj(self.basePath, self.paramFileDict, self.cellFileDict)

	@mock.patch("gen_basis_helpers.castep.castep_creator.parseCastep")
	def testCorrectParamFileStrDictPassedToWriter(self, mockedParserModule):
		self.testObjA.writeFile()
		expParamFile = self.basePath + ".param"
		mockedParserModule.writeCastepParamFileFromDict.assert_called_once_with(expParamFile, self.paramFileDict)

	@mock.patch("gen_basis_helpers.castep.castep_creator.parseCastep")
	def testCorrectCellFileStrPassedToWriter(self, mockedParserModule):
		self.testObjA.writeFile()
		expCellFile = self.basePath + ".cell"
		mockedParserModule.writeCastepCellFileFromTokens.assert_called_once_with(expCellFile,self.cellFileDict)

	def testBasePathGetterAndSetterIgnoresExtensions(self):
		expBaseFile = "fileWithoutExt"
		self.testObjA.basePath = expBaseFile + ".test_extension"
		actBaseFile = self.testObjA.basePath
		self.assertEqual(expBaseFile, actBaseFile)

	@mock.patch("gen_basis_helpers.castep.castep_creator.parseCastep")
	def testParsedFileCallsCorrectParser(self, mockedParserModule):
		mockedParserModule.parseCastepOutfile.side_effect = lambda *args,**kwargs: {"unitCell":mock.Mock()}
		unused = self.testObjA.parsedFile
		expCastepFile = self.basePath + ".castep"
		mockedParserModule.parseCastepOutfile.assert_called_once_with(expCastepFile)








