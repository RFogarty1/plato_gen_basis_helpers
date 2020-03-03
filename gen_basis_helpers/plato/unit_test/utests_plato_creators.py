
import unittest
import unittest.mock as mock

import gen_basis_helpers.plato.plato_creator as tCode

class TestPlatoCalcObj(unittest.TestCase):

	def setUp(self):
		self.basePath = "fake/path/to/nowhere"
		self.strDict = mock.Mock()
		self.pathToRunCommFunction = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.PlatoCalcObj(self.basePath, self.strDict, self.pathToRunCommFunction)

	def testRunCommAsExpected(self):
		testCommFormat = "fake command with path: {} included"
		expCommand = testCommFormat.format(self.basePath+".in")
		self.pathToRunCommFunction.side_effect = lambda x: testCommFormat.format(self.basePath+".in")
		actCommand = self.testObjA.runComm
		self.assertEqual(expCommand,actCommand)

	def testBasePathRemovesExtensionAtInit(self):
		self.basePath="some_file.with_extension"
		expBasePath = "some_file"
		self.createTestObjs()
		actPath = self.testObjA.basePath
		self.assertEqual(expBasePath,actPath)

	@mock.patch("gen_basis_helpers.plato.plato_creator.modPlatoInp.writePlatoOutFileFromDict")
	def testWriteFileCallsExpectedFunct(self, mockedWriter):
		expPath = self.basePath + ".in"
		self.testObjA.writeFile()
		mockedWriter.assert_called_once_with(expPath,self.strDict)

	def testOutFilePathAsExpected(self):
		expPath = self.basePath + ".out"
		actPath = self.testObjA.outFilePath
		self.assertEqual(expPath, actPath)

	@mock.patch("gen_basis_helpers.plato.plato_creator.parsePlatoOut")
	def testParsedFileCallsCorrectParser(self,mockedParser):
		parsedFile = self.testObjA.parsedFile
		mockedParser.parsePlatoOutFile_energiesInEv.assert_called_once_with( self.basePath + ".out" )



class TestPlatoMethodObjCreator(unittest.TestCase):

	def setUp(self):
		pass

	@unittest.skip("")
	def testSomething(self):
		self.assertTrue(False)



