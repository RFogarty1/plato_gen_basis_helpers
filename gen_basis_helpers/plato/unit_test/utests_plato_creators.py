
import os
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
		self.workFolder = "fake/work/folder"
		self.fileName = "fake_filename"
		self.kPts = [20,20,20]
		self.gridVals = [50,40,40] #The method itself determines which type of grid we use
		self.methodStr = "dft2_fake_method_key" #The start actually determines how some kwargs are set
		self.dataSet = "fakeDataSet"
		self.geom = "fake_geom"
		self.createTestObjects()

	def createTestObjects(self):
		self.testObjA = tCode.PlatoCalcObjFactoryStandard( workFolder=self.workFolder, fileName=self.fileName,
		                                                   kPts=self.kPts, methodStr=self.methodStr,
		                                                   gridVals=self.gridVals, dataSet=self.dataSet, geom=self.geom )

	@mock.patch("gen_basis_helpers.plato.plato_creator.calcMeth")
	def testBaseOptDictGetterCalled(self, mockedCalcMethods):
		self.testObjA.create()
		mockedCalcMethods.createPlatoMethodObj.assert_called_with(self.methodStr)

	@mock.patch("gen_basis_helpers.plato.plato_creator.PlatoCalcObjFactoryStandard._getMethodObj")
	def testExpectedValuesSetOnMethodObj(self, mockedGetMethodObj):
		methodObjMock = mock.Mock()
		mockedGetMethodObj.side_effect = lambda : methodObjMock
		self.testObjA.create()
		self.assertEqual(self.kPts,methodObjMock.kpts)
		self.assertEqual(self.gridVals, methodObjMock.integGrid)
		self.assertEqual(self.dataSet, methodObjMock.dataSet)

	@mock.patch("gen_basis_helpers.plato.plato_creator.PlatoCalcObjFactoryStandard._getMethodObj")
	def testGeomPassedToStrDictCreator(self, mockedGetMethodObj):
		methodObjMock = mock.Mock()
		mockedGetMethodObj.side_effect = lambda : methodObjMock
		self.testObjA.create()
		methodObjMock.getStrDictWithStruct.assert_called_with(self.geom)


	@mock.patch("gen_basis_helpers.plato.plato_creator.PlatoCalcObjFactoryStandard._getStrDictFromMethodObj")
	@mock.patch("gen_basis_helpers.plato.plato_creator.PlatoCalcObjFactoryStandard._getMethodObj")
	@mock.patch("gen_basis_helpers.plato.plato_creator.PlatoCalcObj")
	def testExpectedCallToPlatoCalcObj(self, mockedCalcObj, mockedGetMethodObj, mockedGetStrDict):
		methodObjMock = mock.Mock()
		strDictMock = mock.Mock()
		fakeOutputObj = "fake_output_obj"
		methodObjMock.runCommFunction = "fake_run_comm_function"
		mockedGetMethodObj.side_effect = lambda : methodObjMock
		mockedGetStrDict.side_effect = lambda x: strDictMock
		mockedCalcObj.side_effect = lambda *args : fakeOutputObj

		fakeStrDict = "fakeStrDict"
		expBasePath = os.path.join(self.workFolder, self.fileName)
		expStrDict = strDictMock
		expRunCommFunction = methodObjMock.runCommFunction

		actOutputObj = self.testObjA.create()
		mockedCalcObj.assert_called_once_with(expBasePath, expStrDict, expRunCommFunction)
		self.assertEqual(fakeOutputObj, actOutputObj)


