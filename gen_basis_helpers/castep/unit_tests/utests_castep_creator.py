
import os

import unittest
import unittest.mock as mock

import gen_basis_helpers.castep.castep_creator as tCode



class TestCastepCreator(unittest.TestCase):

	def setUp(self):
		self.testMethodStrA = "castep_test_object"
		self.workFolder = os.path.join("fake","folder")
		self.fileName = "fake_filename"
		self.kPts = [20,20,12]
		self.mgPot = "mg_pseudo"
		self.zrPot = "zr_pseudo"
		self.cutoffEnergy = 700
		self.symmetryGenerate = True
		self.geom = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		pseudoDict = {"Zr":self.zrPot, "Mg":self.mgPot}
		self.testObjA = tCode.CastepCalcObjFactoryStandard(methodStr=self.testMethodStrA, geom=self.geom, kPts=self.kPts,
		                                                  workFolder=self.workFolder, fileName=self.fileName,
		                                                  symmetryGenerate=self.symmetryGenerate, pseudoPotDict=pseudoDict,
		                                                  cutoffEnergy=self.cutoffEnergy)

	@mock.patch("gen_basis_helpers.castep.castep_creator.methodReg")
	def testMethodStrPassedOn(self,mockedMethodReg):
		unused = self.testObjA.paramFileDict
		mockedMethodReg.createParamDictFromMethodStr.assert_called_once_with(self.testMethodStrA)

	def testBasePathCorrect(self):
		expBasePath = os.path.join(self.workFolder, self.fileName)
		actBasePath = self.testObjA.basePath
		self.assertEqual(expBasePath,actBasePath)

	@mock.patch("gen_basis_helpers.castep.castep_creator.parseCastep")
	def testCreateCellFileDict(self, mockedParser):
		expLattCart, expPosFract = mock.Mock(), mock.Mock()

		mockedParser.getCellGeomDictSectionFromUCell.side_effect = lambda x, **kwargs: { "lattice_cart"  : expLattCart,
		                                                                                 "positions_frac": expPosFract }

		expSpeciesPot = "{} {}\n{} {}".format("Mg", self.mgPot, "Zr", self.zrPot) #Alphabetical force to make it easier to test
		expDict = {"kpoint_mp_grid":"{} {} {}".format(*self.kPts),
		           "lattice_cart":expLattCart,
		           "positions_frac":expPosFract,
		           "species_pot":expSpeciesPot,
		           "symmetry_generate":""}
		actDict = self.testObjA.cellFileDict 

		mockedParser.getCellGeomDictSectionFromUCell.assert_called_once_with(self.geom)		
		for key in expDict:
			self.assertEqual( expDict[key], actDict[key] )

	def testCreateParamFileDict(self):
		expModdedDict = {"cut_off_energy":str(self.cutoffEnergy)}
		actFullDict = self.testObjA.paramFileDict

		for key in expModdedDict:
			self.assertEqual(expModdedDict[key], actFullDict[key])



	@mock.patch("gen_basis_helpers.castep.castep_creator.CastepCalcObj")
	@mock.patch("gen_basis_helpers.castep.castep_creator.CastepCalcObjFactoryStandard.basePath",new_callable=mock.PropertyMock)
	@mock.patch("gen_basis_helpers.castep.castep_creator.CastepCalcObjFactoryStandard.cellFileDict",new_callable=mock.PropertyMock)
	@mock.patch("gen_basis_helpers.castep.castep_creator.CastepCalcObjFactoryStandard.paramFileDict",new_callable=mock.PropertyMock)
	def testPassesCorrectArgsToCalcObj(self, mockedParamFileDict, mockedCellFileDict, mockedBasePath, mockedCalcObj):
		expParamDict, expCellDict, expBasePath, expCalcObj = mock.Mock(), mock.Mock(), mock.Mock(), mock.Mock()
		mockedParamFileDict.return_value = expParamDict
		mockedCellFileDict.return_value = expCellDict
		mockedBasePath.return_value = expBasePath
		mockedCalcObj.side_effect = lambda *args,**kwargs: expCalcObj
		
		actCalcObj = self.testObjA.create()
		mockedCalcObj.assert_called_once_with(expBasePath,expParamDict,expCellDict)
		self.assertEqual(expCalcObj,actCalcObj)




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







