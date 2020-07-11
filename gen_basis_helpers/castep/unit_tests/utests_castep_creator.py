
import os

import unittest
import unittest.mock as mock


import gen_basis_helpers.shared.geom_constraints as geomConstraints
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
		self.geom.cartCoords = list()
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

	def testGeomOptPassedToParamDict(self):
		""" Check that setting runType as a geometry optimisation works """
		self.testObjA.runType = "geomOpt" #Case insensitive
		expParamDictKeyVals = {"task":"GeometryOptimization".lower()}
		actParamDict = self.testObjA.paramFileDict
		for key in expParamDictKeyVals.keys():
			self.assertEqual( expParamDictKeyVals[key], actParamDict[key] )


	@mock.patch("gen_basis_helpers.castep.castep_creator.getIonicConstraintsStrFromAtomicPosConstraintObj")
	@mock.patch("gen_basis_helpers.castep.castep_creator.parseCastep")
	@mock.patch("gen_basis_helpers.castep.castep_creator.getCellConstraintsStrFromGeomConstraintsObj")
	@mock.patch("gen_basis_helpers.shared.geom_constraints.GeomConstraints")
	def testLackOfCellConstraintsPassedToCellDict(self, mockedGeomConstraints, mockedGetConstraintsStr, mockedParser, mockedIonicConstr):
		""" Check that by default we have no cell constraints """
		expConstraintsObj, expConstraintsStr = mock.Mock(), mock.Mock()

		mockedGeomConstraints.initWithNoConstraints.side_effect = [expConstraintsObj]
		mockedGetConstraintsStr.side_effect = [expConstraintsStr]
		expCellDictKeyVals = {"cell_constraints":expConstraintsStr}

		actCellDict = self.testObjA.cellFileDict
		mockedGetConstraintsStr.assert_called_once_with(expConstraintsObj)
		for key in expCellDictKeyVals.keys():
			self.assertEqual( expCellDictKeyVals[key], actCellDict[key] )




class TestGetCellConstraintsFromObj(unittest.TestCase):

	def setUp(self):
		self.anglesToFix = [False, False, False]
		self.lattParamsToFix = [False,False,False]
		self.createTestObjs()

	def createTestObjs(self):
		self.cellConstrA = geomConstraints.CellConstraints(self.anglesToFix, self.lattParamsToFix)
		self.atomicPosConstraints = mock.Mock()
		self.testObjA = geomConstraints.GeomConstraints(self.atomicPosConstraints,self.cellConstrA)

	def _testExpectedStrGiven(self, expStr):
		actStr = tCode.getCellConstraintsStrFromGeomConstraintsObj(self.testObjA)
		self.assertEqual(expStr,actStr)

	def testNoConstraintsCase(self):
		expStr = "1 2 3\n4 5 6"
		self._testExpectedStrGiven(expStr)
		
	def testAllConstrainedCase(self):
		self.anglesToFix = [True,True,True]
		self.lattParamsToFix = [True,True,True]
		self.createTestObjs()
		expStr = "0 0 0\n0 0 0"
		self._testExpectedStrGiven(expStr)

	def testTwoParamsConstrainedCase(self):
		self.lattParamsToFix[0] = True
		self.lattParamsToFix[2] = True
		self.createTestObjs()
		expStr = "0 1 0\n2 3 4"
		self._testExpectedStrGiven(expStr)

	def testAllAnglesConstrainedCase(self):
		self.anglesToFix = [True, True, True]
		self.createTestObjs()
		expStr = "1 2 3\n0 0 0"
		self._testExpectedStrGiven(expStr)


class TestIonicConstraintsFromObj(unittest.TestCase):

	def setUp(self):
		self.eleList = ["X","Y","X"]
		self.atomicConstraintObjs = list()
		self.createTestObjs()

	def createTestObjs(self):
		self.geoConstraints = geomConstraints.GeomConstraints.initWithNoConstraints()
		self.geoConstraints.atomicPostionConstraints.atomicCartConstraints = self.atomicConstraintObjs

	def testReturnsNoneForNoConstraintsCase(self):
		expVal = None
		actVal = tCode.getIonicConstraintsStrFromAtomicPosConstraintObj(self.geoConstraints.atomicPostionConstraints, self.eleList)
		self.assertEqual(expVal,actVal)

	def testReturnsExpectedForSingleAtomCoordFixed(self):
		atomIdxFixed = 1 #Note indexing starts at zero
		atomConstraint = geomConstraints.AtomicCartesianConstraint(atomIdxFixed, fixX=True)
		self.atomicConstraintObjs.append(atomConstraint)
		self.createTestObjs()
		
		expCastepAtomIdx = 1 #Because its the first "Y" element
		expConstraintNumber = 1 #Simply means first constraint applied
		expStr = "{} {} {} 1.0 0.0 0.0".format(1, self.eleList[atomIdxFixed], expCastepAtomIdx)
		actStr = tCode.getIonicConstraintsStrFromAtomicPosConstraintObj(self.geoConstraints.atomicPostionConstraints, self.eleList)
		self.assertEqual(expStr,actStr)
	
	def testReturnsExpectedStrForMultipleAtomConstraints(self):
		atomIdxAFixed, atomIdxBFixed = 1,2
		atomConstraintA = geomConstraints.AtomicCartesianConstraint(atomIdxAFixed, fixY=True)
		atomConstraintB = geomConstraints.AtomicCartesianConstraint(atomIdxBFixed, fixY=True, fixZ=True)
		self.atomicConstraintObjs.append(atomConstraintA)
		self.atomicConstraintObjs.append(atomConstraintB)
		self.createTestObjs()
		expStr = "1 Y 1 0.0 1.0 0.0\n" + "2 X 2 0.0 1.0 0.0\n" + "3 X 2 0.0 0.0 1.0"
		actStr = tCode.getIonicConstraintsStrFromAtomicPosConstraintObj(self.geoConstraints.atomicPostionConstraints, self.eleList)
		self.assertEqual(expStr,actStr)


class TestGetCastepIndicesFromEleList(unittest.TestCase):

	def setUp(self):
		self.eleList = ["X","Y","X","Z"]

	def testForThreeElementCase(self):
		expIndices = [1,1,2,1] #Indices are labelled witihn atomic species in castep
		actIndices = tCode.getCastepIndicesFromEleList(self.eleList)
		self.assertEqual(expIndices,actIndices)

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

	@mock.patch("gen_basis_helpers.castep.castep_creator.baseObjs.StandardParsedOutputFile")
	@mock.patch("gen_basis_helpers.castep.castep_creator.parseCastep")
	def testParsedFileReturnsExpected(self, mockedParserModule, mockedParsedOutputFileClass):
		expParsedFile, expOutObj = {"a":"a"}, mock.Mock()
		mockedParserModule.parseCastepOutfile.side_effect = lambda *args,**kwargs: expParsedFile
		mockedParsedOutputFileClass.fromKwargDict.side_effect = lambda *args, **kwargs: expOutObj
		actOutObj = self.testObjA.parsedFile
		mockedParsedOutputFileClass.fromKwargDict.assert_called_once_with(**expParsedFile)
		expOutObj.unitCell.convAngToBohr.assert_called_once_with()
		self.assertEqual(expOutObj,actOutObj)



