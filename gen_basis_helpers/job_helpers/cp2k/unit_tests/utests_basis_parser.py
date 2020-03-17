
import copy
import itertools as it
import os

import unittest
import unittest.mock as mock

import gen_basis_helpers.job_helpers.cp2k.parse_basis_help as tCode

import plato_pylib.parseOther.parse_cp2k_basis as parseCP2KBasis
import gen_basis_helpers.shared.gaussians as gauHelp



class TestBasisParser(unittest.TestCase):

	def setUp(self):
		self.basisSetFolder = "fake_folder"
		self.basisSetFile = "fake_file"
		self.element = "Mg"
		self.basisNameInFile = "fake-basis-name"
		self.basisAlias = "fake-alias"
		self.coeffsForUnormalisedGaussians = True
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.BasisParser(basisSetFolder=self.basisSetFolder, basisSetFile=self.basisSetFile,
		                                  basisNameInFile=self.basisNameInFile, element=self.element,
		                                  basisAlias=self.basisAlias,
		                                  coeffsForUnormalisedGaussians=self.coeffsForUnormalisedGaussians)

	def testBasisPathAsExpected(self):
		expPath = os.path.join(self.basisSetFolder, self.basisSetFile)
		actPath = self.testObjA._basisPath
		self.assertEqual(expPath,actPath)

	@mock.patch("gen_basis_helpers.job_helpers.cp2k.parse_basis_help.parseCP2KBasis.parseCP2KBasisFile")
	def testGetBasisInNativeFormat(self, mockedParseFunct):
		mockedParsedBasis = mock.Mock()
		expSpecificBasisSet = mock.Mock()
		mockedParseFunct.side_effect = lambda x: mockedParsedBasis
		mockedParsedBasis.getUniqueBasisSet.side_effect = lambda *args: expSpecificBasisSet
		expPath = os.path.join(self.basisSetFolder, self.basisSetFile)
		
		actSpecificBasisSet = self.testObjA._getBasisInNativeFormat()
		mockedParseFunct.assert_called_once_with(expPath)
		mockedParsedBasis.getUniqueBasisSet.assert_called_once_with(self.element, self.basisNameInFile)
		self.assertEqual(expSpecificBasisSet, actSpecificBasisSet)

	@mock.patch("gen_basis_helpers.job_helpers.cp2k.parse_basis_help.getBasisFunctObjsFromCP2KBasisSetNativeFormat")
	@mock.patch("gen_basis_helpers.job_helpers.cp2k.parse_basis_help.BasisParser._getBasisInNativeFormat")
	@mock.patch("gen_basis_helpers.job_helpers.cp2k.parse_basis_help.basObjs.OrbitalBasisSetStandard")
	def testExpectedCallsWhenParsing(self, mockedOrbBasSetClass, mockedGetNativeFmt, mockedGetBasisFuncts):
		expBasisFuncts, expNativeFmt, expOutObj = mock.Mock(), mock.Mock(), mock.Mock()
		mockedGetNativeFmt.side_effect = lambda *args,**kwargs: expNativeFmt
		mockedGetBasisFuncts.side_effect = lambda *args,**kwargs: expBasisFuncts
		mockedOrbBasSetClass.side_effect = lambda *args,**kwargs: expOutObj
		
		outObj = self.testObjA.parse()

		mockedGetNativeFmt.assert_called_once_with() #Probably a superfluos test
		mockedGetBasisFuncts.assert_called_once_with(expNativeFmt,convNormToRawGaussians=self.coeffsForUnormalisedGaussians)
		mockedOrbBasSetClass.assert_called_once_with(expBasisFuncts,self.basisAlias)
		self.assertEqual(expOutObj,outObj)

	def testExpectedAliasWhenSetAsNone(self):
		self.basisAlias = None
		self.createTestObjs()
		expAlias = self.basisNameInFile
		actAlias = self.testObjA.basisAlias
		self.assertEqual(expAlias,actAlias)


class TestSingleExponentSetToCP2KBasis(unittest.TestCase):

	def setUp(self):
		self.exponentsA = [1,2]
		self.coeffsA = [ [1,2], [3,4] ]
		self.lValsA = [0,1]
		self.nValA = 3 #Irrelevant except for creating the test object
		self.convNormToRawGaussians = False
		self.createTestObjs()

	def createTestObjs(self):
		self.testExponentSetA = parseCP2KBasis.ExponentSetCP2K( self.exponentsA, self.coeffsA, self.lValsA, self.nValA)

	def _runTestFunctOnExpSetA(self):
		return tCode._getBasisFunctObjFromExponentSet(self.testExponentSetA, convNormToRawGaussians=self.convNormToRawGaussians)

	@mock.patch("gen_basis_helpers.job_helpers.cp2k.parse_basis_help.basObjs.GauSumOrbitalBasisFunction")
	def testExpCallsToCreateBasisFuncts(self, mockedBasisFunctClass):
		expNValAll = 1 #One of each l-value in this set

		#Manually define the gaussians we expect
		gauPrimAA = gauHelp.GauPrim.fromExpAndCoeffOnly(self.exponentsA[0], self.coeffsA[0][0])
		gauPrimAB = gauHelp.GauPrim.fromExpAndCoeffOnly(self.exponentsA[1], self.coeffsA[0][1])
		gauCompA = gauHelp.GauPrimComposite([gauPrimAA,gauPrimAB])

		gauPrimBA = gauHelp.GauPrim.fromExpAndCoeffOnly(self.exponentsA[0], self.coeffsA[1][0])
		gauPrimBB = gauHelp.GauPrim.fromExpAndCoeffOnly(self.exponentsA[1], self.coeffsA[1][1])
		gauCompB = gauHelp.GauPrimComposite([gauPrimBA,gauPrimBB])

		self._runTestFunctOnExpSetA()

		expCallArgs = [(expNValAll, self.lValsA[0], gauCompA), (expNValAll, self.lValsA[1], gauCompB)]
		actCallArgs = [args for args,kwargs in mockedBasisFunctClass.call_args_list]

		self.assertEqual(expCallArgs,actCallArgs)

	@mock.patch("gen_basis_helpers.job_helpers.cp2k.parse_basis_help.calcNormConstantForCP2KOnePrimitive")
	def testExpectedCallsToCalcNormConstants(self, mockedNormFunct):
		callAA = [ self.exponentsA[0], self.lValsA[0] ]
		callAB = [ self.exponentsA[1], self.lValsA[0] ]
		callBA = [ self.exponentsA[0], self.lValsA[1] ]
		callBB = [ self.exponentsA[1], self.lValsA[1] ]

		expCallArgs = [ callAA, callAB, callBA, callBB ]

		self.convNormToRawGaussians = True
		self._runTestFunctOnExpSetA()

		for x in expCallArgs:
			mockedNormFunct.assert_any_call(*x)


class TestCP2KBasisSetToBasisFuncts(unittest.TestCase):

	def setUp(self):
		self.expSetA = mock.Mock()
		self.expSetB = mock.Mock()
		self.convNormToRawGaussians = False
		self.createTestObjs()

	def createTestObjs(self):
		self.expSets = [self.expSetA, self.expSetB]
		self.basisObjA = mock.Mock()
		self.basisObjA.exponentSets = self.expSets

	def _runTestFunctOnObjA(self):
		return tCode.getBasisFunctObjsFromCP2KBasisSetNativeFormat(self.basisObjA, convNormToRawGaussians=self.convNormToRawGaussians)

	@mock.patch("gen_basis_helpers.job_helpers.cp2k.parse_basis_help._getBasisFunctObjFromExponentSet")
	def testSingleSetConverterCalledForEach(self, mockedSingleSetConverter):
		self._runTestFunctOnObjA()
		for x in self.expSets:
			mockedSingleSetConverter.assert_any_call(x,convNormToRawGaussians=self.convNormToRawGaussians)


	@mock.patch("gen_basis_helpers.job_helpers.cp2k.parse_basis_help._getBasisFunctObjFromExponentSet")
	def testExpectedNValsGiven(self, mockedSingleSetConverter):
		testLValues = [0,1,1] #This will get passed back for EACH call
		mockOutput = [mock.Mock() for x in testLValues] 
		for lVal,mockObj in zip(testLValues, mockOutput):
			mockObj.lVal = lVal
			mockObj.nVal = 1

		mockedSingleSetConverter.side_effect = lambda *args, **kwargs: copy.deepcopy(mockOutput)

		expNVals = [1,1,2] + [2,3,4] #The same testLValues are returned twice; hecne the two lists
		actOutput = self._runTestFunctOnObjA()
		actNVals = [x.nVal for x in actOutput]
		self.assertEqual(expNVals,actNVals)


