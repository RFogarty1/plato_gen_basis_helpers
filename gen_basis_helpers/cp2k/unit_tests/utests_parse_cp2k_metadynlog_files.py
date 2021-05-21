
import copy
import itertools as it
import unittest
import unittest.mock as mock

import numpy as np

import gen_basis_helpers.analyse_md.analyse_metadyn_hills as metadynHillsHelp
import gen_basis_helpers.cp2k.parse_md_files as tCode

class TestParseGenericMetadynLogFile(unittest.TestCase):

	def setUp(self):
		self.inpPath = "fake_path_a"
		self.createTestObjs()

	def createTestObjs(self):
		self.fileAsListA = _getGenericMetadynLogFileStringA().split("\n")

	def _runTestFunct(self):
		return tCode._parseGenericMetadynLogFile(self.inpPath)

	@mock.patch("gen_basis_helpers.cp2k.parse_md_files._readFileIntoList")
	def testExpectedResultA(self, mockReadFileToList):
		mockReadFileToList.side_effect = lambda *args,**kwargs: self.fileAsListA
		expVals = [ [70.0, 20.59949, 1.00000, 0.01840],
		            [72.5, 20.95701, 0.90000, 0.01740],
		            [75.0, 21.33221, 1.00000, 0.01840] ]
		actVals = self._runTestFunct()
		mockReadFileToList.assert_called_with(self.inpPath)
		self._checkExpAndActEqual(expVals,actVals)

	def _checkExpAndActEqual(self, expArr, actArr):
		exp = np.array(expArr)
		act = np.array(actArr)
		self.assertTrue( np.allclose(exp,act) )


def _getGenericMetadynLogFileStringA():
	outStr = ""
	outStr += "        70.0     20.59949      1.00000      0.01840\n"
	outStr += "        72.5     20.95701      0.90000      0.01740\n"
	outStr += "        75.0     21.33221      1.00000      0.01840\n"
	return outStr


class TestParseHillsLogFile(unittest.TestCase):

	def setUp(self):
		self.testPathA = "test_path_a"
		self.createTestObjs()

	def createTestObjs(self):
		self.arrA = _loadExpectedHillsArrayTwoColVars()

	def _runTestFunct(self):
		return tCode.parseMetadynamicsHillsLogFile(self.testPathA)

	def _loadExpValsA(self):
		dictA = { "time":[47.5, 50.0], "position":[17.42807, 17.18365],
		          "scale":[1.0,1.0], "height": [0.01840, 0.01840] }
		dictB = { "time":[47.5, 50.0], "position":[12.96243, 13.39474],
		          "scale":[0.5,0.5], "height": [0.01840, 0.01840] }

		return [dictA, dictB]

	@mock.patch("gen_basis_helpers.cp2k.parse_md_files._parseGenericMetadynLogFile")
	def testExpectedTwoColVars(self, mockParseLogFile):
		mockParseLogFile.side_effect = lambda *args,**kwargs: self.arrA
		expVals = self._loadExpValsA()
		actVals = self._runTestFunct()
		mockParseLogFile.assert_called_with(self.testPathA)

		self.assertEqual( len(expVals), len(actVals) )
		for exp,act in zip(expVals, actVals):
			self._checkExpAndActEqual(exp,act)

	def _checkExpAndActEqual(self, expDict, actDict):
		for key in expDict.keys():
			expList,actList = expDict[key],actDict[key]
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expList,actList)]


def _loadExpectedHillsArrayTwoColVars():
	outVals = [ [47.5, 17.42807, 12.96243, 1.00000, 0.50000, 0.01840],
	            [50.0, 17.18365, 13.39474, 1.00000, 0.50000, 0.01840] ]
	return outVals


class TestParseIterOfMetadynHillsFiles(unittest.TestCase):

	def setUp(self):
		self.inpPaths = ["test_path_a", "test_path_b"]
		self.outDictFileAColVarA = { "time":[47.5, 50.0], "position":[17.42807, 17.18365],
		                             "scale":[1.0,1.0], "height": [0.01840, 0.01840] }

		self.outDictFileAColVarB = { "time":[47.5, 50.0], "position":[12.96243, 13.39474],
		                             "scale":[0.5,0.5], "height": [0.01840, 0.01840] }

		self.outDictFileBColVarA = { "time":[52.5, 55.0], "position":[20, 21],
		                             "scale":[1.0,1.0], "height": [0.01840, 0.01840] }

		self.outDictFileBColVarB = { "time":[52.5, 55.0], "position":[22, 23],
		                             "scale":[0.5,0.5], "height": [0.01840, 0.01840] }

		self.createTestObjs()

	def createTestObjs(self):
		self.outputA = [self.outDictFileAColVarA, self.outDictFileAColVarB]
		self.outputB = [self.outDictFileBColVarA, self.outDictFileBColVarB]

	def _runTestFunct(self):
		return tCode.parseIterOfMetadynamicsHillsLogFiles(self.inpPaths)

	def _loadExpOutputA(self):
		sharedTimes = [47.5, 50.0, 52.5, 55.0]
		sharedHeights = [0.01840 for x in range(4)] 
		expColVarA = {"time":sharedTimes, "position": [17.42807, 17.18365, 20, 21],
		              "scale":[1.0 for x in range(4)], "height": sharedHeights}
		expColVarB = {"time":sharedTimes, "position": [12.96243, 13.39474, 22, 23],
		              "scale":[0.5 for x in range(4)], "height": sharedHeights}
		return [expColVarA, expColVarB]

	@mock.patch("gen_basis_helpers.cp2k.parse_md_files.parseMetadynamicsHillsLogFile")
	def testExpCaseA(self, mockParseSingle):
		mockParseSingle.side_effect = [self.outputA, self.outputB]
		expOutput = self._loadExpOutputA()
		actOutput = self._runTestFunct()
		[mockParseSingle.assert_any_call(x) for x in self.inpPaths]

		self.assertEqual( len(expOutput), len(actOutput) )
		for exp,act in zip(expOutput, actOutput):
			self._checkExpAndActEqual(exp,act)

	def _checkExpAndActEqual(self, expDict, actDict):
		for key in expDict.keys():
			expList,actList = expDict[key],actDict[key]
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expList,actList)]


class TestGetGaussianFunctsFromMetadynHills(unittest.TestCase):

	def setUp(self):
		#First collective variable (meaning each is a 1-dimension Gaussian)
		self.heightsA = [2,3,4] 
		self.scalesA = [5,6,7]
		self.positionsA = [8,9,10]
		self.timeA = [1,2,3]

		#Second collective variable
		self.timeB = [1,2,3]
		self.heightsB = [1,2,3]
		self.scalesB = [4,5,6]
		self.positionsB = [7,8,9]

		self.createTestObjs()

	def createTestObjs(self):
		self.colVarDictA = {"time":self.timeA, "height":self.heightsA, "scale":self.scalesA, "position":self.positionsA}
		self.colVarDictB = {"time":self.timeB, "height":self.heightsB, "scale":self.scalesB, "position":self.positionsB}
		self.colVarDicts = [self.colVarDictA, self.colVarDictB]

	def _runTestFunct(self):
		return tCode.getGroupedMultiDimGaussHillsFromParsedHillsFileOutput(self.colVarDicts)

	def _loadExpectedCaseA(self):
		multiDimGaus = list()
		for idx in range(3):
			currMultiDimGau = metadynHillsHelp.MultiDimGaussHill.fromIters( heights=[self.heightsA[idx], self.heightsB[idx]],
		                                                         scales=[self.scalesA[idx], self.scalesB[idx]],
		                                                         positions=[self.positionsA[idx], self.positionsB[idx]])
			multiDimGaus.append(currMultiDimGau) 

		expObj = metadynHillsHelp.GroupedMultiDimGaussHills(multiDimGaus)
		return expObj

	def testExpectedA(self):
		expObj = self._loadExpectedCaseA()
		actObj = self._runTestFunct()
		self.assertEqual(expObj, actObj)



