
import copy
import os
import unittest
import unittest.mock as mock

import plato_pylib.shared.unit_convs as uConvHelp

import gen_basis_helpers.cp2k.parse_pdos_files as tCode




class TestParsePdosFromCpout(unittest.TestCase):

	def setUp(self):
		self.cpoutPath = "fake_cpout"
		self.parsedAtomKinds = [mock.Mock(), mock.Mock()]
		self.parsedAtomLists = [mock.Mock(), mock.Mock(), mock.Mock()]

		self.createTestObjs()

	def createTestObjs(self):
		pass

	@mock.patch("gen_basis_helpers.cp2k.parse_pdos_files._parsePdosListFilesFromCpoutPath")
	@mock.patch("gen_basis_helpers.cp2k.parse_pdos_files._parsePdosKindFilesFromCpoutPath")
	def testExpectedCallsToParseFileTypes(self, mockedParseKindFiles, mockedParseListFiles):
		mockedParseKindFiles.side_effect = lambda *args,**kwargs: self.parsedAtomKinds
		mockedParseListFiles.side_effect = lambda *args,**kwargs: self.parsedAtomLists

		expOutDict = {"atomKinds":self.parsedAtomKinds, "atomLists":self.parsedAtomLists}
		actOutDict = tCode.parsePdosFromCpoutPath(self.cpoutPath)

		mockedParseKindFiles.assert_called_with( self.cpoutPath )
		mockedParseListFiles.assert_called_with( self.cpoutPath )
		self.assertEqual(expOutDict,actOutDict)

	@mock.patch("gen_basis_helpers.cp2k.parse_pdos_files.parsePdosFromFile")
	@mock.patch("gen_basis_helpers.cp2k.parse_pdos_files.getPdosKindsPathsFromCpoutPath")
	def testGetExpectedPdosCallsToParseAtomKinds(self, mockGetPaths, mockParsePdos):
		#Setup
		expPaths = ["fake_path_a", "fake_path_b"]
		expPdos = expPaths
		mockGetPaths.side_effect  = lambda *args,**kwargs: expPaths
		mockParsePdos.side_effect = lambda inpPath: inpPath 

		#
		actPdos = tCode._parsePdosKindFilesFromCpoutPath(self.cpoutPath)
		mockGetPaths.assert_called_with(self.cpoutPath)
		self.assertEqual(expPdos, actPdos)

	@mock.patch("gen_basis_helpers.cp2k.parse_pdos_files.parsePdosFromFile")
	@mock.patch("gen_basis_helpers.cp2k.parse_pdos_files.getPdosAtomicListsPathsFromCpoutPath")
	def testGetExpectedPdosCallsToParseListKinds(self, mockGetPaths, mockParsePdos):
		#Setup
		expPaths = ["fake_path_a", "fake_path_b"]
		expPdos = expPaths
		mockGetPaths.side_effect  = lambda *args,**kwargs: expPaths
		mockParsePdos.side_effect = lambda inpPath: inpPath 

		#
		actPdos = tCode._parsePdosListFilesFromCpoutPath(self.cpoutPath)
		mockGetPaths.assert_called_with(self.cpoutPath)
		self.assertEqual(expPdos, actPdos)



class TestGetFilePaths(unittest.TestCase):

	def setUp(self):
		self.fakeCpoutFolder = "fake_folder"
		self.fakeCpoutName = "geom_opt"
		self.extensions = ["-k1-1.pdos", ".cpout" , "-k2-1.pdos", "-list1-1.pdos",
		                   "-list2-1.pdos", "-list3-1.pdos"]
		self.createTestObjs()

	def createTestObjs(self):
		self.mockFilesPresent = [ self.fakeCpoutName+x for x in self.extensions ] 
		self.cpoutPath = os.path.join(self.fakeCpoutFolder, self.fakeCpoutName)

	@mock.patch("gen_basis_helpers.cp2k.parse_pdos_files.os.listdir")
	def testGetKindBreakdownPathsA(self, mockedListDir):
		currExts = ["-k1-1.pdos", "-k2-1.pdos"]
		mockedListDir.side_effect = lambda *args,**kwargs: self.mockFilesPresent
		expPaths = [os.path.join(self.fakeCpoutFolder, self.fakeCpoutName+x) for x in currExts]
		actPaths = tCode.getPdosKindsPathsFromCpoutPath(self.cpoutPath)
		mockedListDir.assert_called_with(self.fakeCpoutFolder)
		self.assertEqual(expPaths, actPaths)

	@mock.patch("gen_basis_helpers.cp2k.parse_pdos_files.os.listdir")
	def testGetListBreakdownPathsA(self, mockedListDir):
		currExts = ["-list1-1.pdos", "-list2-1.pdos", "-list3-1.pdos"]
		mockedListDir.side_effect = lambda *args,**kwargs: self.mockFilesPresent
		expPaths = [os.path.join(self.fakeCpoutFolder, self.fakeCpoutName+x) for x in currExts]
		actPaths = tCode.getPdosAtomicListsPathsFromCpoutPath(self.cpoutPath)
		mockedListDir.assert_called_with(self.fakeCpoutFolder)
		self.assertEqual(expPaths, actPaths)
 

class TestParseSinglePdos(unittest.TestCase):

	def setUp(self):
		self.filePath = "fake_path_a"
		self.createTestObjs()

	def createTestObjs(self):
		self.fileAsListA = self._loadFileStrA().split("\n")
		self.expObjA = self._loadExpObjA()

	def _loadFileStrA(self):
		return """# Projected DOS for atomic kind O at iteration step i = 0, E(Fermi) =    -0.260014 a.u.
#     MO Eigenvalue [a.u.]      Occupation                 s                 p                 d
       1         -0.921635        2.000000        0.45700077        0.01994362        0.00076722
       2          0.007051        0.000000        0.15757794        0.08100025        0.02367271"""


	def _loadExpObjA(self):
		haToEv = uConvHelp.RYD_TO_EV*2
		eigenValues = [x*haToEv for x in [-0.921635, 0.007051]]
		occs = [2,0]
		breakdowns = [ [0.45700077, 0.01994362, 0.00076722],
		               [0.15757794, 0.08100025, 0.02367271] ]
		fragName = "O"
		breakdownHeaders = ["s", "p", "d"]
		currKwargs = {"eigenValues":eigenValues, "occs":occs, "fragName":fragName,
		              "breakdowns":breakdowns, "breakdownHeaders":breakdownHeaders}
		return tCode.PdosFragmentStandard(**currKwargs)

	def _runTestFunct(self):
		return tCode.parsePdosFromFile(self.filePath)

	@mock.patch("gen_basis_helpers.cp2k.parse_pdos_files._readFileIntoList")
	def testExpectedValsA(self, mockReadFileIntoList):
		mockReadFileIntoList.side_effect = lambda *args,**kwargs: self.fileAsListA
		actObj = self._runTestFunct()
		mockReadFileIntoList.assert_called_with(self.filePath)
		self.assertEqual(self.expObjA, actObj)


class TestCp2kPdosFragmentClass(unittest.TestCase):

	def setUp(self):
		self.occs      = [2   , 2   , 1.4, 0.8]
		self.eigenVals = [-3.4, -1.4, 0.5, 0.9]
		self.fragName = "testFrag"
		self.breakdowns = [ [ 1, 2],
		                    [ 5, 6],
		                    [ 3, 4],
		                    [ 6, 7] ]

		self.breakdownHeaders = ["a", "b"]
		self.createTestObjs()

	def createTestObjs(self):
		currKwargDict = {"eigenValues":self.eigenVals, "occs":self.occs,
		                 "fragName":self.fragName, "breakdowns":self.breakdowns,
		                 "breakdownHeaders": self.breakdownHeaders}
		self.testObj = tCode.PdosFragmentStandard(**currKwargDict)

	def testTwoEqualObjsCompareEqual(self):
		objA = copy.deepcopy(self.testObj)
		self.createTestObjs()
		self.assertEqual( objA, self.testObj )
		
	def testTwoUnequalCompareUnequal_diffHeaders(self):
		objA = copy.deepcopy(self.testObj)
		self.breakdownHeaders[-1] = "d"
		self.createTestObjs()
		self.assertNotEqual( objA, self.testObj )

	def testTwoUnequalCompareUnequal_diffEigenvals(self):
		objA = copy.deepcopy(self.testObj)
		self.eigenVals[-1] += 2
		self.createTestObjs()
		self.assertNotEqual( objA, self.testObj )

	def testTwoUnequalCompareUnequal_diffBreakdowns(self):
		objA = copy.deepcopy(self.testObj)
		self.breakdowns[-1][0] += 3
		self.createTestObjs()
		self.assertNotEqual( objA, self.testObj )

	def testToAndFromDictConsistent(self):
		dictA = self.testObj.toDict()
		newObj = tCode.PdosFragmentStandard.fromDict(dictA)
		self.assertEqual(self.testObj, newObj)


