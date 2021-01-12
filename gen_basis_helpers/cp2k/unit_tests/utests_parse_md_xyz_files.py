
import copy
import unittest
import unittest.mock as mock

import gen_basis_helpers.cp2k.parse_md_files as tCode


class TestParseCP2KMDXyz(unittest.TestCase):

	def setUp(self):
		self.createTestObjs()

	def createTestObjs(self):
		self.fileStrA = _loadFileStrA().split("\n")
		self.coordsStepA = [ [ 0.0067115483, 0.0018448598, 0.0036115223, "Mg"], 
		                     [-0.0067115636, 1.8434995966, 2.5814191416, "Mg"] ]
		self.coordsStepB = [ [ 0.0083879577, 0.0023044814, 0.0045100468, "Mg"],
		                     [-0.0083879729, 1.8430399749, 2.5805206170, "Mg"] ]
		self.coordsStepC = [ [ 0.0100632990, 0.0027633584, 0.0054058351, "Mg"],
		                     [-0.0100633140, 1.8425810978, 2.5796248286, "Mg"] ]

		self.expOutputA = [ {"coords":self.coordsStepA, "step":4, "time":2.000},
		                  {"coords":self.coordsStepB, "step":5, "time":2.500},
		                  {"coords":self.coordsStepC, "step":6, "time":3.000} ]


	@mock.patch("gen_basis_helpers.cp2k.parse_md_files._readFileIntoList")
	def testExpectedResultsForTestFileA(self, mockedReadFileIntoList):
		fakeFilePath = "fake_file_path"
		mockedReadFileIntoList.side_effect = lambda *args, **kwargs: self.fileStrA
		actOutputA = tCode.parseCp2kMdXyzFile(fakeFilePath)
		mockedReadFileIntoList.assert_called_with(fakeFilePath)

		self.assertEqual( len(self.expOutputA), len(actOutputA) )
		for exp,act in zip(self.expOutputA, actOutputA):
			self._checkExpAndActEntryEqual(exp,act)

	@mock.patch("gen_basis_helpers.cp2k.parse_md_files._readFileIntoList")
	def testExpectedResultsIfJobRunTwice(self, mockedReadFileIntoList):
		""" Want to check we only get parsed results from the LAST run of the job; we can usually tell if one jobs appended to another by adjacent step numbers decreasing """
		self.fileStrA = copy.deepcopy(self.fileStrA) + copy.deepcopy(self.fileStrA)
		mockedReadFileIntoList.side_effect = lambda *args, **kwargs: self.fileStrA
		actOutputA = tCode.parseCp2kMdXyzFile("fake_file_path")

		self.assertEqual( len(self.expOutputA), len(actOutputA) )
		for exp,act in zip(self.expOutputA, actOutputA):
			self._checkExpAndActEntryEqual(exp,act)


	def _checkExpAndActEntryEqual(self, expDict, actDict):

		self.assertEqual( len(expDict["coords"]), len(actDict["coords"]) )
		for expCoords,actCoords in zip( expDict["coords"], actDict["coords"] ):
			[self.assertAlmostEqual(e,a) for e,a in zip(expCoords[:3], actCoords[:3])]
			self.assertEqual(expCoords[-1],actCoords[-1])

		self.assertEqual(expDict["step"],actDict["step"])
		self.assertAlmostEqual( expDict["time"], actDict["time"] )

def _loadFileStrA():
	outStr = """       2
 i =        4, time =        2.000, E =        -1.7610217358
 Mg         0.0067115483        0.0018448598        0.0036115223
 Mg        -0.0067115636        1.8434995966        2.5814191416
       2
 i =        5, time =        2.500, E =        -1.7610194305
 Mg         0.0083879577        0.0023044814        0.0045100468
 Mg        -0.0083879729        1.8430399749        2.5805206170
       2
 i =        6, time =        3.000, E =        -1.7610164892
 Mg         0.0100632990        0.0027633584        0.0054058351
 Mg        -0.0100633140        1.8425810978        2.5796248286"""
	return outStr





