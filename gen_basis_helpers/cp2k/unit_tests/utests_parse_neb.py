
import copy
import os
import itertools as it
import unittest
import unittest.mock as mock


import plato_pylib.shared.ucell_class as uCellHelp 
import plato_pylib.shared.unit_convs as uConvHelp


import gen_basis_helpers.cp2k.parse_neb_files as tCode



class TestParseFullNebInfo(unittest.TestCase):

	def setUp(self):
		self.convAngToBohr = False
		self.nReplicas = 3
		self.energiesObjs = [mock.Mock(), mock.Mock(), mock.Mock()]
		self.energyVals = [6,7,8]
		self.dists = [1,2,3]
		self.cellA = uCellHelp.UnitCell(lattParams=[10,10,10], lattAngles=[90,90,90]) #Angstrom i think

		self.coordsA = [ [2,3,4,"X"] ]
		self.coordsB = [ [3,4,5,"Y"] ]
		self.coordsC = [ [4,5,6,"Z"] ]

		self.cpoutPath = "fake_out_path.cpout"
		self.createTestObjs()

	def createTestObjs(self):
		summaryVals = {"dists":self.dists,  "energies": self.energyVals}
		self.expSummaryDict = {"unitCell":self.cellA, "final_neb_summary": summaryVals}

	def _runTestFunct(self):
		return tCode.parseNudgedBandCalcStandard(self.cpoutPath, convAngToBohr=self.convAngToBohr)

	@mock.patch("gen_basis_helpers.cp2k.parse_neb_files._parseFinalGeomFromIterOfXyzFiles")
	@mock.patch("gen_basis_helpers.cp2k.parse_neb_files._getXyzFilesFromCpoutPathAndNumbReplicas")
	def testGetNebGeoms(self, mockGetPaths, mockParseFinalGeoms):
		#Setup
		coords = [self.coordsA, self.coordsB, self.coordsC]
		expPaths = ["fake_a","fake_b","fake_c"]
		mockGetPaths.side_effect = lambda *args,**kwargs: expPaths
		mockParseFinalGeoms.side_effect = lambda *args,**kwargs: coords

		#Figure out expected
		expGeoms = [copy.deepcopy(self.cellA) for x in range(len(self.dists))]
		for currCell,currXyz in it.zip_longest(expGeoms,coords):
			currCell.cartCoords = currXyz

		#Run + check accuracy
		actGeoms = tCode._getNebGeomsFromCpoutPathAndNumbReplicasAndEmptyCell(self.cpoutPath, self.nReplicas, self.cellA)
		mockGetPaths.assert_called_with(self.cpoutPath, self.nReplicas)
		mockParseFinalGeoms.assert_called_with(expPaths)
		self.assertEqual(expGeoms, actGeoms)

	@mock.patch("gen_basis_helpers.cp2k.parse_neb_files._parseFinalEnergyObjsFromIterOfOutFiles")
	@mock.patch("gen_basis_helpers.cp2k.parse_neb_files._getBandOutFilesFromCpoutPathAndNumbReplicas")
	def testGetEnergyObjs(self, mockGetPaths, mockParseFinalEnergies):
		#Setup
		expPaths = ["fake_e","fake_f","fake_g"]
		mockGetPaths.side_effect = lambda *args,**kwargs: expPaths
		mockParseFinalEnergies.side_effect = lambda *args,**kwargs: self.energyVals
		expEnergies = self.energyVals

		#Run + check accuracy
		actEnergies = tCode._getEnergiesFromCpoutPathAndNumbReplicas(self.cpoutPath, self.nReplicas)
		mockGetPaths.assert_called_with(self.cpoutPath, self.nReplicas)
		mockParseFinalEnergies.assert_called_with(expPaths)

		self.assertEqual(expEnergies, actEnergies)


class TestGetFilePathsFromCpoutPathAndNImages(unittest.TestCase):

	def setUp(self):
		self.cpoutPath = os.path.join( "fake_dir_a","fake_cpout_a.cpout" )
		self.nReplicas = 4

	def testExpectedBandFiles(self):
		exts = ["-BAND1.out", "-BAND2.out", "-BAND3.out", "-BAND4.out"]
		expPaths = [os.path.join("fake_dir_a","fake_cpout_a") + ext for ext in exts]
		actPaths = tCode._getBandOutFilesFromCpoutPathAndNumbReplicas(self.cpoutPath, self.nReplicas)
		self.assertEqual(expPaths,actPaths)

	def testExpectedXyzFiles(self):
		exts = ["-pos-Replica_nr_1-1.xyz", "-pos-Replica_nr_2-1.xyz",
		        "-pos-Replica_nr_3-1.xyz", "-pos-Replica_nr_4-1.xyz"]
		expPaths = [os.path.join("fake_dir_a","fake_cpout_a") + ext for ext in exts]
		actPaths = tCode._getXyzFilesFromCpoutPathAndNumbReplicas(self.cpoutPath, self.nReplicas)
		self.assertEqual(expPaths, actPaths)


class TestParseSummaryNebCalc(unittest.TestCase):

	def setUp(self):
		self.createTestObjs()

	def createTestObjs(self):
		self.fileAsListA = self._loadFileStrA().split("\n")
		self.startIdxA = 1

	def _loadFileStrA(self):
		return """ *******************************************************************************
 BAND TYPE                     =                                          CI-NEB
 BAND TYPE OPTIMIZATION        =                                            DIIS
 STEP NUMBER                   =                                              56
 NUMBER OF NEB REPLICA         =                                               9
 DISTANCES REP =        0.392508        0.390173        0.384187        0.379559
                        0.380281        0.384784        0.390536        0.392665
 ENERGIES [au] =      -37.564946      -37.564941      -37.564170      -37.558733
                      -37.554902      -37.558751      -37.564184      -37.564944
                      -37.564947
 BAND TOTAL ENERGY [au]        =                             -338.06051542100988
 *******************************************************************************

PROGRAM ENDED AT
"""

	def testExpectedDictFromSummarySectionA(self):
		haToEv = uConvHelp.RYD_TO_EV*2
		expEndIdx = 11
		expDists = [0.392508, 0.390173, 0.384187, 0.379559,
		            0.380281, 0.384784, 0.390536, 0.392665]
		expEnergiesHa = [-37.564946, -37.564941, -37.564170, -37.558733,
		                 -37.554902, -37.558751, -37.564184, -37.564944,
		                 -37.564947]
		expEnergies = [x*haToEv for x in expEnergiesHa]
		
		expDict = {"dists":expDists, "energies":expEnergies,
		           "bandEnergy": haToEv*-338.06051542100988,
		           "bandType": "CI-NEB", "stepNumber":56}

		actDict, actEndIdx = tCode._parseNebSummarySection(self.fileAsListA, self.startIdxA)
		self.assertEqual(expEndIdx,actEndIdx)
		self._checkExpDictMatchesAct(expDict, actDict)

	@mock.patch("gen_basis_helpers.cp2k.parse_neb_files._parseNebSummarySection")
	@mock.patch("gen_basis_helpers.cp2k.parse_neb_files._getFileAsListFromInpFile")
	def testExpectedDictFromFullParseFunction(self, mockedGetFileAsList, mockedParseNebSummary):
		expNebSummary = mock.Mock()
		mockedGetFileAsList.side_effect = lambda *args,**kwargs : self.fileAsListA


		def _mockedParser(fileAsList, lineIdx):
			return expNebSummary,lineIdx+1
		mockedParseNebSummary.side_effect = _mockedParser

		actDict = tCode.parseNebSummaryFile( mock.Mock() )
		actNebSummary = actDict["final_neb_summary"] 
		self.assertEqual(expNebSummary, actNebSummary)

	def _checkExpDictMatchesAct(self, expDict, actDict):
		floatIterKeys = ["dists", "energies"]
		otherKeys = ["bandType", "stepNumber", "bandEnergy"]

		for key in otherKeys:
			expVal, actVal = expDict[key], actDict[key]
			self.assertAlmostEqual(expVal,actVal)

		for key in floatIterKeys:
			expVals, actVals = expDict[key], actDict[key]
			[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expVals,actVals)]



