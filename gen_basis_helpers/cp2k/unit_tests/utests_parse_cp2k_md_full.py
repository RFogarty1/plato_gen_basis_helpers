
import copy
import itertools as it
import unittest
import unittest.mock as mock
import types

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.traj_core as trajHelp
import gen_basis_helpers.analyse_md.thermo_data as thermoHelp
import gen_basis_helpers.cp2k.parse_md_files as tCode



class TestParseMultipleCP2kFull(unittest.TestCase):

	def setUp(self):
		self.stepsA = [1,2]
		self.stepsB = [3,4]
		self.testXyzPaths = ["fake_xyz_path_a","fake_xyz_path_b"]
		self.testCpoutPaths = ["fake_cpout_path_a", "fake_cpout_path_b"]
		self.createTestObjs()

	def createTestObjs(self):
		self.thermoDataA = thermoHelp.ThermoDataStandard( {"step":self.stepsA} )
		self.thermoDataB = thermoHelp.ThermoDataStandard( {"step":self.stepsB} )
		self.trajA = trajHelp.TrajectoryInMemory( [trajHelp.TrajStepBase(step=x) for x in self.stepsA] )
		self.trajB = trajHelp.TrajectoryInMemory( [trajHelp.TrajStepBase(step=x) for x in self.stepsB] )
		self.retDictA = {"trajectory":self.trajA, "thermo_data":self.thermoDataA}
		self.retDictB = {"trajectory":self.trajB, "thermo_data":self.thermoDataB}

		self.expThermoData = thermoHelp.ThermoDataStandard( {"step":self.stepsA + self.stepsB} )
		self.expTraj = trajHelp.TrajectoryInMemory(  [trajHelp.TrajStepBase(step=x) for x in self.stepsA+self.stepsB] )
		self.expDictA = {"trajectory":self.expTraj, "thermo_data":self.expThermoData}

	@mock.patch("gen_basis_helpers.cp2k.parse_md_files.parseFullMdInfoFromCpoutAndXyzFilePaths")
	def testExpectedResultGivenA(self, mockedParseCpoutAndXyz):
		def fake_parser(inpCpout, inpXyz, **kwargs):
			if (inpCpout == self.testCpoutPaths[0]) and (inpXyz == self.testXyzPaths[0]):
				return self.retDictA
			elif (inpCpout == self.testCpoutPaths[1]) and (inpXyz == self.testXyzPaths[1]):
				return self.retDictB
			else:
				return None


		mockedParseCpoutAndXyz.side_effect = fake_parser
		actResult = tCode.parseMdInfoFromMultipleCpoutAndXyzPaths(self.testCpoutPaths, self.testXyzPaths)

		mockedParseCpoutAndXyz.assert_any_call(self.testCpoutPaths[0], self.testXyzPaths[0], tempKindPath=None, velocityPath=None, forcePath=None)
		mockedParseCpoutAndXyz.assert_any_call(self.testCpoutPaths[1], self.testXyzPaths[1], tempKindPath=None, velocityPath=None, forcePath=None)

		for key in self.expDictA:
			self.assertEqual( self.expDictA[key], actResult[key] )	

	@mock.patch("gen_basis_helpers.cp2k.parse_md_files.parseFullMdInfoFromCpoutAndXyzFilePaths")
	def testExpectedWhenOrderSwitched(self, mockedParseCpoutAndXyz):
		def fake_parser(inpCpout, inpXyz, **kwargs):
			if (inpCpout == self.testCpoutPaths[0]) and (inpXyz == self.testXyzPaths[0]):
				return self.retDictA
			elif (inpCpout == self.testCpoutPaths[1]) and (inpXyz == self.testXyzPaths[1]):
				return self.retDictB
			else:
				return None

		mockedParseCpoutAndXyz.side_effect = fake_parser
		self.testCpoutPaths = [x for x in reversed(self.testCpoutPaths)]
		self.testXyzPaths = [x for x in reversed(self.testXyzPaths)]
		actResult = tCode.parseMdInfoFromMultipleCpoutAndXyzPaths(self.testCpoutPaths, self.testXyzPaths)

		mockedParseCpoutAndXyz.assert_any_call(self.testCpoutPaths[0], self.testXyzPaths[0], tempKindPath=None, velocityPath=None, forcePath=None)
		mockedParseCpoutAndXyz.assert_any_call(self.testCpoutPaths[1], self.testXyzPaths[1], tempKindPath=None, velocityPath=None, forcePath=None)

		for key in self.expDictA:
			self.assertEqual( self.expDictA[key], actResult[key] )	


#See also TestGetMergedTrajectoryCpoutAndXyz
class TestParseCp2kMdFull(unittest.TestCase):

	def setUp(self):
		self.emptyCellA = uCellHelp.UnitCell(lattParams=[2,2,2], lattAngles=[90,90,90])
		self.emptyCellB = uCellHelp.UnitCell(lattParams=[4,4,4], lattAngles=[90,90,90])
		self.coordsA = [ [1,1,1,"Mg"] ]
		self.coordsB = [ [2,2,2,"Mg"] ]
		self.steps = [1,2]
		self.times = [2.5,5.0]
		self.eKinetic = [4,5]

		self.initMdCell = uCellHelp.UnitCell(lattParams=[3,3,3], lattAngles=[90,90,90])
		self.initThermoDict = {"step":0, "time":0, "eKinetic":1}
		self.initCoords = [ [0,0,0,"Mg"] ]

		self.atomicKindDictA = {"step":self.steps, "time":self.times, "kindTemp":[[200,250]]}

		self.createTestObjs()

	def createTestObjs(self):
		#Create thermo arrays + thermoDataObj
		self.thermalArraysA = {"step":self.steps, "time":self.times, "eKinetic":self.eKinetic}
		self.thermoDataA = thermoHelp.ThermoDataStandard(self.thermalArraysA)

		#thngs to help mocking
		trajSteps = [ trajHelp.TrajStepBase(unitCell=self.emptyCellA, step=self.steps[0], time=self.times[0]),
		              trajHelp.TrajStepBase(unitCell=self.emptyCellB, step=self.steps[1], time=self.times[1]) ]
		self.outTrajA = trajHelp.TrajectoryInMemory(trajSteps)
		self.parsedCpout = {"trajectory":self.outTrajA, "thermo_data":self.thermoDataA,
		                    "init_thermo_dict":self.initThermoDict, "init_md_cell":self.initMdCell}
		self.parsedXyz = [ {"coords":self.coordsA, "step":self.steps[0], "time":self.times[0]},
		                   {"coords":self.coordsB, "step":self.steps[1], "time":self.times[1]} ]

		#Expected output
		self.expCellA, self.expCellB = copy.deepcopy(self.emptyCellA), copy.deepcopy(self.emptyCellB)
		self.expCellA.cartCoords = self.coordsA
		self.expCellB.cartCoords = self.coordsB

		expTrajSteps = copy.deepcopy(trajSteps)
		expTrajSteps[0].unitCell.cartCoords = self.coordsA
		expTrajSteps[1].unitCell.cartCoords = self.coordsB
		self.expDict = dict()
		self.expDict["trajectory"] = trajHelp.TrajectoryInMemory(expTrajSteps)
		self.expDict["thermo_data"] = copy.deepcopy(self.thermoDataA)

	@mock.patch("gen_basis_helpers.cp2k.parse_md_files.parseCp2kMdXyzFile")
	@mock.patch("gen_basis_helpers.cp2k.parse_md_files.parseCpoutForMDJob")
	def testExpValsA(self, mockedParseCpout, mockedParseXyz):
		#Set mocks
		expCpoutPath, expXyzPath = "fake_cpout_path", "fake_xyz_path"
		mockedParseCpout.side_effect = lambda *args,**kwargs: self.parsedCpout
		mockedParseXyz.side_effect = lambda *args,**kwargs: self.parsedXyz

		#Run + test
		actDict = tCode.parseFullMdInfoFromCpoutAndXyzFilePaths(expCpoutPath, expXyzPath)
		mockedParseCpout.assert_called_with(expCpoutPath)
		mockedParseXyz.assert_called_with(expXyzPath)

		self.assertEqual(self.expDict, actDict)


	@mock.patch("gen_basis_helpers.cp2k.parse_md_files.parseCp2kMdXyzFile")
	@mock.patch("gen_basis_helpers.cp2k.parse_md_files.parseCpoutForMDJob")
	def testExpVals_step0InXyz(self, mockedParseCpout, mockedParseXyz):
		#Modify the xyz output
		stepZeroXyz = {"coords":self.initCoords, "step":0, "time":0}
		self.parsedXyz.insert(0, stepZeroXyz)

		#modify the expected dict
		trajStepZero = trajHelp.TrajStepBase(unitCell=copy.deepcopy(self.initMdCell), step=0, time=0)
		self.expDict["trajectory"].trajSteps = [trajStepZero] + self.expDict["trajectory"].trajSteps
		self.expDict["trajectory"].trajSteps[0].unitCell.cartCoords = self.initCoords
		for key in self.expDict["thermo_data"].dataDict.keys():
			currVal = self.initThermoDict[key]
			self.expDict["thermo_data"].dataDict[key].insert(0, currVal)

		#Set mocks
		expCpoutPath, expXyzPath = "fake_cpout_path", "fake_xyz_path"
		mockedParseCpout.side_effect = lambda *args,**kwargs: self.parsedCpout
		mockedParseXyz.side_effect = lambda *args,**kwargs: self.parsedXyz

		#Run + test
		actDict = tCode.parseFullMdInfoFromCpoutAndXyzFilePaths(expCpoutPath, expXyzPath)
		mockedParseCpout.assert_called_with(expCpoutPath)
		mockedParseXyz.assert_called_with(expXyzPath)

		self.assertEqual(self.expDict, actDict)

	@unittest.skip("Not sure theres any reason to raise assert error here; since the trajectory/xyz will only be merged when step numbers match anyway")
	@mock.patch("gen_basis_helpers.cp2k.parse_md_files.parseCp2kMdXyzFile")
	@mock.patch("gen_basis_helpers.cp2k.parse_md_files.parseCpoutForMDJob")
	def testRaisesWhenNumberOfStepsDiffer(self, mockedParseCpout, mockedParseXyz):
		#Set mocks
		expCpoutPath, expXyzPath = "fake_cpout_path", "fake_xyz_path"
		mockedParseCpout.side_effect = lambda *args,**kwargs: self.parsedCpout
		mockedParseXyz.side_effect = lambda *args,**kwargs: self.parsedXyz
#		self.parsedXyz.append( copy.deepcopy(self.parsedXyz[0]) )
#		self.parsedCpout.append( copy.deepcopy(self.parsedCpout[0]) )
		self.parsedCpout["trajectory"].trajSteps.append( copy.deepcopy(self.parsedCpout["trajectory"].trajSteps[0]) )
		with self.assertRaises(AssertionError):
			tCode.parseFullMdInfoFromCpoutAndXyzFilePaths(expCpoutPath, expXyzPath)

	@mock.patch("gen_basis_helpers.cp2k.parse_md_files.parseAtomTempFile")
	@mock.patch("gen_basis_helpers.cp2k.parse_md_files.parseCp2kMdXyzFile")
	@mock.patch("gen_basis_helpers.cp2k.parse_md_files.parseCpoutForMDJob")
	def testExpVals_tempKindPathSet(self, mockedParseCpout, mockedParseXyz, mockedParseTempKind):
		#Set mocks
		expCpoutPath, expXyzPath, expAtomTempPath = "fake_cpout_path", "fake_xyz_path", "fake_atomic_temp_path"
		mockedParseCpout.side_effect = lambda *args,**kwargs: self.parsedCpout
		mockedParseXyz.side_effect = lambda *args,**kwargs: self.parsedXyz
		mockedParseTempKind.side_effect = lambda *args,**kwargs: self.atomicKindDictA

		#Mod expected dict to include atomic kind temperatures
		expDict = copy.deepcopy(self.expDict)
		expDict["thermo_data"].dataDict["kindTemp_Mg"] = [200,250]

		#Run and test
		actDict = tCode.parseFullMdInfoFromCpoutAndXyzFilePaths(expCpoutPath, expXyzPath, tempKindPath=expAtomTempPath)
		mockedParseCpout.assert_called_with(expCpoutPath)
		mockedParseXyz.assert_called_with(expXyzPath)
		mockedParseTempKind.assert_called_with(expAtomTempPath)
		self.assertEqual(expDict, actDict)

class TestGetMergedTrajectoryCpoutAndXyz(unittest.TestCase):

	def setUp(self):
		self.stepsCpout = [1, 2  , 3  , 4  ,5]
		self.timesCpout = [2, 2.5, 4  , 5  ,6]
		self.stepsXyz = [1,3,5]
		self.timesXyz = [2,4,6]
		self.createTestObjs()

	def createTestObjs(self):
		self.cartXyz   = [mock.Mock() for x in range(len(self.stepsXyz))]
		self.cellCpout = [mock.Mock() for x in range(len(self.stepsCpout))]
		self.cpoutTSteps = [types.SimpleNamespace(step=s, time=t, unitCell=g) for s,t,g in zip(self.stepsCpout, self.timesCpout, self.cellCpout)]
		self.parsedCpoutA = { "trajectory": trajHelp.TrajectoryInMemory(self.cpoutTSteps) }
		self.parsedXyz = [ {"step":s, "time":t, "coords":g} for s,t,g in zip(self.stepsXyz, self.timesXyz, self.cartXyz) ]

	def testExpected_parsedXyzHasLessSteps(self):
		outObj = tCode._getMergedTrajectoryFromParsedCpoutAndXyz(self.parsedCpoutA, self.parsedXyz)

		self.assertEqual( len(outObj.trajSteps), len(self.stepsXyz) )
		expModified = [self.cpoutTSteps[idx] for idx in [0,2,4]]
		for trajStep, cart in it.zip_longest(expModified, self.cartXyz):
			self.assertEqual(trajStep.unitCell.cartCoords,cart)


class TestGetKindIdxToSymbolDict(unittest.TestCase):

	def setUp(self):
		self.eleList = ["Mg", "Mg","X","Y","X"]

	def _runTestFunct(self):
		return tCode._getKindIdxToSymbolDict(self.eleList)

	def testExpectedValsA(self):
		expDict = {0: "Mg", 1: "X", 2: "Y"}
		actDict = self._runTestFunct()
		self.assertEqual(expDict,actDict)


#Probably mock out a function which parses the file like it was an xyz co-ord file
class TestParseVelocitiesAndAddToTrajFile(unittest.TestCase):

	def setUp(self):
		self.fakePathA = "fake_path_a"
		self.createTestObjs()

	def createTestObjs(self):
		self.fileAsListA = _getFileStrForTwoXyzStepsA().split("\n")
		self.trajSteps = _loadStartTrajStepsArrayDataTwoStepsA()

	@mock.patch("gen_basis_helpers.cp2k.parse_md_files._readFileIntoList")
	def testExpectedCaseA(self, mockedReadFileIntoList):
		mockedReadFileIntoList.side_effect = lambda *args,**kwargs: self.fileAsListA

		#Figure out what we expect
		expSteps = copy.deepcopy(self.trajSteps)
		expDataArrays = _loadExpectedArrayDataTwoXyzStepsA()
		expSteps[0].addExtraAttrDict( { "velocities":{"value":expDataArrays[0], "cmpType": "numericalArray"} } )
		expSteps[1].addExtraAttrDict( { "velocities":{"value":expDataArrays[1], "cmpType": "numericalArray"} } )

		#Run
		tCode._parseVelocitiesAndAddToTrajSteps(self.fakePathA, self.trajSteps)
		actSteps = self.trajSteps

		#Test
		mockedReadFileIntoList.assert_called_with(self.fakePathA)
		self.assertEqual(expSteps, actSteps)

	@mock.patch("gen_basis_helpers.cp2k.parse_md_files._readFileIntoList")
	def testExpectedWhenSecondStepIdxDoesntMatch(self, mockedReadFileIntoList):
		mockedReadFileIntoList.side_effect = lambda *args,**kwargs: self.fileAsListA
		self.trajSteps[1].step += 1

		#Figure out what we expect
		expSteps = copy.deepcopy(self.trajSteps)
		expDataArrays = _loadExpectedArrayDataTwoXyzStepsA()
		expSteps[0].addExtraAttrDict( { "velocities":{"value":expDataArrays[0], "cmpType": "numericalArray"} } )
		expSteps[1].addExtraAttrDict( { "velocities":{"value":None, "cmpType": "numericalArray"} } )

		#Run
		tCode._parseVelocitiesAndAddToTrajSteps(self.fakePathA, self.trajSteps)
		actSteps = self.trajSteps

		#Test
		mockedReadFileIntoList.assert_called_with(self.fakePathA)
		self.assertEqual(expSteps, actSteps)

	@mock.patch("gen_basis_helpers.cp2k.parse_md_files._readFileIntoList")
	def testExpectedForForcesA(self, mockedReadFileIntoList):
		mockedReadFileIntoList.side_effect = lambda *args,**kwargs: self.fileAsListA

		#Figure out what we expect
		expSteps = copy.deepcopy(self.trajSteps)
		expDataArrays = _loadExpectedArrayDataTwoXyzStepsA()
		expSteps[0].addExtraAttrDict( { "forces":{"value":expDataArrays[0], "cmpType": "numericalArray"} } )
		expSteps[1].addExtraAttrDict( { "forces":{"value":expDataArrays[1], "cmpType": "numericalArray"} } )

		#Run
		tCode._parseForcesXyzAndAddToTrajSteps(self.fakePathA, self.trajSteps)
		actSteps = self.trajSteps

		#Test
		mockedReadFileIntoList.assert_called_with(self.fakePathA)
		self.assertEqual(expSteps, actSteps)



def _getFileStrForTwoXyzStepsA():
	return """      3
 i =     1684, time =      842.000, E =      -551.1892718614
  O         0.0000273977        0.0001497837       -0.0000059469
  O         0.0003008562       -0.0000448785        0.0005541464
  O         0.0000779939       -0.0003770666        0.0001145407
      3
 i =     1687, time =      843.500, E =      -551.1958582717
  O         0.0000092244        0.0000987555       -0.0000157442
  O         0.0003614182       -0.0000414387        0.0005957346
  O         0.0001496509       -0.0004363731        0.0001015894"""

def _loadExpectedArrayDataTwoXyzStepsA():

	stepA = [ [0.0000273977,  0.0001497837, -0.0000059469],
	          [0.0003008562, -0.0000448785,  0.0005541464],
	          [0.0000779939, -0.0003770666,  0.0001145407] ]

	stepB = [ [0.0000092244,  0.0000987555, -0.0000157442],
	          [0.0003614182, -0.0000414387,  0.0005957346],
	          [0.0001496509, -0.0004363731,  0.0001015894] ]

	return [stepA, stepB]

def _loadStartTrajStepsArrayDataTwoStepsA():
	stepA = trajHelp.TrajStepFlexible(step=1684, time=842.0)
	stepB = trajHelp.TrajStepFlexible(step=1687, time=843.5)
	return stepA, stepB



