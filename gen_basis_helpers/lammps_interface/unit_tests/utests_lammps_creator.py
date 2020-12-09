
import collections
import unittest
import unittest.mock as mock

import gen_basis_helpers.lammps_interface.lammps_creator_simple as tCode

class TestLammpsCalcObjFactory(unittest.TestCase):

	def setUp(self):
		self.mockGeomObj = mock.Mock()
		self.mockPotentialObj = None
		self.mockEnsembleObj = None
		self.units = "real"
		self.atomStyle = "full"
		self.fileName = "fake_file"
		self.workFolder = "fake_path"
		self.thermoStyle = None
		self.printThermoEveryNSteps = None
		self.timeStep = None
		self.velocityObj = None
		self.nSteps = None
		self.dumpOpts = None
		self.boundaries = None
		self.walls = None
		self.eleToTypeIdx = {"O":1}
		self.createTestObjs()

	def createTestObjs(self):
		kwargs = {"lammpsGeomObj": self.mockGeomObj, "potentialObj":self.mockPotentialObj,
		          "ensembleObj": self.mockEnsembleObj, "units":self.units, "atomStyle":self.atomStyle,
		          "fileName":self.fileName, "timeStep":self.timeStep, "thermoStyle":self.thermoStyle,
		          "printThermoEveryNSteps":self.printThermoEveryNSteps, "velocityObj":self.velocityObj,
		          "nSteps":self.nSteps, "dumpOptions":self.dumpOpts, "workFolder":self.workFolder,
		          "boundaries":self.boundaries, "walls":self.walls}
		self.testObjA = tCode.LammpsCalcObjFactorySimple(**kwargs)

	def testExpectedDataDict(self):
		expDataDict = {"keyA":"valA"}
		self.mockGeomObj.getDataDictFunct.side_effect = lambda *args: expDataDict
		actDataDict = self.testObjA._getDataFileDict()
		self.mockGeomObj.getDataDictFunct.assert_called_with(self.mockGeomObj)
		self.assertEqual(expDataDict,actDataDict)

	@mock.patch("gen_basis_helpers.lammps_interface.lammps_creator_simple.LammpsCalcObjFactorySimple._getInitCommandDict")
	def testExpectedScriptDictPotentialOnly(self, mockedGetInitDict):
		mockedGetInitDict.side_effect = lambda *args: collections.OrderedDict()
		self.mockPotentialObj = mock.Mock()
		self.createTestObjs()
		expDict = collections.OrderedDict([["keyPotObj","valPotObj"]])
		self.mockPotentialObj.commandDict = expDict
		actDict = self.testObjA._getScriptFileDict()
		self.assertEqual(expDict,actDict)

	@mock.patch("gen_basis_helpers.lammps_interface.lammps_creator_simple.LammpsCalcObjFactorySimple._getInitCommandDict")
	def testExpectedScriptDictEnsembleOnly(self, mockedGetInitDict):
		mockedGetInitDict.side_effect = lambda *args: collections.OrderedDict()
		self.mockEnsembleObj = mock.Mock()
		mockFixStr = "fake fix string"
		self.mockEnsembleObj.fixStr = mockFixStr
		self.createTestObjs()
		expDict = collections.OrderedDict([["fix","1 " + mockFixStr]])
		actDict = self.testObjA._getScriptFileDict()
		self.assertEqual(expDict, actDict)

	def testExpectedScriptDictInitOptsOnly(self):
		expDictArgs = [ ["units",self.units], ["atom_style",self.atomStyle], ["read_data", self.fileName + ".data"] ]
		expDict = collections.OrderedDict( expDictArgs )
		actDict = self.testObjA._getScriptFileDict()
		self.assertEqual(expDict,actDict)

	@mock.patch("gen_basis_helpers.lammps_interface.lammps_creator_simple.LammpsCalcObjFactorySimple._getInitCommandDict")
	def testExpectedScriptDictSettingsOptsOnly(self, mockedGetInitDict):
		mockedGetInitDict.side_effect = lambda *args: collections.OrderedDict()
		expVelocityStr = "fake_velocity_comm_str"
		self.velocityObj = mock.Mock()
		self.velocityObj.commandStr = expVelocityStr
		self.thermoStyle = "fake_thermo_style"
		self.printThermoEveryNSteps = 50
		self.timeStep = 2
		self.createTestObjs()

		expDictArgs = [ ["velocity",expVelocityStr], ["timestep", "{:.1f}".format(self.timeStep)],
		                ["thermo_style",self.thermoStyle], ["thermo",str(self.printThermoEveryNSteps)] ] 
		expDict = collections.OrderedDict(expDictArgs)
		actDict = self.testObjA._getScriptFileDict()

		self.assertEqual(expDict, actDict)

	@mock.patch("gen_basis_helpers.lammps_interface.lammps_creator_simple.LammpsCalcObjFactorySimple._getInitCommandDict")
	def testExpectedScriptDictNStepsOnly(self, mockedGetInitDict):
		mockedGetInitDict.side_effect = lambda *args: collections.OrderedDict()
		self.nSteps = 300
		self.createTestObjs()
		expDict = collections.OrderedDict([ ["run", str(self.nSteps)] ])
		actDict = self.testObjA._getScriptFileDict()
		self.assertEqual(expDict,actDict)

	@mock.patch("gen_basis_helpers.lammps_interface.lammps_creator_simple.LammpsCalcObjFactorySimple._getInitCommandDict")
	def testExpectedScriptDictDumpOptsOnly(self, mockedGetInitDict):
		mockedGetInitDict.side_effect = lambda *args: collections.OrderedDict()
		expDumpOpts = collections.OrderedDict( [["keyA","fake_val_a"]] )
		self.dumpOpts = mock.Mock()
		self.dumpOpts.commandDict = expDumpOpts
		self.createTestObjs()
		actDict = self.testObjA._getScriptFileDict()
		self.assertEqual(expDumpOpts, actDict)

	def testExpectedScriptDictBoundariesOnly(self):
		self.units,self.atomStyle = None, None
		self.boundaries = ["p","p","p"]
		self.createTestObjs()
		expDict = collections.OrderedDict([["boundary","p p p"], ["read_data", self.fileName + ".data"]])
		actDict = self.testObjA._getScriptFileDict()
		self.assertEqual(expDict,actDict)

	@mock.patch("gen_basis_helpers.lammps_interface.lammps_creator_simple.LammpsCalcObjFactorySimple._getInitCommandDict")
	def testExpectedScriptDirWallsOnly(self, mockedGetInitDict):
		expCommA, expCommB = "comm_a", "comm_b"
		wallA, wallB = mock.Mock(), mock.Mock()
		wallA.fixStr = expCommA
		wallB.fixStr = expCommB
		self.walls = [wallA,wallB]
		mockedGetInitDict.side_effect = lambda *args: collections.OrderedDict()
		self.createTestObjs()
		expOutDict = collections.OrderedDict( [["fix","1 comm_a\nfix 2 comm_b"]] )
		actOutDict = self.testObjA._getScriptFileDict()
		self.assertEqual(expOutDict,actOutDict)


	@mock.patch("gen_basis_helpers.lammps_interface.lammps_creator_simple.calcObjHelp.LammpsCalcObjStandard")
	@mock.patch("gen_basis_helpers.lammps_interface.lammps_creator_simple.LammpsCalcObjFactorySimple._getDataFileDict")
	@mock.patch("gen_basis_helpers.lammps_interface.lammps_creator_simple.LammpsCalcObjFactorySimple._getScriptFileDict")
	def testCreateFromSelf(self, mockedGetScriptDict, mockedGetDataDict, mockedCalcObj):
		expDataFileDict, expScriptFileDict = "fake_dict_data_file", "fake_dict_script_file"
		expCalcObj = mock.Mock()
		mockedGetScriptDict.side_effect = lambda *args:expScriptFileDict
		mockedGetDataDict.side_effect = lambda *args:expDataFileDict
		mockedCalcObj.side_effect = lambda *args,**kwargs: expCalcObj
		self.mockGeomObj.eleToTypeIdx = self.eleToTypeIdx

		actCalcObj = self.testObjA.create()

		mockedCalcObj.assert_called_with(self.workFolder, self.fileName, expDataFileDict, expScriptFileDict,typeIdxToEle={1: 'O'})
		self.assertEqual(expCalcObj, actCalcObj)



