
import os
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp
import gen_basis_helpers.analyse_md.traj_core as trajHelp
import gen_basis_helpers.db_help.query_help as tCode

class TestGetMDTrajFromRecords(unittest.TestCase):

	def setUp(self):
		self.stepA, self.stepB = 0, 50
		self.timeA, self.timeB = 0, 500
		self.unitCellA = uCellHelp.UnitCell(lattParams=[10,10,10], lattAngles=[90,90,90])
		self.unitCellB = uCellHelp.UnitCell(lattParams=[1,2,3], lattAngles=[90,90,90])
		self.folderPathA = "/fake/path/to/folder"
		self.fileNameA = "file_name_a"
		self.folderKeyA = "md_traj_folder"
		self.fileKeyA = "md_traj_file"
		self.createTestObjs()

	def createTestObjs(self):
		self.trajStepA = trajHelp.TrajStepBase(unitCell=self.unitCellA, step=self.stepA, time=self.timeA)
		self.trajStepB = trajHelp.TrajStepBase(unitCell=self.unitCellB, step=self.stepB, time=self.timeB)
		self.testTrajA = trajHelp.TrajectoryInMemory( [self.trajStepA, self.trajStepB] )
		self.testRecord = { self.folderKeyA:self.folderPathA, self.fileKeyA:self.fileNameA }

	@mock.patch("gen_basis_helpers.db_help.query_help.trajHelp.readTrajObjFromFileToTrajectoryInMemory")
	def testGetFullTraj(self, mockedReadTrajFromFile):
		mockedReadTrajFromFile.side_effect = lambda *args,**kwrgs: self.testTrajA
		expObj = self.testTrajA
		actObj = tCode.getTrajectoryFromRecord(self.testRecord)
		mockedReadTrajFromFile.assert_called_with( os.path.join(self.folderPathA,self.fileNameA) )
		self.assertEqual(expObj, actObj)

	@mock.patch("gen_basis_helpers.db_help.query_help.trajHelp.readLastTrajStepFromFile")
	def testGetLastTrajStepFromRecord(self, mockedReadLastTraj):
		expObj = self.trajStepB
		mockedReadLastTraj.side_effect = lambda *args,**kwargs: expObj
		actObj = tCode.getFinalTrajStepFromRecord(self.testRecord)
		mockedReadLastTraj.assert_called_with( os.path.join(self.folderPathA,self.fileNameA) )
		self.assertEqual(expObj, actObj)


