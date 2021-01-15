
import os
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.traj_core as trajHelp
import gen_basis_helpers.analyse_md.traj_io as tCode

class TestDumpExtendedXYZ(unittest.TestCase):

	def setUp(self):
		self.cellA = uCellHelp.UnitCell(lattParams=[7,8,9],lattAngles=[90,90,90])
		self.cellB = uCellHelp.UnitCell(lattParams=[4,5,6],lattAngles=[90,90,90])
		self.coordsA = [ [1,2,3,"Mg"], [4,5,6,"Mg"] ]
		self.coordsB = [ [2,2,2,"Mg"], [3,3,3,"Mg"], [1,1,1,"Mg"] ]
		self.stepA, self.stepB = 4, 8
		self.timeA, self.timeB = 12, 24
		self.tempPath = "temp_file.exyz"
		self.createTestObjs()

	def tearDown(self):
		try:
			os.remove(self.tempPath)
		except FileNotFoundError:
			pass

	def createTestObjs(self):
		self.cellA.cartCoords = self.coordsA
		self.cellB.cartCoords = self.coordsB
		self.trajStepA = trajHelp.TrajStepBase(unitCell=self.cellA, step=self.stepA, time=self.timeA)
		self.trajStepB = trajHelp.TrajStepBase(unitCell=self.cellB, step=self.stepB, time=self.timeB)
		self.testTrajA = trajHelp.TrajectoryInMemory([self.trajStepA, self.trajStepB])

	def testWriteAndReadLeadToSameTraj(self):
		tCode.writeTrajToSimpleExtendedXyzFormat(self.testTrajA, self.tempPath)
		expTraj = self.testTrajA
		actTraj = tCode.readTrajFromSimpleExtendedXyzFormat(self.tempPath)
		self.assertEqual(expTraj,actTraj)


class TestParseExtendedXYZ(unittest.TestCase):

	def setUp(self):
		self.inpPathA = "fake_path"
		self.fileListA = _loadTestExtendedXYZStrA().split("\n")

	@mock.patch("gen_basis_helpers.analyse_md.traj_io._readFileIntoList")
	def testExpectedTrajReturnedA(self, mockGetFileList):
		#Set out mock
		mockGetFileList.side_effect = lambda *args,**kwargs: self.fileListA

		#Figure out expected trajectory
		cellA = uCellHelp.UnitCell.fromLattVects([[6.06,0,0],[-3, 5.2, 0], [0, 0, 9.8]])
		cellB = uCellHelp.UnitCell.fromLattVects([[5,0,0], [0,6,0], [0,0,7]])
		coordsA = [ [1,1,1,"Mg"], [3,3,3,"X"] ]
		coordsB = [ [2,3,4,"Mg"], [4,5,6,"X"] ]
		cellA.cartCoords, cellB.cartCoords = coordsA, coordsB
		stepA = trajHelp.TrajStepBase(unitCell=cellA)
		stepB = trajHelp.TrajStepBase(unitCell=cellB)
		expTraj = trajHelp.TrajectoryInMemory([stepA,stepB])

		#Run + test
		actTraj = tCode.readTrajFromSimpleExtendedXyzFormat(self.inpPathA)
		mockGetFileList.assert_called_with(self.inpPathA)
		self.assertEqual(expTraj,actTraj)



def _loadTestExtendedXYZStrA():
	return """2
Lattice="6.06 0.0 0.0 -3.0 5.2 0.0 0.0 0.0 9.8" Properties=species:S:1:pos:R:3
 Mg  1 1 1
  X  3 3 3
2
Lattice="5 0 0 0 6 0 0 0 7"
 Mg 2 3 4
  X 4 5 6"""

