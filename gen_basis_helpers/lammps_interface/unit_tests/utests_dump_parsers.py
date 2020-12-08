
import plato_pylib.shared.ucell_class as uCellHelp
import unittest
import unittest.mock as mock

import gen_basis_helpers.analyse_md.traj_core as trajObjHelp
import gen_basis_helpers.lammps_interface.lammps_parsers as tCode

class TestParseLammpsDumpFile(unittest.TestCase):

	def setUp(self):
		self.timeStep = None
		self.typeIdxToEle = None
		self.createTestObjs()

	def createTestObjs(self):
		self.dumpFileAsListA = _loadSmallFullFileStrA().split("\n")
		self._loadExpectedObjsForFullStrA()

	@mock.patch("gen_basis_helpers.lammps_interface.lammps_parsers._getFileAsListFromInpPath")
	def testExpectedTrajA(self, mockGetFileAsList):
		mockGetFileAsList.side_effect = lambda *args: self.dumpFileAsListA
		actTrajObj = tCode.getTrajectoryFromLammpsDumpFile("fake_file_path")
		expTrajObj = self.expTrajObjA
		self.assertEqual(expTrajObj,actTrajObj)

	@mock.patch("gen_basis_helpers.lammps_interface.lammps_parsers._getFileAsListFromInpPath")
	def testExpectedTrajA_withTimeSteps(self, mockGetFileAsList):
		self.timeStep = 1
		self.createTestObjs()
		mockGetFileAsList.side_effect = lambda *args: self.dumpFileAsListA
		actTrajObj = tCode.getTrajectoryFromLammpsDumpFile("fake_file_path", timeStep=self.timeStep)
		expTrajObj = self.expTrajObjA
		self.assertEqual(expTrajObj, actTrajObj)

	@mock.patch("gen_basis_helpers.lammps_interface.lammps_parsers._getFileAsListFromInpPath")
	def testExpectedTrajA_withEleKeys(self, mockedGetFileAsList):
		self.typeIdxToEle = { "1":"Mg", "2": "O"}
		self.createTestObjs()
		mockedGetFileAsList.side_effect = lambda *args: self.dumpFileAsListA
		actTrajObj = tCode.getTrajectoryFromLammpsDumpFile("fake_file_path", typeIdxToEle=self.typeIdxToEle)
		expTrajObj = self.expTrajObjA
		self.assertEqual(expTrajObj, actTrajObj)
	
	def _loadExpectedObjsForFullStrA(self):
		lattVectsBoth = [ [9.63    , 0      , 0       ],
		                  [-4.81018, 8.33982, 0       ],
		                  [0       , 0      , 17.87895] ]

		cartCoordsA =           [ [0.40125, 2.08496, 1.11743, "1"],
		                          [0.0226 , 3.0105 , 1.11743, "2"],
		                          [1.39212, 2.21981, 1.11743, "2"],
		                          [3.61125, 2.08496, 1.11743, "1"],
		                          [3.2326 , 3.0105 , 1.11743, "2"],
		                          [4.60212, 2.21981, 1.11743, "2"] ]
		cartCoordsB =           [ [-0.426578, 1.73786, 0.96724 ,"1"],
		                          [0.148568 , 2.48075, 1.30372 ,"2"],
		                          [-1.43515 , 1.8717 , 1.23895 ,"2"],
		                          [4.47632  , 2.423  , 0.945189,"1"],
		                          [3.6439   , 2.54601, 1.57979 ,"2"],
		                          [4.44324  , 1.45043, 0.686697,"2"] ]

		if self.typeIdxToEle is not None:
			for idx,vals in enumerate(cartCoordsA):
				cartCoordsA[idx][-1] = self.typeIdxToEle[ cartCoordsA[idx][-1] ]
			for idx,vals in enumerate(cartCoordsB):
				cartCoordsB[idx][-1] = self.typeIdxToEle[ cartCoordsB[idx][-1] ]


		self.cellA = uCellHelp.UnitCell.fromLattVects(lattVectsBoth)
		self.cellB = uCellHelp.UnitCell.fromLattVects(lattVectsBoth)

		self.cellA.cartCoords = cartCoordsA
		self.cellB.cartCoords = cartCoordsB

		if self.timeStep is None:
			self.trajStepA = trajObjHelp.TrajStepBase(unitCell=self.cellA, step=0 , time=None) 
			self.trajStepB = trajObjHelp.TrajStepBase(unitCell=self.cellB, step=50, time=None)
		else:
			self.trajStepA = trajObjHelp.TrajStepBase(unitCell=self.cellA, step=0 , time=0) 
			self.trajStepB = trajObjHelp.TrajStepBase(unitCell=self.cellB, step=50, time=50*self.timeStep)

		self.expTrajObjA = trajObjHelp.TrajectoryInMemory([self.trajStepA, self.trajStepB])


def _loadSmallFullFileStrA():
	outStr = """ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
6
ITEM: BOX BOUNDS xy xz yz pp pp pp
-4.8101799999999999e+00 9.6300000000000008e+00 -4.8101799999999999e+00
0.0000000000000000e+00 8.3398199999999996e+00 0.0000000000000000e+00
0.0000000000000000e+00 1.7878950000000000e+01 0.0000000000000000e+00
ITEM: ATOMS id type x y z
1 1 0.40125 2.08496 1.11743
2 2 0.0226 3.0105 1.11743
3 2 1.39212 2.21981 1.11743
4 1 3.61125 2.08496 1.11743
5 2 3.2326 3.0105 1.11743
6 2 4.60212 2.21981 1.11743
ITEM: TIMESTEP
50
ITEM: NUMBER OF ATOMS
6
ITEM: BOX BOUNDS xy xz yz pp pp pp
-4.8101799999999999e+00 9.6300000000000008e+00 -4.8101799999999999e+00
0.0000000000000000e+00 8.3398199999999996e+00 0.0000000000000000e+00
0.0000000000000000e+00 1.7878950000000000e+01 0.0000000000000000e+00
ITEM: ATOMS id type x y z
1 1 -0.426578 1.73786 0.96724
2 2 0.148568 2.48075 1.30372
3 2 -1.43515 1.8717 1.23895
4 1 4.47632 2.423 0.945189
5 2 3.6439 2.54601 1.57979
6 2 4.44324 1.45043 0.686697
"""
	return outStr

