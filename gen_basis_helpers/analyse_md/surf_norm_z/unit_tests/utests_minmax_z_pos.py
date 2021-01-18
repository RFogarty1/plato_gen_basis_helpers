
import itertools as it
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.traj_core as trajHelp
import gen_basis_helpers.analyse_md.surf_norm_z.get_minmax_positions as tCode

class TestExpectedMinMaxValsStandard(unittest.TestCase):

	def setUp(self):
		self.coordsA = [ [1,1,2,"X"], [1,1,4,"Y"], [1,1,5,"Y"],
		                 [1,1,7,"X"], [1,1,8,"Y"], [1,1,9,"Y"] ]
		self.coordsB = [ [1,1,3,"X"], [1,1,5,"Y"], [1,1,9,"Y"],
		                 [1,1,6,"X"], [1,1,8,"Y"], [1,1,7,"Y"] ]

		self.lattParams = [10,10,10]
		self.lattAngles = [90,90,90]
		self.groupX = [0,3]
		self.groupY = [1,2,4,5]
		self.currAtomGroup = self.groupX
		self.minZ, self.maxZ = None, None
		self.createTestObjs()

	def createTestObjs(self):
		#Create the unit cells
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellB = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coordsA
		self.cellB.cartCoords = self.coordsB

		#Create traj steps
		self.trajStepA = trajHelp.TrajStepBase(unitCell=self.cellA)
		self.trajStepB = trajHelp.TrajStepBase(unitCell=self.cellB)

		#create trajectories
		self.trajA = trajHelp.TrajectoryInMemory([self.trajStepA,self.trajStepB])

	def _runTestFunct(self):
		return tCode.getMinMaxZPositionsAtomGroups(self.trajA, self.currAtomGroup, minZ=self.minZ, maxZ=self.maxZ)

	def testMinMaxSimpleGroupX(self):
		self.currAtomGroup = self.groupX
		expVals = [ [2,7], [3,6] ]
		actVals = self._runTestFunct()
		self._checkExpAndActEqual(expVals,actVals)

	def testMinMaxGroupY_restrictedMin(self):
		self.currAtomGroup = self.groupY
		self.minZ = 4.5
		expVals = [ [5,9], [5,9] ]
		actVals = self._runTestFunct()
		self._checkExpAndActEqual(expVals,actVals)

	def testMinMaxGroupY_restrictedMax(self):
		self.currAtomGroup = self.groupY
		self.maxZ = 8.5
		expVals = [ [4,8], [5,8] ]
		actVals = self._runTestFunct()
		self._checkExpAndActEqual(expVals, actVals)

	def testMinMaxGroupY_restrictedMinMax(self):
		self.currAtomGroup = self.groupY
		self.maxZ = 8.5
		self.minZ = 4.5
		expVals = [ [5,8], [5,8] ]
		actVals = self._runTestFunct()
		self._checkExpAndActEqual(expVals, actVals)

	def _checkExpAndActEqual(self,expVals,actVals):
		self.assertEqual( len(expVals), len(actVals) )
		for exp,act in zip(expVals, actVals):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp,act)]


