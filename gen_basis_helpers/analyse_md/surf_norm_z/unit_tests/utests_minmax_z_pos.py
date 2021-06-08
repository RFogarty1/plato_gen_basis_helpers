
import itertools as it
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.traj_core as trajHelp
import gen_basis_helpers.analyse_md.surf_norm_z.get_minmax_positions as tCode


class TestGetSurfaceThicknessesFromMinMaxData(unittest.TestCase):

	def setUp(self):
		self.valsA = [ [0,4], [9,12] ]

	def _runTestFunct(self):
		return tCode.getLayerThicknessDataFromMinMaxData(self.valsA)

	def testSimpleCaseA(self):
		expVals = [4,3]
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)


class TestGetAverageMinMaxForSubsetOfGroups(unittest.TestCase):

	def setUp(self):
		self.traj = createTestTrajA()
		self.groupIndices = [x for x in range(6)]
		self.nMinZ, self.nMaxZ = None, None

	def _runTestFunct(self):
		args   = [self.traj, self.groupIndices]
		kwargs = {"nMinZ":self.nMinZ, "nMaxZ":self.nMaxZ}
		return tCode.getAverageMinMaxZPositionsForSubsetOfAtomGroups(*args, **kwargs)

	def testExpCaseA_fullGroup(self):
		expValA = sum([2,4,5,7,8,9])/6
		expValB = sum([3,5,9,6,8,7])/6
		expVals = [ [expValA, expValA], [expValB,expValB] ]
		actVals = self._runTestFunct()
		self._checkExpAndActEqual(expVals, actVals)
	
	def testExpCaseA_fullGroup_nMax2(self):
		self.nMaxZ = 2

		#Figure out expected
		expMinA = sum([2,4,5,7,8,9])/6
		expMinB = sum([3,5,9,6,8,7])/6
		expMaxA = (8+9)/2
		expMaxB = (8+9)/2
		expVals = [ [expMinA, expMaxA], [expMinB, expMaxB] ]

		#Test
		actVals = self._runTestFunct()
		self._checkExpAndActEqual(expVals, actVals)

	def testExpCase_fullGroup_nMin2(self):
		self.nMinZ = 2

		#Figure out expected
		expMaxA, expMaxB = sum([2,4,5,7,8,9])/6, sum([3,5,9,6,8,7])/6
		expMinA, expMinB = (2+4)/2, (3+5)/2
		expVals = [ [expMinA,expMaxA], [expMinB,expMaxB] ]

		#Test
		actVals = self._runTestFunct()
		self._checkExpAndActEqual(expVals, actVals)

	def testHalfGroup_nMinAndMax2(self):
		self.groupIndices = [0,1,2]
		self.nMinZ, self.nMaxZ = 2,2

		#Figure out expected
		expMinA, expMinB = (2+4)/2, (3+5)/2
		expMaxA, expMaxB = (4+5)/2, (5+9)/2
		expVals = [ [expMinA,expMaxA], [expMinB,expMaxB] ]

		#Test
		actVals = self._runTestFunct()
		self._checkExpAndActEqual(expVals,actVals)

	def _checkExpAndActEqual(self,expVals,actVals):
		self.assertEqual( len(expVals), len(actVals) )
		for exp,act in zip(expVals, actVals):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp,act)]


class TestExpectedMinMaxValsStandard(unittest.TestCase):

	def setUp(self):
		self.lattParams = [10,10,10]
		self.lattAngles = [90,90,90]
		self.groupX = [0,3]
		self.groupY = [1,2,4,5]
		self.currAtomGroup = self.groupX
		self.minZ, self.maxZ = None, None
		self.createTestObjs()

	def createTestObjs(self):
		self.trajA = createTestTrajA()

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


class TestGetAverageZPosForGroupOfAtoms(unittest.TestCase):

	def setUp(self):
		self.atomGroup = [1,3]
		self.minZ, self.maxZ = None, None
		self.trajA = createTestTrajA()

	def _runTestFunct(self):
		return tCode.getAverageZPositionForAtomGroup( self.trajA, self.atomGroup )

	def testExpectedValA(self):
		expVals = [ (4+7)/2, (5+6)/2 ]
		actVals = self._runTestFunct()
		[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expVals, actVals)]


class TestGetZPositionsForAtomIndices(unittest.TestCase):

	def setUp(self):
		self.inpIndices = [0,2,3]
		self.times = [10,20]
		self.retTimeVsZ = False
		self.deltaZ = False
		self.createTestObjs()

	def createTestObjs(self):
		self.trajA = createTestTrajA()
		for idx,step in enumerate(self.trajA):
			step.time = self.times[idx]

	def _runTestFunct(self):
		args = [self.trajA, self.inpIndices]
		kwargs = {"retTimeVsZ":self.retTimeVsZ, "deltaZ":self.deltaZ}
		return tCode.getZPositionsForAtomIndices(*args, **kwargs)

	def testExpectedValsA(self):
		expCoords = [ [2, 3],
		              [5, 9],
		              [7, 6] ]
		actCoords = self._runTestFunct()
		self.assertEqual( len(expCoords), len(actCoords) )
		for exp,act in zip(expCoords,actCoords):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp,act)]

	def testExpectedTimeVsZPos(self):
		self.retTimeVsZ = True
		expCoords = [ [ [10,2], [20,3] ],
		              [ [10,5], [20,9] ],
		              [ [10,7], [20,6] ] ]
		actCoords = self._runTestFunct()
		self.assertEqual( len(expCoords), len(actCoords) )

		for exp, act in zip(expCoords, actCoords):
			for expA,actA in it.zip_longest(exp,act):
				self.assertAlmostEqual(expA[0],actA[0])
				self.assertAlmostEqual(expA[1],actA[1])

	def testExpectedDeltaZPos(self):
		self.deltaZ = True
		expCoords = [ [0, 1],
		              [0, 4],
		              [0,-1] ]
		actCoords = self._runTestFunct()
		self.assertEqual( len(expCoords), len(actCoords) )
		for exp,act in zip(expCoords,actCoords):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp,act)]


def createTestTrajA():

	#Cells
	coordsA = [ [1,1,2,"X"], [1,1,4,"Y"], [1,1,5,"Y"],
	            [1,1,7,"X"], [1,1,8,"Y"], [1,1,9,"Y"] ]

	coordsB = [ [1,1,3,"X"], [1,1,5,"Y"], [1,1,9,"Y"],
	            [1,1,6,"X"], [1,1,8,"Y"], [1,1,7,"Y"] ]

	lattParams = [10,10,10]
	lattAngles = [90,90,90]

	#Create the unit cells
	cellA = uCellHelp.UnitCell(lattParams=lattParams, lattAngles=lattAngles)
	cellB = uCellHelp.UnitCell(lattParams=lattParams, lattAngles=lattAngles)
	cellA.cartCoords = coordsA
	cellB.cartCoords = coordsB

	#Create traj steps
	trajStepA = trajHelp.TrajStepBase(unitCell=cellA)
	trajStepB = trajHelp.TrajStepBase(unitCell=cellB)

	#create trajectories
	trajA = trajHelp.TrajectoryInMemory([trajStepA,trajStepB])

	return trajA

