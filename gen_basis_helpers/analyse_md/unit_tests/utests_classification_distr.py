

import unittest

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.classification_distr_opt_objs as classDistrOptObjs
import gen_basis_helpers.analyse_md.traj_core as trajCoreHelp

import gen_basis_helpers.analyse_md.classification_distr as tCode

class TestGetClassDistrCountsOverTraj(unittest.TestCase):

	def setUp(self):
		#Geometry options
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.coordsA = [ [0,0,-1,"X"], [0,0,1,"O"], [0,1,3,"H"], [0,-1,3,"H"] ]
		self.coordsB = [ [0,0,-1,"X"], [0,0,6,"O"], [0,1,6,"H"], [0,-1,6,"H"] ] 

		#Options object
		self.binResObjA, self.binResObjB = [None, None] #Not actually needed for this; but non-optional arg for opts object 
		self.oxyIndices, self.hyIndices = [1], [[2,3]]
		self.distFilterIndices = [0]
		self.distRangeA = [0,2.1] #1 in first, 0 in second
		self.distRangeB = [2.9,8] #0 in first, 1 in second
		self.distRangeC = [0,10] #1 in both

		self.createTestObjs()

	def createTestObjs(self):
		#Geom/Trajectory
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellB = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords, self.cellB.cartCoords = self.coordsA, self.coordsB

		trajStepA = trajCoreHelp.TrajStepFlexible(unitCell=self.cellA, step=0, time=10)
		trajStepB = trajCoreHelp.TrajStepFlexible(unitCell=self.cellB, step=1, time=20)

		self.traj = trajCoreHelp.TrajectoryInMemory([trajStepA,trajStepB])

		#Option object
		currArgs = [self.binResObjA, self.oxyIndices, self.hyIndices, self.distFilterIndices, [self.distRangeA, self.distRangeB,self.distRangeC]]
		self.optObj = classDistrOptObjs.WaterCountTypesMinDistAndHBondSimpleOpts(*currArgs)

	def _runTestFunct(self):
		currArgs = [self.traj, self.optObj]
		return tCode.getClassDistrCountsOverTraj(*currArgs)

	def testExpectedValsA(self):
		expVals = [  (1,0,1), (0,1,1) ] 
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)




