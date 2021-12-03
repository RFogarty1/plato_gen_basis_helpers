
import itertools as it
import unittest

import numpy as np

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.classification_distr_opt_objs as classDistrOptObjs
import gen_basis_helpers.analyse_md.distr_opt_objs as distrOptObjHelp
import gen_basis_helpers.analyse_md.filtered_atom_combo_opt_objs as filteredOptObjHelp
import gen_basis_helpers.analyse_md.traj_core as trajCoreHelp

import gen_basis_helpers.shared.plane_equations as planeEqnHelp

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


class TestFilteredClassDistrScatterData(unittest.TestCase):

	def setUp(self):
		#Geometry options; coordsA has two hydroxyl while coordsB has only 1
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.coordsA = [ [0,0,0,"X"], [4,4,4,"O"], [5,4,4,"H"],
		                 [7,7,7.1,"O"], [7,7,8,"H"], [3,4,5,"X"] ]

		#Something about this geometry stops the classifier working????
		self.coordsB = [ [0,0,0,"X"], [5,5,5,"O"], [6,5,5,"H"],
		                 [7,7,7,"O"], [2,2,2,"H"], [5,7,8,"X"] ]

		self.timeA, self.timeB = 10, 20

		#Not really neeeded but a non-optional argument somewhere
		self.binResObjA, self.binResObjB = [None,None]

		#Classifier options
		self.oxyIndices = [1,3] 
		self.hyIndices = [2,4]
		self.maxOOHBond = 1.5
		self.nNebs = [1]

		#Options for the planar distribution
		currArgs = [ self.binResObjA, self.hyIndices] #Neither arg should matter; hence i set the wrong value for the 2nd on purpose
		self.distrOptObjs = [distrOptObjHelp.CalcPlanarDistOptions(*currArgs, planeEqn=planeEqnHelp.ThreeDimPlaneEquation(0,0,1,1))]

		#Options for the filtered object thing
		self.useGroups = [ [0] ]
		self.useNonHyIdx, self.useIdxEach = True, 0

		#Options for the actual function
		self.foldOneDimData = True


		self.createTestObjs()

	def createTestObjs(self):
		#Geom/Trajectory
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellB = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)

		self.cellA.cartCoords, self.cellB.cartCoords = self.coordsA, self.coordsB

		trajStepA = trajCoreHelp.TrajStepFlexible(unitCell=self.cellA, step=0, time=self.timeA)
		trajStepB = trajCoreHelp.TrajStepFlexible(unitCell=self.cellB, step=1, time=self.timeB)

		self.traj = trajCoreHelp.TrajectoryInMemory([trajStepA,trajStepB])

		#Create the classification object options
		currArgs = [ [self.binResObjA], self.oxyIndices, self.hyIndices]
		currKwargs = {"maxOHDist":self.maxOOHBond,"nNebs":self.nNebs}
		hydroxylClassifier = classDistrOptObjs.WaterDerivativeBasedOnDistanceClassifierOptsObj(*currArgs, **currKwargs)


		#
		currArgs = [self.oxyIndices, self.hyIndices, hydroxylClassifier, self.distrOptObjs, self.useGroups ]
		currKwargs = {"useNonHyIdx":self.useNonHyIdx, "useIdxEach":self.useIdxEach}
		self.filteredOptObj = filteredOptObjHelp.WaterDerivativeFilteredOptsObj_simple(*currArgs, **currKwargs)

	def _runTestFunct(self):
		return tCode.getFilteredScatterDataOverTime(self.traj, self.filteredOptObj, foldOneDimData=self.foldOneDimData)

	#NOTE: Make all values n-dimensional I guess; can add a fold 1-d data keyword for convenience
	#(It should just mean we return floats in y-columns rather than a len-1 list)
	def testExpectedValsA(self):
		expVals = [ [self.timeA,3], [self.timeA,3.9], [self.timeB,4] ]
		actVals = self._runTestFunct()
		self.assertTrue( np.allclose( np.array(expVals), np.array(actVals) ) )

	def testExpected_twoIdenticalDistrOpts(self):
		self.useGroups = [ [0], [0] ]
		currArgs = [ self.binResObjA, self.hyIndices] #Neither arg should matter; hence i set the wrong value for the 2nd on purpose
		self.distrOptObjs = [distrOptObjHelp.CalcPlanarDistOptions(*currArgs, planeEqn=planeEqnHelp.ThreeDimPlaneEquation(0,0,1,1)),
		                     distrOptObjHelp.CalcPlanarDistOptions(*currArgs, planeEqn=planeEqnHelp.ThreeDimPlaneEquation(0,0,1,1)) ]
		self.createTestObjs()

		expVals = [ [self.timeA, (3,3)], [self.timeA, (3.9,3.9)], [self.timeB,(4,4)] ]
		actVals = self._runTestFunct()

		for currExpVal, currActVal in it.zip_longest(expVals, actVals):
			expTime, actTime = currExpVal[0], currActVal[0]
			expVals, actVals = currExpVal[1], currActVal[1]
			self.assertAlmostEqual(expTime,actTime)
			[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expVals,actVals)]	


if __name__ == '__main__':
	unittest.main()
