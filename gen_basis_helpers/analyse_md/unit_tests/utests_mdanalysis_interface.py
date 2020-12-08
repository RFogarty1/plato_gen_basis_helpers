
import unittest
import unittest.mock as mock

import numpy as np

import plato_pylib.shared.ucell_class as uCellHelp
import gen_basis_helpers.analyse_md.traj_core as trajHelp
import gen_basis_helpers.analyse_md.mdanalysis_interface as tCode


class TestGetAtomicTopologyObjFromTrajObj(unittest.TestCase):

	def setUp(self):
		self.stepA = 0
		self.timeA = None
		self.lattParamsA, self.lattAnglesA = [10,10,10], [90,90,90]
		self.cartCoordsA = [ [5,5,5,"X"], [7,7,7,"Y"], [8,8,8,"Z"] ]
		self.createTestObjs()
		
	def createTestObjs(self):
		self._createCells()	
		self.trajStepA = trajHelp.TrajStepBase(unitCell=self.cellA, step=self.stepA, time=self.timeA)
		trajSteps = [self.trajStepA]
		self.testObjA = trajHelp.TrajectoryInMemory(trajSteps)

	def _createCells(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParamsA, lattAngles=self.lattAnglesA)
		self.cellA.cartCoords = self.cartCoordsA
	
	#Mock topology obj and AtomAttr
	@mock.patch("gen_basis_helpers.analyse_md.mdanalysis_interface.mdAnalysisLib.core.topologyattrs.Atomnames")
	@mock.patch("gen_basis_helpers.analyse_md.mdanalysis_interface.mdAnalysisLib.core.topology.Topology")
	def testExpectedCallsForSimpleInputA(self, mockedTopologyCls, mockedAtomNamesAttr):
		#Deal with mock objects
		expTopologyObj = mock.Mock()
		expMockAtomNameAttr = mock.Mock()
		mockedTopologyCls.side_effect = lambda *args,**kwargs: expTopologyObj
		mockedAtomNamesAttr.side_effect = lambda *args,**kwargs: expMockAtomNameAttr

		#Get expected vals
		expNAtoms, expNResids, expNSegments = [3,3,1]
		expResIDs = [0,1,2]
		expSegIDs = [0,0,0]
		expAttrs = [expMockAtomNameAttr]

		#Run + check expected calls made
		actTopologyObj = tCode.getSimpleAtomicTopologyFromTrajObj(self.testObjA)
		mockedAtomNamesAttr.assert_called_with(["X","Y","Z"])
		mockedTopologyCls.assert_called_with(expNAtoms, n_res=expNResids, n_seg=expNSegments, attrs=expAttrs,
		                                     atom_resindex=expResIDs, residue_segindex=expSegIDs)
		self.assertEqual(expTopologyObj, actTopologyObj)


class TestGetMemoryReaderFromTrajObj(unittest.TestCase):


	def setUp(self):
		self.stepA, self.stepB = 0, 50
		self.timeA, self.timeB = None, None
		self.lattParamsA, self.lattAnglesA = [10,10,10], [90,90,90]
		self.cartCoordsA = [ [5,5,5,"X"], [7,7,7,"X"] ]
		self.cartCoordsB = [ [6,6,6,"X"], [8,8,8,"X"] ]
		self.createTestObjs()

	def createTestObjs(self):
		self._createCells()	
		self.trajStepA = trajHelp.TrajStepBase(unitCell=self.cellA, step=self.stepA, time=self.timeA)
		self.trajStepB = trajHelp.TrajStepBase(unitCell=self.cellB, step=self.stepB, time=self.timeB)
		trajSteps = [self.trajStepA, self.trajStepB]
		self.testObjA = trajHelp.TrajectoryInMemory(trajSteps)

	def _createCells(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParamsA, lattAngles=self.lattAnglesA)
		self.cellB = uCellHelp.UnitCell(lattParams=self.lattParamsA, lattAngles=self.lattAnglesA)
		self.cellA.cartCoords = self.cartCoordsA
		self.cellB.cartCoords = self.cartCoordsB

	@mock.patch("gen_basis_helpers.analyse_md.mdanalysis_interface.mdAnalysisLib.coordinates.memory.MemoryReader")
	def testExpectedArgsPassedSimpleCaseA(self, mockedMemReader):
		#Setup
		expOutput = mock.Mock()
		mockedMemReader.side_effect = lambda *args, **kwargs: expOutput
		expCoords = self._getExpectedCoordinateArrayA()
		expDims = [ self.cellA.getLattParamsList() + self.cellA.getLattAnglesList(),
		            self.cellA.getLattParamsList() + self.cellA.getLattAnglesList() ]

		#Run
		actOutput = tCode.getMDAMemoryReaderTrajFromTrajInMemoryObj(self.testObjA)
		actCoords, kwargs = mockedMemReader.call_args
		actDims = kwargs["dimensions"]

		#Check equality
		self.assertTrue( np.allclose(np.array(expDims), np.array(actDims)) )
		self.assertTrue( np.allclose(expCoords,actCoords) )
		self.assertEqual(expOutput, actOutput)

	@unittest.skip("")
	@mock.patch("gen_basis_helpers.analyse_md.mdanalysis_interface.mdAnalysisLib.coordinates.memory.MemoryReader")
	def testExpectedArgsPassedWhenTimeStepsDefined(self, mockedMemReader):
		self.assertTrue(False)


	def _getExpectedCoordinateArrayA(self):
		cartListA = [x[:3] for x in self.cartCoordsA]
		cartListB = [x[:3] for x in self.cartCoordsB]
		outArray = np.array( [np.array(cartListA), np.array(cartListB)] )
		return outArray
