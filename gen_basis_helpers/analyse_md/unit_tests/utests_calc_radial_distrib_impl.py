
import copy
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.binned_res as binResHelp
import gen_basis_helpers.analyse_md.traj_core as trajCoreHelp
import gen_basis_helpers.shared.plane_equations as planeEqnHelp

import gen_basis_helpers.analyse_md.calc_radial_distrib_impl as tCode

class TestPopulateBinsWithPlanarRdfVals(unittest.TestCase):

	def setUp(self):
		#The cell
		self.lattParamsA, self.lattAnglesA = [10,10,10], [90,90,90]
		self.cartCoordsA = [ [5,5,9,"X"],
		                     [5,5,1,"Y"],
		                     [6,6,3,"Z"] ]

		#The bins
		self.binEdges = [0, 2, 4]

		#Options
		self.indicesA = [0,1,2]
		self.inpVolumes = None
		self.planeEqn = None

		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParamsA, lattAngles=self.lattAnglesA)
		self.cellA.cartCoords = self.cartCoordsA
		trajStepA = trajCoreHelp.TrajStepFlexible(unitCell=self.cellA)

		self.trajA = trajCoreHelp.TrajectoryInMemory([trajStepA])
		self.binResObjA = binResHelp.BinnedResultsStandard.fromBinEdges(self.binEdges)

	def _runTestFunct(self):
		currArgs = [ self.trajA, [self.binResObjA], [self.indicesA] ]
		currKwargs = {"volumes":self.inpVolumes, "planeEqn":self.planeEqn}
		return tCode._populateBinsWithPlanarRdfVals(*currArgs, **currKwargs)

	def _loadExpBinsCaseA(self):
		expVols = [200,200]
		expRdfVals = [10/3, 5/3]
		expBinObj = copy.deepcopy(self.binResObjA)

		expBinObj.binVals["counts"] = [2,1]
		expBinObj.binVals["rdf"], expBinObj.binVals["volume"] = expRdfVals, expVols
		return expBinObj		

	def testExpectedCaseA(self):
		expVols = [200,200]
		expRdfVals = [10/3, 5/3]
		expBinObj = copy.deepcopy(self.binResObjA)

		expBinObj.binVals["counts"] = [2,1]
		expBinObj.binVals["rdf"], expBinObj.binVals["volume"] = expRdfVals, expVols
		self._runTestFunct()
		actBinObj = self.binResObjA

		self.assertEqual(expBinObj, actBinObj)

	def testExpectedCaseA_commandInterface(self):
		optsObj = tCode.CalcPlanarRdfOptions(self.binResObjA, self.indicesA, planeEqn=self.planeEqn, volume=self.inpVolumes)
		expBinObj = self._loadExpBinsCaseA()
		tCode.populatePlanarRdfsFromOptionsObjs(self.trajA, [optsObj])
		actBinObj = self.binResObjA
		self.assertEqual(expBinObj, actBinObj)


class TestPlanarRdfMultiBinnerFixedIndices(unittest.TestCase):

	def setUp(self):
		#The cell
		self.lattParamsA, self.lattAnglesA = [10,10,10], [90,90,90]
		self.cartCoordsA = [ [5,5,9,"X"],
		                     [5,5,1,"Y"],
		                     [6,6,3,"Z"] ]

		#The indices to use for each binner
		self.indicesA = [0,1,2]
		self.indicesB = [0,2]

		#Bin edges;
		self.binEdgesA = [0, 2, 4]
		self.binEdgesB = [0, 2, 4]

		#Plane equations to use for each binner
		self.planeEqnA = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,0)
		self.planeEqnB = planeEqnHelp.ThreeDimPlaneEquation(0,0,1,0)

		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParamsA, lattAngles=self.lattAnglesA)
		self.cellA.cartCoords = self.cartCoordsA
		self.trajStepA = trajCoreHelp.TrajStepFlexible(unitCell=self.cellA)

		#Sort the bins
		self.binsA = binResHelp.BinnedResultsStandard.fromBinEdges(self.binEdgesA)
		self.binsB = binResHelp.BinnedResultsStandard.fromBinEdges(self.binEdgesB)

		#Sort the single binners
		self.binnerA = tCode._PlanarRdfBinnerFixedIndices(resBins=self.binsA, indices=self.indicesA, planeEqn=self.planeEqnA)
		self.binnerB = tCode._PlanarRdfBinnerFixedIndices(resBins=self.binsB, indices=self.indicesB, planeEqn=self.planeEqnB)

		#Create the test object
		self.testObj = tCode._PlanarRdfMultiBinnerFixedIndices([self.binnerA, self.binnerB])

	def _runTestFunct(self):
		return self.testObj.updateCountsFromTrajStep(self.trajStepA)

	def testExpectedCaseA(self):
		"""  Two separate with different indices """
		expBinsA, expBinsB = copy.deepcopy(self.binsA), copy.deepcopy(self.binsB)
		expBinsA.binVals["counts"], expBinsB.binVals["counts"] = [2,1], [1,1]

		self._runTestFunct()
		actBinsA, actBinsB = self.binsA, self.binsB

		self.assertEqual(expBinsA, actBinsA)
		self.assertEqual(expBinsB, actBinsB)

	def testExpected_middleIdxInNeither(self):
		self.indicesA = [0]
		self.indicesB = [0,2]
		self.createTestObjs()
		expBinsA, expBinsB = copy.deepcopy(self.binsA), copy.deepcopy(self.binsB)
		expBinsA.binVals["counts"], expBinsB.binVals["counts"] = [1,0], [1,1]

		self._runTestFunct()
		actBinsA, actBinsB = self.binsA, self.binsB

		self.assertEqual(expBinsA, actBinsA)
		self.assertEqual(expBinsB, actBinsB)

	def testRaisesForNonEqualParallelPlanes(self):
		self.planeEqnB = planeEqnHelp.ThreeDimPlaneEquation(0,1,1,2)
		self.createTestObjs()
		with self.assertRaises(AssertionError):
			self._runTestFunct()





