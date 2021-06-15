
import copy
import itertools as it
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.binned_res as binResHelp
import gen_basis_helpers.analyse_md.traj_core as trajCoreHelp

import gen_basis_helpers.analyse_md.calc_distrib_core as tCode


class TestPopulateBinsWithRdfBetweenAtomGroups(unittest.TestCase):

	def setUp(self):
		#The cell
		self.lattParamsA, self.lattAnglesA = [10,10,10], [90,90,90]
		self.cartCoordsA = [ [5,5,5,"X"],
		                     [5,5,8,"Y"],
		                     [5,5,9,"Z"] ]

		#The bins
		self.binEdges = [0, 2.5, 5]

		#Options
		self.indicesA = [0]
		self.indicesB = [1,2]
		self.inpVolumes = None

		self.createTestObjs()

	def createTestObjs(self):
		#Create very simple trajectory [single step]
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParamsA, lattAngles=self.lattAnglesA)
		self.cellA.cartCoords = self.cartCoordsA
		trajStepA = trajCoreHelp.TrajStepFlexible(unitCell=self.cellA)

		self.trajA = trajCoreHelp.TrajectoryInMemory([trajStepA])
		self.binResObjA = binResHelp.BinnedResultsStandard.fromBinEdges(self.binEdges)

	def _runTestFunct(self):
		currArgs = [ self.trajA, [self.binResObjA], [self.indicesA], [self.indicesB] ]
		currKwargs = {"volumes":self.inpVolumes}  
		return tCode._populateBinsWithRdfBetweenAtomGroups(*currArgs, **currKwargs)

	def _loadExpectedValsA(self):
		expVols = [49.0873852123405, 441.786466911065]
		expRdfVals = [0,2.26353696841807] #Figured out in excel

		expBinObj = copy.deepcopy(self.binResObjA)
		expBinObj.binVals["counts"] = [0,2]
		expBinObj.binVals["rdf"], expBinObj.binVals["volume"] = expRdfVals, expVols 
		return expBinObj	
	
	def testExpectedValsForTrimerSingleTrajStep(self):
		expBinObj = self._loadExpectedValsA()
		self._runTestFunct()
		actBinObj = self.binResObjA
		self.assertEqual(expBinObj, actBinObj)

	def testExpectedUsingRdfOptionInterface(self):
		expBinObj = self._loadExpectedValsA()
		optsObj = tCode.CalcRdfOptions(self.binResObjA, self.indicesA, self.indicesB)
		tCode.populateRdfValsOnOptionObjs(self.trajA, [optsObj])
		self.assertEqual(expBinObj, self.binResObjA)

	@unittest.skip("")
	def testExpectedReverseIndices(self):
		self.assertTrue(False)

	@unittest.skip("")
	def testExpectedWhenVolumeReadIn(self):
		self.assertTrue(False)

class TestAddRdfToBinVals(unittest.TestCase):

	def setUp(self):
		self.nA = 2
		self.nB = 4
		self.binEdges = [0,2,4]
		self.nSteps = 10
		self.counts = [40, 100]
		self.volTotal = 1000
		self.countKey = "counts"
		self.createTestObjs()

	def createTestObjs(self):
		binVals = {self.countKey:self.counts}
		self.binRes = binResHelp.BinnedResultsStandard.fromBinEdges(self.binEdges, binVals=binVals)

	def _runTestFunct(self):
		args = [self.binRes, self.volTotal, self.nA, self.nB, self.nSteps]
		return tCode._addRdfToBinValsForBinsWithCounts(*args, countKey=self.countKey)

	def testExpectedValsA(self):
		#Calculated in excel; Note the ridic. wide bins make these values sorta silly
		expVolumes = [25.1327412287183, 226.194671058465]
		expRdfVals =[19.8943678864869, 5.52621330180192]
		expBinVals = copy.deepcopy(self.binRes)
		expBinVals.binVals["rdf"] = expRdfVals
		expBinVals.binVals["volume"] = expVolumes
		self._runTestFunct()
		self.assertEqual(expBinVals, self.binRes)


class TestMultiRdfBinnerFixedIndices(unittest.TestCase):

	def setUp(self):
		self.binEdgesA = [0,3]
		self.binEdgesB = [0,6]

		self.distMatrix = [ [0,2,4],
		                    [2,0,5],
		                    [4,5,0] ]
		self.indicesA = [0,1,2]
		self.indicesB = [0,1,2]
		self.inpTrajStep = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		#Create the test object
		self.binResA = binResHelp.BinnedResultsStandard.fromBinEdges(self.binEdgesA)
		self.binResB = binResHelp.BinnedResultsStandard.fromBinEdges(self.binEdgesB)
		self.binnerA = tCode._RdfBinnerFixedIndices(resBins=self.binResA, indicesA=self.indicesA, indicesB=self.indicesB)
		self.binnerB = tCode._RdfBinnerFixedIndices(resBins=self.binResB, indicesA=self.indicesA, indicesB=self.indicesB)
		self.testObj = tCode._MultiRdfBinnerFixedIndices([self.binnerA, self.binnerB])

	def _runTestFunct(self):
		self.testObj.updateCountsFromTrajStep(self.inpTrajStep)

	@mock.patch("gen_basis_helpers.analyse_md.calc_distrib_core.calcDistsHelp.calcDistanceMatrixForCell_minImageConv")
	def testExpCaseA(self, mockGetDistMatrix):
		mockGetDistMatrix.side_effect = lambda *args,**kwargs:self.distMatrix

		expBinResA, expBinResB = copy.deepcopy(self.binResA), copy.deepcopy(self.binResB)
		expBinResA.binVals = {"counts":[1]}
		expBinResB.binVals = {"counts":[3]}
		self._runTestFunct()
		actBinResA, actBinResB = self.binResA, self.binResB

		self.assertEqual(expBinResA, actBinResA)
		self.assertEqual(expBinResB, actBinResB)


class TestRdfBinnerFixedIndices(unittest.TestCase):

	def setUp(self):
		self.distMatrix = [ [0,2,4],
		                    [2,0,5],
		                    [4,5,0] ]
		self.binEdges = [0, 3.0, 6.1]
		self.binVals = None
		self.indicesA = [0]
		self.indicesB = [1,2]
		self.createTestObjs()

	def createTestObjs(self):
		self.binResA = binResHelp.BinnedResultsStandard.fromBinEdges(self.binEdges, binVals=self.binVals)
		currKwargs = {"resBins":self.binResA, "indicesA":self.indicesA, "indicesB":self.indicesB} 
		self.calcObjA = tCode._RdfBinnerFixedIndices(**currKwargs)

	def _runTestFunct_fromDistMatrix(self):
		self.calcObjA.updateCountsFromDistMatrix(self.distMatrix)

	def testExpectedResA_bothIndices_fromBlank(self):
		expBinVals = {"counts":[1,1]}
		expBinRes = copy.deepcopy(self.binResA)
		expBinRes.binVals = expBinVals
		self._runTestFunct_fromDistMatrix()
		actBinRes = self.binResA
		self.assertEqual(expBinRes, actBinRes)

	def testExpectedResB_allIndicesBoth(self):
		self.indicesA = [0,1,2]
		self.indicesB = [0,1,2]
		self.createTestObjs()
		expBinVals = {"counts":[1,2]}
		expBinRes = copy.deepcopy(self.binResA)
		expBinRes.binVals = expBinVals
		self._runTestFunct_fromDistMatrix()
		actBinRes = self.binResA
		self.assertEqual(expBinRes, actBinRes)


class TestExpectedBinValsFromFullDistMatrix(unittest.TestCase):

	def setUp(self):
		self.indicesA = None
		self.indicesB = None
		self.distMatrix = [ [0,2,3],
		                    [2,0,6],
		                    [3,6,0] ]

	def _runTestFunct(self):
		currKwargs = {"indicesA":self.indicesA, "indicesB":self.indicesB}
		return tCode._getRadialToBinValsFromFullDistMatrix(self.distMatrix, **currKwargs)

	def testExpectedNoIndicesGiven(self):
		expVals = [2,3,6]
		actVals = self._runTestFunct()
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(sorted(expVals),sorted(actVals))]

	def testBothSetsOfIndicesSpecified(self):
		self.indicesA = [0]
		self.indicesB = [1,2]
		expVals = [2,3]
		actVals = self._runTestFunct()
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(sorted(expVals),sorted(actVals))]

	def testOneSetOfIndicesSet(self):
		self.indicesA = [0]
		expVals = [2,3]
		actVals = self._runTestFunct()
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(sorted(expVals),sorted(actVals))]

	def testRaisesForNonSquareMatrix(self):
		self.distMatrix = [ [2,3], [4,5], [6,7] ]
		with self.assertRaises(AssertionError):
			self._runTestFunct()


