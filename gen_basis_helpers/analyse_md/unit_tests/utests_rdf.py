
import math
import copy
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.traj_core as trajHelp
import gen_basis_helpers.analyse_md.calc_rdfs as tCode

class TestEleToEleRdfWaterDimer(unittest.TestCase):

	def setUp(self):
		self.oxyDistA = 2
		self.oxyDistB = 4
		self.nBins = 2
		self.range = [0,6]
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = self._getCellFromOxyDist(self.oxyDistA)
		self.cellB = self._getCellFromOxyDist(self.oxyDistB)
		trajA = trajHelp.TrajStepBase(unitCell=self.cellA, step=0)
		trajB = trajHelp.TrajStepBase(unitCell=self.cellB, step=50)
		self.trajObj = trajHelp.TrajectoryInMemory([trajA,trajB])

	def _getCellFromOxyDist(self, oxyDist):
		outCell = uCellHelp.UnitCell(lattParams=[30,30,30], lattAngles=[90,90,90])
		outCoords = [ [14, 14, 14        ,"H"],
		              [13, 13, 13        ,"H"],
		              [16, 16, 16        ,"H"],
		              [17, 17, 17        ,"H"],
		              [15, 15, 15        ,"O"],
		              [15, 15, 15+oxyDist,"O"] ]
		outCell.cartCoords = outCoords
		return outCell

	def testExpectedRdfForWaterDimerInLargeBoxTwoSnapshots_twoBins(self):
		#Construct the expected object (need to 2x check we know what boxes MDAnalysis will choose)
		nRelevantAtomsInCell = 2
		expDists = [1.5,4.5]
		expEdges = [0,3,6]
		expCounts = [2,2]
		volumeSpheres = [ (4/3)*math.pi*3**3, (4/3)*math.pi*( (6**3)-(3**3) )]
		totalVolume = self.cellA.volume
		expNormFactor = (totalVolume/nRelevantAtomsInCell) * (1/nRelevantAtomsInCell) #First term related to AVERAGE volume per atom, second term is to account for the fact we're getting an AVERAGE over all the oxygen-oxygen rds
		expRdf = [(x/v)*expNormFactor for x,v in zip(expCounts,volumeSpheres)]
		expOutObj = tCode.RdfBinnedResultsSimple(expDists, expRdf, binEdges=expEdges, counts=expCounts)
		actOutObj = tCode.getSimpleEleEleRdf(self.trajObj, "O", "O", self.range, nBins=self.nBins)
		self.assertEqual(expOutObj, actOutObj)

	def testRaisesIfRangeIsTooLarge(self):
		""" Test we get an error if trying to calculate rdf at ranges > L/2 """
		self.range = [3,16] #Max values if greater than L/2
		with self.assertRaises(ValueError):
			outObj = tCode.getSimpleEleEleRdf(self.trajObj, "O", "O", self.range)


class TestRdfBinnedResultsClass(unittest.TestCase):

	def setUp(self):
		self.distsA = [1, 2, 3]
		self.binEdgesA = [0, 1.5, 2.5, 3.5]
		self.countsA = [2, 4, 6]
		self.rdf = [0.5, 1, 1.5]
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.RdfBinnedResultsSimple(self.distsA, self.rdf, binEdges=self.binEdgesA, counts=self.countsA)

	def testEqualObjsCompareEqual(self):
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = self.testObjA
		self.assertEqual(objA,objB)		
		self.assertEqual(objB,objA)		

	def testEqualObjsCompareEqual_countsSetToNone(self):
		self.countsA = None
		self.createTestObjs()
		objA = copy.deepcopy(self.testObjA)
		self.createTestObjs()
		objB = self.testObjA
		self.assertEqual(objA, objB)

	def testUnequalObjsCompareUnequal_diffRdfVals(self):
		objA = copy.deepcopy(self.testObjA)
		self.rdf[-1]+=1
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)
		self.assertNotEqual(objB, objA)

	def testUnequalObjsCompareUnequal_diffLengthDists(self):
		objA = copy.deepcopy(self.testObjA)
		self.distsA.append(4)
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)
		self.assertNotEqual(objB, objA)

	def testUnequalObjsCompareUnequal_diffCounts(self):
		objA = copy.deepcopy(self.testObjA)
		self.countsA[-1] += 1
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)
		self.assertNotEqual(objB, objA)

	def testUnequalObjsCompareUnequal_countsNoneOnOne(self):
		objA = copy.deepcopy(self.testObjA)
		self.assertTrue(objA.counts is not None)
		self.countsA = None
		self.createTestObjs()
		objB = self.testObjA
		self.assertNotEqual(objA, objB)
		self.assertNotEqual(objB, objA)


