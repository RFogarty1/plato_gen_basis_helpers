
import copy
import itertools as it
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.binned_res as binResHelp
import gen_basis_helpers.analyse_md.traj_core as trajCoreHelp
import gen_basis_helpers.analyse_md.water_rotations as tCode


class TestStandardWaterDistribOpts(unittest.TestCase):

	def setUp(self):
		self.binEdges = [0,50,90]
		self.checkEdges = True
		self.waterIndices = [ [1,2,3], [4,5,6] ]
		self.angleType = "roll"
		self.createTestObjs()

	def createTestObjs(self):
		self.binObjA = binResHelp.BinnedResultsStandard.fromBinEdges(self.binEdges)
		args = [self.binObjA, self.waterIndices]
		kwargs = {"angleType":self.angleType, "checkEdges":self.checkEdges}
		self.testObjA = tCode.CalcStandardWaterOrientationDistribOptions(*args,**kwargs)

	def testExpectedDomainsUponChanging(self):
		angleTypes = ["roll", "pitch", "azimuth"]
		expDomains = [ [-90,90], [-90,90], [-180,180] ]
		for aType,expDomain in it.zip_longest(angleTypes,expDomains):
			self.testObjA.angleType = aType
			actDomain = self.testObjA.domain
			self.assertAlmostEqual(expDomain[0], actDomain[0])
			self.assertAlmostEqual(expDomain[1], actDomain[1])

	def testRaisesIfCheckEdgesAndWeSetToAngleTypeWithLowerDomain(self):
		#Check error on initiation
		self.binEdges = [0,50,92]
		with self.assertRaises(ValueError):
			self.createTestObjs()

		#Check error when changin from azimuth (large domain) to pitch (smaller domain)
		self.angleType="azimuth"
		self.createTestObjs()
		with self.assertRaises(ValueError):
			self.testObjA.angleType = "roll"


class TestGetWaterStandardRotationsFromInpCell(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		#First water should be (0,0,90) while second should be in a roll=90 (90,0,180)
		self.coordsA = [ [0,0,4,"X"],
		                 [2,2,2,"O"],
		                 [3,4,2,"H"],
		                 [1,4,2,"H"],
		                 [7,9,7,"H"],
		                 [9,9,9,"O"],
		                 [7,9,1,"H"],
		                 [3,4,5,"X"] ] #Min image conv probably needed here

		self.waterIndices = [ [1,2,3], [6,5,4] ] #Internal ordering shouldnt matter
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coordsA 

	def _runTestFunct(self):
		return tCode.getWaterStandardRotationAnglesForInpCell(self.cellA, self.waterIndices)

	def _checkExpAndActAnglesEqual(self, expAngles, actAngles):
		for exp,act in it.zip_longest(expAngles, actAngles):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp,act)]

	def testExpectedRotationsA(self):
		expAngles = [ [0,0,90], [90,0,180] ]
		actAngles = self._runTestFunct()

		for exp,act in it.zip_longest(expAngles, actAngles):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp,act)]


	def testStandardOrientationWithDiffOrderings(self):
		""" Checking that roll is zero for standard orientation regardless of ordering; 180 degrees is another value it MIGHT have ended up taking """
		#All these should have roll=0, but could also have roll = 180
		waterA = [ [2,2,2,"O"],
		           [3,3,2,"H"],
		           [3,1,2,"H"] ]

		waterB = [ [2,2,2,"O"],
		           [3,1,2,"H"],
		           [3,3,2,"H"] ]

		waterC = [ [2,2,2,"O"],
		           [1,3,2,"H"],
		           [1,1,2,"H"] ]

		waterD = [ [2,2,2,"O"],
		           [1,1,2,"H"],
		           [1,3,2,"H"] ]

		self.coordsA = waterA + waterB + waterC + waterD

		self.waterIndices = [ [0,1,2], [3,4,5], [6,7,8], [9,10,11] ]
		self.createTestObjs()

		expAngles = [ [0,0,0], [0,0,0], [0,0,180], [0,0,180] ]
		actAngles = self._runTestFunct() 

		self._checkExpAndActAnglesEqual(expAngles,actAngles)


class TestPopulateOrientationAngularDistribs(unittest.TestCase):

	def setUp(self):
		#coords to use
		water_azi90   = [ [0,0,0,"O"], [-1,1,0,"H"], [1, 1, 0,"H"] ]
		water_azi120  = [ [0,0,0,"O"], [-1.37, 0.37, 0,"H"], [0.37,1.37,0,"H"] ]
		water_roll70  = [ [0,0,0,"O"], [1,0.34,0.94,"H"] , [1,-0.34,-0.94,"H"] ]

#		water_roll70  = [ [0,0,0,"O"], [1,0.34,0.94,"H"] , [1,-0.34,-0.94,"H"] ]
#		water_rollm70 = [ [0,0,0,"O"], [1,0.34,-0.94,"H"], [1,-0.34,0.94,"H"] ]
		water_roll70 = [ [0,0,0,"O"], [1,0.34,-0.94,"H"], [1,-0.34,0.94,"H"] ]



		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.cartCoords = water_azi90 + water_azi120 + water_roll70
		
		#bin options
		self.angleTypeA = "azimuth"
		self.angleTypeB = "roll"
		self.indicesA = [ [0,1,2],[3,4,5],[6,7,8] ]
		self.indicesB = [ [0,1,2],[3,4,5],[6,7,8] ]
		self.binEdgesA = [-10,80,180]
		self.binEdgesB = [-10,40,90]

		self.createTestObjs()

	def createTestObjs(self):
		#Create the trajectory
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords
		self.trajStepA = trajCoreHelp.TrajStepFlexible(unitCell=self.cellA)
		self.traj = trajCoreHelp.TrajectoryInMemory([self.trajStepA])

		#Create two options objs
		self.binsA = binResHelp.BinnedResultsStandard.fromBinEdges(self.binEdgesA)
		self.binsB = binResHelp.BinnedResultsStandard.fromBinEdges(self.binEdgesB)

		self.optsA = tCode.CalcStandardWaterOrientationDistribOptions(self.binsA, self.indicesA, angleType=self.angleTypeA)
		self.optsB = tCode.CalcStandardWaterOrientationDistribOptions(self.binsB, self.indicesB, angleType=self.angleTypeB)

	def _runTestFunct(self):
		tCode.populateWaterOrientationDistribsFromOptionsObjs(self.traj, [self.optsA, self.optsB])

	#Mainly trying to check adf values are as expected; rather than the counts
	def testExpectedValsSingleTrajStep(self):
		#Figure out what to expect
		expCountsA, expCountsB = [1,2], [2,1]
		expPdfA, expPdfB = [1/3, 2/3], [2/3, 1/3]
		expAdfA = [ (1/3)*(360/90), (2/3)*(360/100)]
		expAdfB = [ (2/3)*(180/50) , (1/3)*(180/50)]

		#Create the expected bins
		expBinA, expBinB = copy.deepcopy(self.binsA), copy.deepcopy(self.binsB)
		expBinA.binVals = {"adf":expAdfA, "pdf":expPdfA, "counts":expCountsA}
		expBinB.binVals = {"adf":expAdfB, "pdf":expPdfB, "counts":expCountsB}

		#Test and compare
		self._runTestFunct()

		self.assertEqual(expBinA, self.binsA)
		self.assertEqual(expBinB, self.binsB)

class TestOrientationMultiBinner(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		#NOTE: Not sure about the labelling on -70 and +70 roll
		#setting up water with an angle of 90 degrees and sqrt(2) bondlengths...because
		water_stdRot  = [ [0,0,0,"O"], [ 1,1,0,"H"], [1,-1, 0,"H"] ]
		water_azi90   = [ [0,0,0,"O"], [-1,1,0,"H"], [1, 1, 0,"H"] ]
		water_azi120  = [ [0,0,0,"O"], [-1.37, 0.37, 0,"H"], [0.37,1.37,0,"H"] ]
		water_rollm70  = [ [0,0,0,"O"], [1,0.34,0.94,"H"] , [1,-0.34,-0.94,"H"] ]
		water_roll70 = [ [0,0,0,"O"], [1,0.34,-0.94,"H"], [1,-0.34,0.94,"H"] ]

		self.cartCoords = ( [ [5,5,5,"X"], [8,8,8,"Y"] ] + water_stdRot + water_azi90 + water_azi120 + 
		                    water_roll70 + water_rollm70 )

		#Setup indices to probe
		self.indicesA = [ [5,6,7], [8,9,10], [11,12,13], ]
		self.indicesB = [  [11,12,13], [8,9,10], [2,3,4], [5,6,7],  [14,15,16] ]

		#Setup bin edges
		self.binEdgesA = [-0.1,89,180] #Azimuth domain is -180 to 180
		self.binEdgesB = [-90,40,90]

		#Angle types to look for
		self.angleIdxA = 2 #Azimuth
		self.angleIdxB = 0 #Roll

		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords
		self.trajStepA = trajCoreHelp.TrajStepFlexible(unitCell=self.cellA)

		#Create the bin objs
		self.binObjA = binResHelp.BinnedResultsStandard.fromBinEdges(self.binEdgesA)
		self.binObjB = binResHelp.BinnedResultsStandard.fromBinEdges(self.binEdgesB)

		#Create single binners
		self.binnerA = tCode._WaterRotationAngleBinnerFixedIndices(resBins=self.binObjA, indices=self.indicesA, angleIdx=self.angleIdxA)
		self.binnerB = tCode._WaterRotationAngleBinnerFixedIndices(resBins=self.binObjB, indices=self.indicesB, angleIdx=self.angleIdxB)

		#Create the multi binner (the test obj)
		self.testObj = tCode._WaterRotationAngleMultiBinnerFixedIndices([self.binnerA, self.binnerB])

	def _runTestFunct(self):
		self.testObj.updateCountsFromTrajStep(self.trajStepA)

	def testExpectedCountsTwoBinObjs(self):
		expBinA, expBinB = copy.deepcopy(self.binObjA), copy.deepcopy(self.binObjB)
		expBinA.binVals["counts"], expBinB.binVals["counts"] = [1,2], [4,1]
		self._runTestFunct()
		actBinA, actBinB = self.binObjA, self.binObjB
		self.assertEqual(expBinA, actBinA)
		self.assertEqual(expBinB, actBinB)


