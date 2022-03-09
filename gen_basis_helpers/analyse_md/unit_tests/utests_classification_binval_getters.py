
import copy
import unittest

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.atom_combo_opts_obj_maps as optObjMaps
import gen_basis_helpers.analyse_md.classification_distr_opt_objs as classDistrOptObjHelp
import gen_basis_helpers.analyse_md.binned_res as binResHelp
import gen_basis_helpers.analyse_md.filtered_atom_combo_obj_maps as filteredObjMapHelp


class TestCountAtomClassify(unittest.TestCase):

	def setUp(self):
		# Geometry
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.coords = [  [0,0,1,"X"],
		                 [0,0,2,"Y"],
		                 [0,0,3,"Y"],
		                 [0,0,4,"Y"],
		                 [0,0,5,"Y"] ]

		# Options
		binEdges = [0,1,2,3,4,5,6]
		self.binResObjs = [ binResHelp.BinnedResultsStandard.fromBinEdges(binEdges) for x in range(2) ]
		self.atomIndices = [1,2,3,4]
		self.distFilterIndices = [0]
		self.distFilterRanges = [ [0,3.5], [3.5,10] ]

		self.createTestObjs()

	def createTestObjs(self):
		#Geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams,lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Options
		currArgs = [self.binResObjs, self.atomIndices, self.distFilterIndices, self.distFilterRanges]
		self.classifierOpts = classDistrOptObjHelp.AtomClassifyBasedOnDistsFromIndicesSimpleOpts(*currArgs)

		#Matrix populator
		self.sparseMatrixCalculator = optObjMaps.getSparseMatrixCalculatorFromOptsObjIter([self.classifierOpts])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

		#Create the binner object
		self.testObj = optObjMaps.getMultiDimBinValGetterFromOptsObjs([self.classifierOpts])

	def _runTestFunct(self):
		return self.testObj.getValsToBin(self.sparseMatrixCalculator)

	def testExpectedA(self):
		expVals = [ (3,1) ]
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)


class TestWaterCountBasedOnAdsSiteHozDists(unittest.TestCase):

	def setUp(self):
		#1) All geometric parameters for testing
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		#roll 90,pitch=45, OH len~1, HOH angle 104.5. Then just added translation vectors
		self.waterACoords = [ [0,0,0,"O"], [-0.13,0,0.99,"H"], [0.99,0,-0.13,"H"] ]
		#translation = [2,0,0]; no rotations
		self.waterBCoords = [ [2,0,0,"O"], [2.61, 0.79, 0, "H"], [2.61,-0.79,0,"H"] ]
		#translation = [2+(2*0.61), 2*0.79, 0]
		self.waterCCoords = [ [3.22, 1.58, 0, "O"],  [3.83, 2.37, 0, "H"], [3.83, 0.79, 0, "H"] ]
		#translation = [2+(2*0.61), -2*0.79,0]
		self.waterDCoords = [ [3.22, -1.58, 0, "O"], [3.83, -0.79, 0, "H"], [3.83, -2.37, 0, "H"] ]

		self.adsCoords = [ [0.1,0,0.1,"X"], [2.1,0,0.1,"X"], [ 3,0,0.1,"X" ] ]
		self.coords = self.waterACoords + self.waterBCoords + self.waterCCoords + self.waterDCoords + self.adsCoords


		#Options
		self.oxyIndices = [0,3,6,9]
		self.hyIndices = [ [1,2], [4,5], [7,8], [10,11] ]
		self.distFilterIndices = [12,13,14]
		self.distFilterVals =   [ [0,1]    , [0,1] ] #These dist filter ranges mean we ignore the case of a shared adsorption site
		self.adsHozDistRanges = [ [0,2.5], [4,10] ] 
		self.binResObjs = [None,None] #Irrelevent to these tests so....


		self.createTestObjs()

	def createTestObjs(self):
		#Sort the geometry out
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams,lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Create an options object
		currArgs = [ self.binResObjs, self.oxyIndices, self.hyIndices, self.distFilterIndices, self.distFilterVals ]
		currKwargs = {"adsSiteMinHozToOtherAdsSiteRanges": self.adsHozDistRanges}
		self.optObj = classDistrOptObjHelp.WaterAdsorbedClassifier_usingMinHozDistsBetweenAdsorptionSitesOptsObj(*currArgs, **currKwargs)

		#Get sparse matrix populator + populate it
		self.sparseMatrixCalculator = optObjMaps.getSparseMatrixCalculatorFromOptsObjIter([self.optObj])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

		#Get the binner object
		self.testObj = optObjMaps.getMultiDimBinValGetterFromOptsObjs([self.optObj])

	def testExpectedA(self):
		expBinVals = [ (2,0) ]
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)
		self.assertEqual(expBinVals, actBinVals)

		#Calling it twice; second time it should be using various cached values
		#(So it COULD theoretically fail second time for certain types of cache-error)
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)
		self.assertEqual(expBinVals, actBinVals)

	#Need a special check for shared adsorption sites;
	def testExpected_sharedAdsSiteImportant(self):
		""" Shouldnt automatically be zero distance """
		self.distFilterVals = [ [0,3], [0,3] ]
		self.adsHozDistRanges = [ [0,1], [1,2.5] ]
		self.createTestObjs()
		expBinVals = [ (3,1) ] #1st ads site is 2 away from nearest; others are each 0.9 away
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)
		self.assertEqual(expBinVals, actBinVals)

	def testExpected_hozDistRequriedToFilterOut(self):
		self.distFilterVals = [ [0,3] ]
		self.adsHozDistRanges = [ [1,2.5] ]
		self.createTestObjs()
		expBinVals = [(1,)]
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)
		self.assertEqual(expBinVals, actBinVals)


class TestWaterCountTypesMinDistAndHBond(unittest.TestCase):

	#Taken from TestDiscHBondCounterBetweenGroupsOxyDistFilter mostly
	def setUp(self):
		#1) All geometric parameters for testing
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		#roll 90,pitch=45, OH len~1, HOH angle 104.5. Then just added translation vectors
		self.waterACoords = [ [0,0,0,"O"], [-0.13,0,0.99,"H"], [0.99,0,-0.13,"H"] ]
		#translation = [2,0,0]; no rotations
		self.waterBCoords = [ [2,0,0,"O"], [2.61, 0.79, 0, "H"], [2.61,-0.79,0,"H"] ]
		#translation = [2+(2*0.61), 2*0.79, 0]
		self.waterCCoords = [ [3.22, 1.58, 0, "O"],  [3.83, 2.37, 0, "H"], [3.83, 0.79, 0, "H"] ]
		#translation = [2+(2*0.61), -2*0.79,0]
		self.waterDCoords = [ [3.22, -1.58, 0, "O"], [3.83, -0.79, 0, "H"], [3.83, -2.37, 0, "H"] ]

		self.xCoord = [[0,0,0,"X"]]
		self.coords = self.waterACoords + self.waterBCoords + self.waterCCoords + self.waterDCoords + self.xCoord

		#Options
		self.oxyIndices = [0,3,6,9]
		self.hyIndices = [ [1,2], [4,5], [7,8], [10,11] ]
		self.distFilterIndices = [12]
		self.distFilterVals = [ [0,3], [3,5] ]#AB should be one group, with CD as the other. + easy to flip this
		self.maxOO = 3 #AC h-bond would be possible iff this was set high enough i suspect
		self.maxAngle = 35

		self.binResObjs = [None,None] #Irrelevent to these tests so....

		self.nDonorFilterRanges = [ [-0.5,1.5], [0.5,2.5] ]
		self.nAcceptorFilterRanges = [ [-0.5,0.5], [0.5,2.5] ]
		self.nTotalFilterRanges = [ [-0.5,1.5], [0.5,2.5] ]
		self.checkInputConsistent = True

		self.createTestObjs()

	def createTestObjs(self):
		#Sort the geometry out
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams,lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Create an options object
		currArgs = [self.binResObjs, self.oxyIndices, self.hyIndices, self.distFilterIndices, self.distFilterVals]
		currKwargs = {"nDonorFilterRanges": self.nDonorFilterRanges, "nAcceptorFilterRanges": self.nAcceptorFilterRanges,
		              "nTotalFilterRanges": self.nTotalFilterRanges, "maxOOHBond": self.maxOO, "maxAngleHBond": self.maxAngle,
		              "checkInputConsistent": self.checkInputConsistent}
		self.optObj = classDistrOptObjHelp.WaterCountTypesMinDistAndHBondSimpleOpts(*currArgs,**currKwargs)

		#Get sparse matrix populator + populate it
		self.sparseMatrixCalculator = optObjMaps.getSparseMatrixCalculatorFromOptsObjIter([self.optObj])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

		#Create the binner object
		self.testObj = optObjMaps.getMultiDimBinValGetterFromOptsObjs([self.optObj])

	def testExpectedA(self):
		expBinVals = [ (1,0) ]
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)
		self.assertEqual(expBinVals, actBinVals)

	def testExpectedB_diffRanges(self):
		#Put first 2 in first bin, second 2 in 2nd bin
		self.nDonorFilterRanges    =  [  [-0.5,3.5], [-0.5,0.5] ]
		self.nAcceptorFilterRanges =  [  [-0.5,1.5], [0.5,1.5 ] ]
		self.nTotalFilterRanges    =  [  [0.5,3.5 ], [0.5,1.5 ] ]
		self.createTestObjs()

		expBinVals = [ (2,2) ]
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)
		self.assertEqual(expBinVals,actBinVals)


class TestChainedClassifyUsingMinHozDists(unittest.TestCase):

	def setUp(self):
		#
		self.lattParams, self.lattAngles = [13,13,13], [90,90,90]

		self.hydroxylA = [ [0,0,0,"O"], [1,0,0,"H"] ]
		self.hydroxylB = [ [1,0,1,"O"], [2,0,1,"H"] ]
		self.hydroxylC = [ [2,0,2,"O"], [3,0,2,"H"] ]
		self.hydroxylD = [ [3,0,3,"O"], [4,0,3,"H"] ]
		self.hydroxylE = [ [4,0,4,"O"], [5,0,4,"H"] ]

		self.cartCoords = self.hydroxylA + self.hydroxylB + self.hydroxylC + self.hydroxylD + self.hydroxylE

		#Options for classifiers
		self.binResObjs = [None]
		self.fromNonHyA, self.fromNonHyB = [ [0], [2] ], [ [0], [2] ]
		self.fromHyA, self.fromHyB       = [ [1], [3] ], [ [1], [3] ]
		self.toNonHyA, self.toNonHyB = [ [4], [6], [8] ], [ [4], [6], [8] ]
		self.toHyA, self.toHyB       = [ [5], [7], [9] ], [ [5], [7], [9] ]

		#
		self.minHozDistRangesA, self.minHozDistRangesB = [ [-0.1,2.5] ], [ [0.5,2.5] ] #Essentially only the narrower one matters here
		self.useIndicesFromA, self.useIndicesFromB = "all", "all"
		self.useIndicesToA, self.useIndicesToB = "all", "all"
		self.minDistValA, self.minDistValB = -0.01, -0.01

		self.createTestObjs()

	def createTestObjs(self):
		#Geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

		#Create both options objects
		currArgs = [self.binResObjs, self.fromNonHyA, self.fromHyA, self.toNonHyA, self.toHyA, self.minHozDistRangesA]
		currKwargs = {"useIndicesFrom":self.useIndicesFromA, "useIndicesTo":self.useIndicesToA, "minDistVal":self.minDistValA}
		self.optObjA = classDistrOptObjHelp.ClassifyNonHyAndHyBasedOnMinHozDistsToAtomGroups(*currArgs, **currKwargs)

		currArgs = [self.binResObjs, self.fromNonHyB, self.fromHyB, self.toNonHyB, self.toHyB, self.minHozDistRangesB]
		currKwargs = {"useIndicesFrom":self.useIndicesFromB, "useIndicesTo":self.useIndicesToA, "minDistVal":self.minDistValB}
		self.optObjB = classDistrOptObjHelp.ClassifyNonHyAndHyBasedOnMinHozDistsToAtomGroups(*currArgs, **currKwargs)

		#
		self.optObjCombo = classDistrOptObjHelp.ClassifyNonHyAndHyChainedAllCommon([self.optObjA,self.optObjB])

		#Get the sparse matrix calculator and populate it
		self.sparseMatrixCalculator = optObjMaps.getSparseMatrixCalculatorFromOptsObjIter([self.optObjCombo])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

		#Create the binner object
		self.testObj = optObjMaps.getMultiDimBinValGetterFromOptsObjs([self.optObjCombo])

	def testExpectedCaseA(self):
		expBinVals = [(1,)]
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)
		self.assertEqual(expBinVals, actBinVals)


class TestNonHyClassifyByMinHozDist(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [13,13,13], [90,90,90]

		self.hydroxylA = [ [0,0,0,"O"], [1,0,0,"H"] ]
		self.hydroxylB = [ [1,0,1,"O"], [2,0,1,"H"] ]
		self.hydroxylC = [ [2,0,2,"O"], [3,0,2,"H"] ]
		self.hydroxylD = [ [3,0,3,"O"], [4,0,3,"H"] ]
		self.hydroxylE = [ [4,0,4,"O"], [5,0,4,"H"] ]

		self.cartCoords = self.hydroxylA + self.hydroxylB + self.hydroxylC + self.hydroxylD + self.hydroxylE

		#Options for the classifier
		self.binResObjs = [ None, None] #Irrelevant for these tests so....
		self.fromNonHyIndices = [ [0], [2] ]
		self.fromHyIndices    = [ [1], [3] ]
		self.toNonHyIndices   = [ [4], [6], [8] ]
		self.toHyIndices      = [ [5], [7], [9] ] 
		self.minHozDistRanges = [ [-0.1,1.5], [1.5,7.5] ] #Need to use <0 at lower end for edge case of exactly zero
		self.useIndicesFrom = "all"
		self.useIndicesTo = "all"
		self.minDistVal = -0.01

		self.createTestObjs()

	def createTestObjs(self):
		#Geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

		#Create an options object
		currArgs = [self.binResObjs, self.fromNonHyIndices, self.fromHyIndices, self.toNonHyIndices, self.toHyIndices,
		            self.minHozDistRanges]
		currKwargs = {"useIndicesFrom":self.useIndicesFrom, "useIndicesTo":self.useIndicesTo, "minDistVal":self.minDistVal}
		self.optObj = classDistrOptObjHelp.ClassifyNonHyAndHyBasedOnMinHozDistsToAtomGroups(*currArgs, **currKwargs)

		#Get sparse matrix populator + populate it
		self.sparseMatrixCalculator = optObjMaps.getSparseMatrixCalculatorFromOptsObjIter([self.optObj])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

		#Create the binner object
		self.testObj = optObjMaps.getMultiDimBinValGetterFromOptsObjs([self.optObj])

	def testExpected_fromAllToAll(self):
		""" Test we get expected when using all indices on both source and target molecules """
		expBinVals = [ (2,0) ]
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)
		self.assertEqual(expBinVals, actBinVals)

	def testExpected_fromOxyToAll(self):
		self.useIndicesFrom = "nonHy"
		self.createTestObjs()

		expBinVals = [ (1,1) ]
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)
		self.assertEqual(expBinVals, actBinVals)

	def testExpected_fromOxyToHy(self):
		self.useIndicesFrom = "nonHy"
		self.useIndicesTo = "hy"
		self.createTestObjs()

		expBinVals = [ (0,2) ]
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)
		self.assertEqual(expBinVals, actBinVals)

	def testExpected_fromOxyToNonHy(self):
		self.useIndicesFrom, self.useIndicesTo = "nonHy", "nonHy"
		self.createTestObjs()

		expBinVals = [ (1,1) ]
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)
		self.assertEqual(expBinVals, actBinVals)

	def testExpected_fromHyToAll(self):
		self.useIndicesFrom, self.useIndicesTo = "hy", "all"
		self.createTestObjs()

		expBinVals = [ (2,0) ]
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)
		self.assertEqual(expBinVals, actBinVals)

	def testExpected_reverseOrder(self):
		self.useIndicesFrom = "nonHy"
		self.fromNonHyIndices = [x for x in reversed([ [0], [2] ])]
		self.fromHyIndices    = [x for x in reversed([ [1], [3] ])]
		self.createTestObjs()

		expBinVals = [ (1,1) ]
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)
		self.assertEqual(expBinVals, actBinVals)

	def testExpected_overlappingGroups_minDistValOff(self):
		self.useIndicesFrom, self.indicesTo = "nonHy", "nonHy"
		self.fromNonHyIndices = [ [0], [2] ]
		self.fromHyIndices    = [ [1], [3] ]
		self.toNonHyIndices   = [ [0], [4], [6], [8] ]
		self.toHyIndices      = [ [1], [5], [7], [9] ] 
		self.createTestObjs()

		expBinVals = [ (2,0) ]
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)
		self.assertEqual(expBinVals, actBinVals)

	def testExpected_overlappingGroups_minDistValOn(self):
		self.useIndicesFrom, self.useIndicesTo = "nonHy", "nonHy"
		self.fromNonHyIndices = [ [0], [2] ]
		self.fromHyIndices    = [ [1], [3] ]
		self.toNonHyIndices   = [ [0], [4], [6], [8] ]
		self.toHyIndices      = [ [1], [5], [7], [9] ] 
		self.minDistVal = 0.02
		self.createTestObjs()

		expBinVals = [ (1,1) ]
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)
		self.assertEqual(expBinVals, actBinVals)


class TestClassifyByHBondsToGenericGroup(unittest.TestCase):

	def setUp(self):
		#The geometry(left OH donates 1-hbond to water; right OH accepts 1)
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		self.hydroxylA = [ [0,0,0,"O"], [1,0,0,"H"] ]
		self.waterA = [ [2,1,0,"O"], [3,0,0,"H"], [3,2,0,"H"] ]
		self.hydroxylB = [ [4,0,0,"O"], [5,0,0,"H"] ]
		self.cartCoords = self.hydroxylA + self.waterA + self.hydroxylB

		#Options for classifier
		self.binResObjs = [None,None] #Irrelevent to these tests so....
		self.nonHyFromIndices = [ [2] ] #Water
		self.hyFromIndices = [ [3,4] ]
		self.nonHyToIndices = [ [0], [5] ]
		self.hyToIndices = [ [1], [6] ]
		self.maxOOHBond = 4 #Mainly the angle thats important in this geometry
		self.maxAngleHBond = 45

		#Filter options for classifier
		self.nDonorFilterRanges = [ [0.5,1.5], [-0.5,0.5] ]
		self.nAcceptorFilterRanges = [ [0.5,1.5], [-0.5,2.5] ]
		self.nTotalFilterRanges = [ [0.5,4.5], [0.5,4.5] ] #Basically unused

		self.createTestObjs()

	def createTestObjs(self):
		#Geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

		#Create an options object
		currArgs = [self.binResObjs, self.nonHyFromIndices, self.hyFromIndices, self.nonHyToIndices, self.hyToIndices]
		currKwargs = {"nDonorFilterRanges":self.nDonorFilterRanges, "nAcceptorFilterRanges":self.nAcceptorFilterRanges,
		              "nTotalFilterRanges":self.nTotalFilterRanges, "maxOOHBond":self.maxOOHBond, "maxAngleHBond":self.maxAngleHBond}

		self.optObj = classDistrOptObjHelp.ClassifyBasedOnHBondingToGroup_simple(*currArgs, **currKwargs)

		#Get sparse matrix populator + populate it
		self.sparseMatrixCalculator = optObjMaps.getSparseMatrixCalculatorFromOptsObjIter([self.optObj])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

		#Create the binner object
		self.testObj = optObjMaps.getMultiDimBinValGetterFromOptsObjs([self.optObj])

	def testExpectedCaseA_filterWater(self):
		expBinVals = [ (1,0) ] #The water molecule goes into the first bin since it donates/accepts 1 hbond from our group
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)
		self.assertEqual(expBinVals,actBinVals)

	def testExpectedCaseB_filterHydroxylA(self):
		#Sort options
		self.nonHyFromIndices = [ [0], [5] ]
		self.hyFromIndices = [ [1], [6] ]
		self.nonHyToIndices = [ [2] ]
		self.hyToIndices = [ [3,4] ]

		#
		self.nDonorFilterRanges = [ [0.5,1.5], [-0.5,0.5] ] #[1 donor, 0 donor] 
		self.nAcceptorFilterRanges = [ [-0.5,0.5], [-0.5,2.5] ] # [0 acceptor, 0/1/2 acceptor] 
		self.nTotalFilterRanges = [ [-0.5,4.5], [-0.5,4.5] ] #Basically unused
		self.createTestObjs()

		#Run + test
		expBinVals = [ (1,1) ]
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)
		self.assertEqual(expBinVals, actBinVals)

class TestClassifyByHBondingToDynamicGroup(unittest.TestCase):

	def setUp(self):
		#The geometry(left OH donates 1-hbond to middle water; middle water donates one to right-hand water)
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.hydroxylA =  [ [2,2,2,"O"], [3,2,2,"H"] ]
		self.waterA    = [ [4,2,2,"O"], [5,2,2,"H"], [3.5,1.5,2,"H"] ]
		self.waterB    = [ [6,2,2,"O"], [7,2,2,"H"], [5.5,1.5,2,"H"] ]

		self.cartCoords = self.hydroxylA + self.waterA + self.waterB

		#Options for the classifiers [some options shared]
		self.binResObjs = [None, None] #Irrelevant to these tests so...
		self.fromNonHyIndices, self.fromHyIndices = [ [2], [5] ], [ [3,4], [6,7]  ] #Used in both the objects
		self.toNonHyIndices, self.toHyIndices = [ [0] ], [ [1] ]
		self.maxOOHBond = 2.1
		self.maxAngleHBond = 30 #Angle doesnt matter in original 3-molecule geom
		self.nTotalFilterRangesStatic  = [ [-1,0.1], [-1 ,0.1] ] 
		self.nTotalFilterRangesDynamic = [ [-1,0.1], [0.1,1.1] ] #Count those with 0 h-bonds and those with 1 h-bond

		self.dynToNonHyIndices = copy.deepcopy(self.fromNonHyIndices)
		self.dynToHyIndices = copy.deepcopy(self.fromHyIndices)

		#Any extra options for the combined classifier options object
		self.mutuallyExclusive = True
		self.firstClassifierObjs = None

		self.createTestObjs()

	def createTestObjs(self):
		#Geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

		#Create the first classifier options object [For the group with static indices]
		currArgs = [self.binResObjs, self.fromNonHyIndices, self.fromHyIndices, self.toNonHyIndices, self.toHyIndices]
		currKwargs = {"nTotalFilterRanges":self.nTotalFilterRangesStatic, "maxOOHBond":self.maxOOHBond,
		              "maxAngleHBond":self.maxAngleHBond}
		self.staticGroupOptObj = classDistrOptObjHelp.ClassifyBasedOnHBondingToGroup_simple(*currArgs, **currKwargs)

		#Create the second classifier options object
		#Note we pass the .fromIndices of classifier 1 as "toIndices" here; since they ARE our target
		currArgs = [self.binResObjs, self.fromNonHyIndices, self.fromHyIndices, self.dynToNonHyIndices, self.dynToHyIndices]
		currKwargs = {"nTotalFilterRanges":self.nTotalFilterRangesDynamic, "maxOOHBond":self.maxOOHBond,
		              "maxAngleHBond":self.maxAngleHBond}
		self.dynamicGroupOptObj = classDistrOptObjHelp.ClassifyBasedOnHBondingToGroup_simple(*currArgs, **currKwargs)


		#Create the overall classifier options
		currArgs = [self.staticGroupOptObj, self.dynamicGroupOptObj]
		currKwargs = {"mutuallyExclusive":self.mutuallyExclusive, "firstClassifierObjs":self.firstClassifierObjs}
		self.optObj = classDistrOptObjHelp.ClassifyBasedOnHBondingToDynamicGroup(*currArgs, **currKwargs)

		#Get sparse matrix populator + populate it
		self.sparseMatrixCalculator = optObjMaps.getSparseMatrixCalculatorFromOptsObjIter([self.optObj])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

		#Create the binner object
		self.testObj = optObjMaps.getMultiDimBinValGetterFromOptsObjs([self.optObj])

	def testExpectedCaseA(self):
		expBinVals = [ (0,1) ] #Second bin; accepts 1 h-bond from a species with 1 h-bond to the hydroxyl 
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)
		self.assertEqual(expBinVals,actBinVals)

	#First water has 0 h-bonds to first group (since its the ONLY member); hence goes in first bin
	#Second water has 1 h-bond to the first group
	def testExpectedNotMutuallyExclusive(self):
		self.mutuallyExclusive = False
		self.createTestObjs()
		expBinVals = [ (1,1) ]
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)
		self.assertEqual(expBinVals, actBinVals)

	def testWhenClassifierOptObjPassed(self):
		self.firstClassifierObjs = [filteredObjMapHelp.getClassifiersFromOptsObj(self.staticGroupOptObj)[0] for x in range(len(self.nTotalFilterRangesStatic))]
		self.createTestObjs()

		expBinVals = [ (0,1) ]
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)
	
		self.assertEqual(expBinVals, actBinVals)
		for classifierObj in self.firstClassifierObjs:
			self.assertEqual( classifierObj.execCount, 1 )

	def testWhenDynamicallyAssignedGroupEmpty(self):
		self.maxOOHBond = 0.1
		self.mutuallyExclusive = False
		self.createTestObjs()

		expBinVals = [ (2,0) ] #Both have 0 h-bonds; means first bin in both cases. And since we allow them to be in both groups its [(2,0)]
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)

		self.assertEqual(expBinVals, actBinVals)

	def testRaisesForInconsitentIndices(self):
		""" Slightly flawed; could have failed for breaking EITHER 1 of 2 (so only half tests what i want really) """
		self.dynToNonHyIndices[0][0] += 1
		self.dynToHyIndices[0][0] += 1
		with self.assertRaises(ValueError):
			self.createTestObjs()

	def testRaisesForInconsistentFilterRanges(self):
		""" Different length filter indices should raise an error; i only test ONE type explicitly though """
		self.nTotalFilterRangesDynamic.append( [0,2] )
		with self.assertRaises(ValueError):
			self.createTestObjs()


class TestWaterDerivDistanceBasedCounter(unittest.TestCase):

	def setUp(self):
		#Geometry
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		self.water = [ [0,0,0,"O"], [0.5,0.5,0,"H"], [0.5,-0.5,0,"H"] ]
		self.hydroxylA = [ [5,5,0,"O"], [5,6,0,"H"] ]
		self.hydroxylB = [ [6,0,0,"O"], [7,0,0,"H"] ]
		self.cartCoords = self.water + self.hydroxylA + self.hydroxylB

		#Options for classifier object
		self.binResObjs = [None,None] #Irrelevent to these tests so....
		self.oxyIndices = [0,3,5]
		self.hyIndices = [1,2,4,6]
		self.maxOHDist = 1.5
		self.nNebs = [2,1]

		self.createTestObjs()

	def createTestObjs(self):
		#Geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

		#Create an options object
		currArgs = [self.binResObjs,self.oxyIndices, self.hyIndices]
		currKwargs = {"maxOHDist":self.maxOHDist, "nNebs":self.nNebs}
		self.optObj = classDistrOptObjHelp.WaterDerivativeBasedOnDistanceClassifierOptsObj(*currArgs, **currKwargs)
		
		#Get sparse matrix populator + populate it
		self.sparseMatrixCalculator = optObjMaps.getSparseMatrixCalculatorFromOptsObjIter([self.optObj])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

		#Create the binner object
		self.testObj = optObjMaps.getMultiDimBinValGetterFromOptsObjs([self.optObj])


	def testExpected_waterAndHydroxyl(self):
		expBinVals = [(1,2)]
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)
		self.assertEqual(expBinVals, actBinVals)

class TestClassifyByNumberNebsWithinDistance(unittest.TestCase):

	def setUp(self):
		#Geometry
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.water = [ [0,0,0,"O"], [0.5,0.5,0,"H"], [0.5,-0.5,0,"H"] ]
		self.freeH = [ [5,5,5,"H"] ]

		self.cartCoords = self.water + self.freeH

		#Options for classifier object
		self.binResObjs = [None,None] #Irrelevent to these tests so....
		self.fromIndices = [1,2,3]
		self.toIndices = [0,1]
		self.minDist = 0.1
		self.maxDist = 1
		self.nNebs = [ [-0.1,0.1],[0.9,1.1] ]

		self.createTestObjs()

	def createTestObjs(self):
		#Geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

		#Create an options object
		currArgs = [self.binResObjs,self.fromIndices, self.toIndices]
		currKwargs = {"maxDist":self.maxDist, "minDist":self.minDist,"nebRanges":self.nNebs}
		self.optObj = classDistrOptObjHelp.ClassifyByNumberNebsWithinDistanceOptsObj(*currArgs, **currKwargs)
		
		#Get sparse matrix populator + populate it
		self.sparseMatrixCalculator = optObjMaps.getSparseMatrixCalculatorFromOptsObjIter([self.optObj])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

		#Create the binner object
		self.testObj = optObjMaps.getMultiDimBinValGetterFromOptsObjs([self.optObj])

	def testExpectedGeomA(self):
		expBinVals = [(1,2)] #1 has zero neighbours; 2 have one neighbour
		actBinVals = self.testObj.getValsToBin(self.sparseMatrixCalculator)
		self.assertEqual(expBinVals, actBinVals)


