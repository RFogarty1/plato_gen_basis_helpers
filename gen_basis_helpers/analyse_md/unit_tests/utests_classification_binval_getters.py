

import unittest

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.atom_combo_opts_obj_maps as optObjMaps
import gen_basis_helpers.analyse_md.classification_distr_opt_objs as classDistrOptObjHelp
import gen_basis_helpers.analyse_md.binned_res as binResHelp



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



