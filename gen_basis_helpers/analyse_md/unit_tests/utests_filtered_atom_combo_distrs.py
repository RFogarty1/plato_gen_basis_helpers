
import copy
import math
import itertools as it
import unittest

import numpy as np

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.binned_res as binResHelp
import gen_basis_helpers.analyse_md.classification_distr_opt_objs as clsDistrOptObjs
import gen_basis_helpers.analyse_md.distr_opt_objs as distrOptObjHelp
import gen_basis_helpers.analyse_md.filtered_atom_combo_opt_objs as filteredComboOptObjHelp
import gen_basis_helpers.analyse_md.filtered_atom_combo_obj_maps as filteredAtomComboObjMaps
import gen_basis_helpers.analyse_md.atom_combo_opts_obj_maps as optsObjMapHelp
import gen_basis_helpers.analyse_md.classifier_objs as classifierObjsHelp
import gen_basis_helpers.analyse_md.water_derivative_species_classifiers as waterDerivSpeciesClassifierHelp

import gen_basis_helpers.shared.plane_equations as planeEqnHelp


class TestGetBinValsAtomFilteredVariousDistribs(unittest.TestCase):

	def setUp(self):
		#Geometry
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.coords = [  [0,0,1,"X"],
		                 [0,0,2,"A"],
		                 [0,0,3,"B"],
		                 [0,0,4,"C"],
		                 [0,0,5,"D"] ]

		#General options
		self.classBinResObj = binResHelp.getEmptyBinResultsFromMinMaxAndWidthStandard(-0.1,20,1,extremesAtCentre=False)
		self.fromIndices = [1,2,3,4]

		#Overall options object specific
		self.useGroups = [ [0] ]

		#Classification Options object specific
		self.distFilterIndices = [0]
		self.distFilterRanges = [ [-0.1,3.5], [3.5,5.5] ] #Gets [1,2,3], [4] by default

		#Distribution opts obj specific
		self.distrBinResObj = binResHelp.BinnedResultsStandard.fromBinEdges([-0.1,20])

		#Create the standard distr opts here
		dudIndices = [20] #Pick a stupid value; since it should NOT be used anyway
		self.distrOpts = [ distrOptObjHelp.CalcPlanarDistOptions(self.distrBinResObj, dudIndices) ]

		self.createTestObjs()


	def createTestObjs(self):
		#Geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams,lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Create the classification options
		currBins = [self.classBinResObj for x in self.distFilterRanges]
		currArgs = [currBins, self.fromIndices, self.distFilterIndices, self.distFilterRanges]
		self.classifierOpts = clsDistrOptObjs.AtomClassifyBasedOnDistsFromIndicesSimpleOpts(*currArgs)

		#
		currArgs = [self.fromIndices, self.classifierOpts, self.distrOpts, self.useGroups]
		self.filterOptObj = filteredComboOptObjHelp.FilteredAtomComboOptsObjGeneric(*currArgs)

		#create the sparse matrix calculator + populate it
		self.sparseMatrixCalculator = optsObjMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.filterOptObj])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

		#Create the binval getter
		self.binValGetter = optsObjMapHelp.getMultiDimBinValGetterFromOptsObjs([self.filterOptObj])

	def _runTestFunct(self):
		return self.binValGetter.getValsToBin(self.sparseMatrixCalculator)

	def testPlanarDistrA(self):
		expVals = [ (2,), (3,), (4,) ]
		actVals = self._runTestFunct()
		for expIter,actIter in it.zip_longest(expVals, actVals):
			[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expIter,actIter)]



class TestGetBinValsAtomFilteredVaryClassifier(unittest.TestCase):

	def setUp(self):
		#Geometry
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		waterA = [ [0.0,0.0,0.0,"O"], [0.5,0.5,0.1,"H"], [0.5,-0.5,-0.2,"H"] ]
		freeH = [ [6,6,6,"H"] ] 
		self.coords = waterA + freeH

		#The default classification options
		unusedBinObjs = [ None ]
		self.fromIndices, self.toIndices = [1,2,3,], [0,1,2,3]
		nebRanges = [ [0.9,2.5] ]
		minDist, maxDist = 0.1, 1
		currArgs = [unusedBinObjs, self.fromIndices, self.toIndices]
		currKwargs = {"minDist":minDist, "maxDist":maxDist, "nebRanges":nebRanges}
		self.classOptObj = clsDistrOptObjs.ClassifyByNumberNebsWithinDistanceOptsObj(*currArgs, **currKwargs)

		#Options for the overall classifier
		self.useGroups = [ [0] ]

		#Distribution options
		dudBinObj = None
		dudIndices = [20] #Pick a stupid value; since it should NOT be used anyway
		self.distrOpts = [ distrOptObjHelp.CalcPlanarDistOptions(dudBinObj, dudIndices) ]

		self.createTestObjs()

	def createTestObjs(self):
		#Geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Create the filter option object
		currArgs = [self.fromIndices, self.classOptObj, self.distrOpts, self.useGroups]
		self.filterOptObj = filteredComboOptObjHelp.FilteredAtomComboOptsObjGeneric(*currArgs)

		#Create the sparse matrix calcualtor and populate it
		self.sparseMatrixCalculator = optsObjMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.filterOptObj])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

		#Create the binval getter
		self.binValGetter = optsObjMapHelp.getMultiDimBinValGetterFromOptsObjs([self.filterOptObj])

	def testNumberOfNebsClassificationOpts(self):
		expVals = [(0.1,),(0.2,)]
		actVals = self.binValGetter.getValsToBin(self.sparseMatrixCalculator)
		self.assertTrue( np.allclose(expVals,actVals) )


#TODO: Move thiss into the TestGetBinValsAtomFilteredVariousDistribs
class TestGetBinValsAtomFilteredRdfAndHozDists(unittest.TestCase):

	def setUp(self):
		#Geometry
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.coords = [  [0,0,1,"X"],
		                 [0,0,2,"A"],
		                 [0,0,3,"B"],
		                 [0,0,4,"C"],
		                 [0,0,5,"D"] ]

		#General options
		self.classBinResObj = binResHelp.getEmptyBinResultsFromMinMaxAndWidthStandard(-0.1,20,1,extremesAtCentre=False)
		self.fromIndices = [1,2,3,4]

		#Overall options object specific
		self.useGroups = [ [0,1] ]

		#Classification Options object specific
		self.distFilterIndices = [0]
		self.distFilterRanges = [ [-0.1,3.5], [3.5,5.5] ]

		#Distribution opts obj specific
		self.distrBinResObj = binResHelp.BinnedResultsStandard.fromBinEdges([-0.1,20])

		self.createTestObjs()

	def createTestObjs(self):
		#Geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams,lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Create the template distribution options
		currArgs = [self.distrBinResObj, self.fromIndices, self.fromIndices]
		self.distrOpts = [distrOptObjHelp.CalcRdfOptions(*currArgs) for x in self.useGroups]

		#Create the classification options
		currBins = [self.classBinResObj for x in self.distFilterRanges]
		currArgs = [currBins, self.fromIndices, self.distFilterIndices, self.distFilterRanges]
		self.classifierOpts = clsDistrOptObjs.AtomClassifyBasedOnDistsFromIndicesSimpleOpts(*currArgs)

		#
		currArgs = [self.fromIndices, self.classifierOpts, self.distrOpts, self.useGroups]
		self.filterOptObj = filteredComboOptObjHelp.FilteredAtomComboOptsObjGeneric(*currArgs)

		#create the sparse matrix calculator + populate it
		self.sparseMatrixCalculator = optsObjMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.filterOptObj])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

		#Create the binval getter
		self.binValGetter = optsObjMapHelp.getMultiDimBinValGetterFromOptsObjs([self.filterOptObj])

	def testExpectedA(self):
		expVals = [ (3,), (2,), (1,) ]
		actVals = self.binValGetter.getValsToBin(self.sparseMatrixCalculator)
		for expIter,actIter in it.zip_longest(expVals, actVals):
			[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expIter,actIter)]

	def testExpected_selfDistr(self):
		self.useGroups = [ [0,0] ]
		self.createTestObjs()

		distAA, distBB, distCC = 0,0,0
		distAB, distAC, distBC = 1,2,1
		distBA, distCA, distCB = distAB, distAC, distBC
		expVals = [ (distAA,), (distAB,), (distAC,), (distBA,), (distBB,), (distBC,),
		            (distCA,), (distCB,), (distCC,) ]
		actVals = self.binValGetter.getValsToBin(self.sparseMatrixCalculator)


		for expIter, actIter in it.zip_longest(expVals, actVals):
			[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expIter,actIter)]


	def testExpected_a_and_self_distr_byRef_classifiers(self):
		""" Test we get the expected distribution when using "byReference" classifiers """
		#Get first filter opts obj + set the classifier objects specifically
		filterObjA = self.filterOptObj #First set
		classifiersA = filteredAtomComboObjMaps.getClassifiersFromOptsObj(self.classifierOpts)
		filterObjA.classificationObjs = classifiersA

		#Get second filter opts obj; use byReference classifiers
		self.useGroups = [ [0,0] ]
		self.createTestObjs()
		filterObjB = self.filterOptObj
		filterObjB.classificationOpts = None #Force to use the objects
		classifiersB = classifierObjsHelp.getByReferenceClassifiers(classifiersA)
		filterObjB.classificationObjs = classifiersB

		#Run the functions - binValGetterA must always be run first
		binValGetterA = optsObjMapHelp.getMultiDimBinValGetterFromOptsObjs([filterObjA])
		binValGetterB = optsObjMapHelp.getMultiDimBinValGetterFromOptsObjs([filterObjB])

		actValsA = binValGetterA.getValsToBin(self.sparseMatrixCalculator)
		actValsB = binValGetterB.getValsToBin(self.sparseMatrixCalculator)

		#Compare actual and expected
		distAA, distBB, distCC = 0,0,0
		distAB, distAC, distBC = 1,2,1
		distBA, distCA, distCB = distAB, distAC, distBC

		expValsA = [ (3,), (2,), (1,) ]
		expValsB = [ (distAA,), (distAB,), (distAC,), (distBA,), (distBB,), (distBC,),
		            (distCA,), (distCB,), (distCC,) ]

		for expIter,actIter in it.zip_longest(expValsA, actValsA):
			[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expIter,actIter)]

		for expIter,actIter in it.zip_longest(expValsB, actValsB):
			[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expIter,actIter)]


	def testExpectedHozDistSelfDistr(self):
		self.useGroups = [ [0,0] ]
		self.coords = [  [0,1,1,"X"],
		                 [0,2.2,2,"A"],
		                 [0,3,3,"B"],
		                 [0,4,4,"C"],
		                 [0,6,5,"D"] ]
		self.distFilterRanges = [ [-0.1,5.7], [5.7,9] ] #First bin goes to just >sqrt(32)' meaning ABC encompassed
		self.createTestObjs()

		#Recreate the sparse matrix calculator/binval getter to use horizontal distances
		currArgs = [self.distrBinResObj, self.fromIndices, self.fromIndices]
		self.distrOpts = [distrOptObjHelp.CalcHozDistOptions(*currArgs, minDistVal=0.01) for x in self.useGroups]
		self.filterOptObj.distrOpts = self.distrOpts
		self.sparseMatrixCalculator = optsObjMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.filterOptObj])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)
		self.binValGetter = optsObjMapHelp.getMultiDimBinValGetterFromOptsObjs([self.filterOptObj])

		#Create the expected values
		distAA, distBB, distCC = 0, 0,0 
		distAB, distAC, distBC = 0.8, 1.8, 1
		distBA, distCA, distCB = distAB, distAC, distBC

		expVals = [ (distAA,), (distAB,), (distAC,), (distBA,), (distBB,), (distBC,),
		            (distCA,), (distCB,), (distCC,) ]
		actVals = self.binValGetter.getValsToBin(self.sparseMatrixCalculator)

		for expIter, actIter in it.zip_longest(expVals, actVals):
			[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expIter,actIter)]


	def testExpectedMinHozDistSelfDistr(self):
		self.useGroups = [ [0,0] ]
		self.coords = [  [0,1,1,"X"],
		                 [0,2.2,2,"A"],
		                 [0,3,3,"B"],
		                 [0,4,4,"C"],
		                 [0,6,5,"D"] ]
		self.distFilterRanges = [ [-0.1,5.7], [5.7,9] ] #First bin goes to just >sqrt(32)' meaning ABC encompassed
		self.createTestObjs()

		#Recreate the sparse matrix calculator/binval getter to use horizontal distances
		currArgs = [self.distrBinResObj, self.fromIndices, self.fromIndices]
		self.distrOpts = [distrOptObjHelp.CalcHozDistOptions(*currArgs, minDistVal=0.01, minDistAToB=True) for x in self.useGroups]
		self.filterOptObj.distrOpts = self.distrOpts
		self.sparseMatrixCalculator = optsObjMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.filterOptObj])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)
		self.binValGetter = optsObjMapHelp.getMultiDimBinValGetterFromOptsObjs([self.filterOptObj])

		#Create the expected values
		expVals = [ (0.8,), (0.8,), (1,) ]
		actVals = self.binValGetter.getValsToBin(self.sparseMatrixCalculator)

		for expIter, actIter in it.zip_longest(expVals, actVals):
			[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expIter,actIter)]


class TestGetBinValsWaterWaterHozDistsWithAdsHozDistFilter(unittest.TestCase):

	def setUp(self):
		#1) All geometric parameters for testing
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		#roll 90,pitch=45, OH len~1, HOH angle 104.5. Then just added translation vectors
		self.waterACoords = [ [0,0,0,"O"], [-0.13,0,0.99,"H"], [0.99,0,-0.13,"H"] ]
		#translation = [2,0,0]; no rotations
		self.waterBCoords = [ [2,0,0,"O"], [2.61, 0.79, 0, "H"], [2.61,-0.79,0,"H"] ]
		#
		self.waterCCoords = [ [3.0, 1.58, 0, "O"],  [3.83, 2.37, 0, "H"], [3.83, 0.79, 0, "H"] ]
		#
		self.waterDCoords = [ [3.0, -1.58, 0, "O"], [3.83, -0.79, 0, "H"], [3.83, -2.37, 0, "H"] ]

		self.adsCoords = [ [0.1,0,0.1,"X"], [2.1,1,0.1,"X"], [ 3,1,0.1,"X" ] ]
		self.coords = self.waterACoords + self.waterBCoords + self.waterCCoords + self.waterDCoords + self.adsCoords

		#General options
		self.oxyIndices = [0,3,6,9]
		self.hyIndices = [ [1,2], [4,5], [7,8], [10,11] ]
		self.minDist = True

		#Classification options
		self.distFilterIndices = [12,13,14]
		self.distFilterVals =   [ [0,1.5]    , [0,1.5] ] 
		self.adsHozDistRanges = [ [0,2.5], [4,10] ] 
		self.binResObjs = [None,None] #Irrelevent to these tests so....

		#filtered distr opts
		self.useGroups = [ [0,0] ]
		self.toIdxType = ["O"]

		self.createTestObjs()

	def createTestObjs(self):
		#Sort the geometry out
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams,lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Create a classification options object
		currArgs = [ self.binResObjs, self.oxyIndices, self.hyIndices, self.distFilterIndices, self.distFilterVals ]
		currKwargs = {"adsSiteMinHozToOtherAdsSiteRanges": self.adsHozDistRanges}
		self.classificationOpts = clsDistrOptObjs.WaterAdsorbedClassifier_usingMinHozDistsBetweenAdsorptionSitesOptsObj(*currArgs, **currKwargs)

		#Get the distribution option objects
		currArgs = [self.binResObjs, self.oxyIndices, self.oxyIndices]
		self.minHozDistOpts = [ distrOptObjHelp.CalcHozDistOptions(*currArgs, minDistAToB=self.minDist) for x in self.useGroups ]
		
		#Create the main options object
		currArgs = [self.oxyIndices, self.hyIndices, self.toIdxType, self.classificationOpts, self.minHozDistOpts, self.useGroups]
		self.optsObjFilteredCombo = filteredComboOptObjHelp.WaterToWaterFilteredAtomComboOptsObjGeneric(*currArgs)

		#create the sparse matrix calculator + populate it
		self.sparseMatrixCalculator = optsObjMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.optsObjFilteredCombo])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

		#Create the binval getter
		self.binValGetter = optsObjMapHelp.getMultiDimBinValGetterFromOptsObjs([self.optsObjFilteredCombo])


	def testExpectedA(self):
		minDistA, minDistB, minDistC = 2, 1.8698663053812163, 1.8698663053812163
		expVals = [(minDistA,), (minDistB,), (minDistC,)]
		actVals = self.binValGetter.getValsToBin(self.sparseMatrixCalculator)

		for expIter, actIter in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(exp,act,places=6) for exp,act in it.zip_longest(expIter,actIter)]

	def testExpectedA_runTwice_using_refClassiferObjs(self):
		#Make the first object
		filterOptsA = self.optsObjFilteredCombo 
		classifierObjsA = filteredAtomComboObjMaps.getClassifiersFromOptsObj(self.classificationOpts)
		filterOptsA.classificationOpts = None
		filterOptsA.classificationObjs = classifierObjsA

		#Second object using reference classifiers
		self.createTestObjs()
		filterOptsB = self.optsObjFilteredCombo
		classifierObjsB = classifierObjsHelp.getByReferenceClassifiers(classifierObjsA)
		filterOptsB.classificationOpts = None
		filterOptsB.classificationObjs = classifierObjsB

		#Load expected values for each
		minDistA, minDistB, minDistC = 2, 1.8698663053812163, 1.8698663053812163
		expVals = [(minDistA,), (minDistB,), (minDistC,)]

		#Run each binval getter in turn
		binValGetterA = optsObjMapHelp.getMultiDimBinValGetterFromOptsObjs([filterOptsA])
		binValGetterB = optsObjMapHelp.getMultiDimBinValGetterFromOptsObjs([filterOptsB])

		actValsA = binValGetterA.getValsToBin(self.sparseMatrixCalculator)
		actValsB = binValGetterB.getValsToBin(self.sparseMatrixCalculator)

		for expIter,actIterA, actIterB in it.zip_longest(expVals, actValsA, actValsB):
			[self.assertAlmostEqual(exp,act, places=6) for exp,act in it.zip_longest(expIter,actIterA)]
			[self.assertAlmostEqual(exp,act, places=6) for exp,act in it.zip_longest(expIter,actIterB)]


	def testExpectedHozRdfBinning(self):
		#Create the opts objs 
		self.minDist = False
		self.createTestObjs()

		#Figure out expected
		distAA, distBB, distCC = 0, 0, 0
		distAB, distAC, distBC = 2, 3.3906341589738047, 1.8698663053812163
		distBA, distCA, distCB = distAB, distAC, distBC

		expVals = [ (distAA,), (distAB,), (distAC,), (distBA,), (distBB,), (distBC,),
		            (distCA,), (distCB,), (distCC,) ]
		actVals = self.binValGetter.getValsToBin(self.sparseMatrixCalculator)

		#Run + compare
		for expIter, actIter in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(exp,act, places=6) for exp,act in it.zip_longest(expIter,actIter)]


class TestGetBinValsWaterMinDists(unittest.TestCase):


	def setUp(self):
		#Define the geometry
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		#roll 90,pitch=45, OH len~1, HOH angle 104.5. Then just added translation vectors
		self.waterACoords = [ [1,0,0,"O"], [-1.13,0,0.99,"H"], [1.99,0,-0.13,"H"] ]
		#translation = [2,0,0]; no rotations
		self.waterBCoords = [ [2,0,0,"O"], [2.61, 0.79, 0, "H"], [2.61,-0.79,0,"H"] ]
		#translation = [2+(2*0.61), 2*0.79, 0]
		self.waterCCoords = [ [3.22, 1.58, 0, "O"],  [3.83, 2.37, 0, "H"], [3.83, 0.79, 0, "H"] ]
		#translation = [2+(2*0.61), -2*0.79,0]
		self.waterDCoords = [ [3.22, -1.58, 0, "O"], [3.83, -0.79, 0, "H"], [3.83, -2.37, 0, "H"] ]

		self.xCoord = [[0,0,0,"X"]]
		self.coords = self.xCoord + self.waterACoords + self.waterBCoords + self.waterCCoords + self.waterDCoords 

		#General options
		self.oxyIndices = [1,4,7,10] 
		self.hyIndices = [ [2,3], [5,6], [8,9], [11,12] ]

		#Filter/classification options
		self.distFilterIndices = [0]
		self.filterObjA = clsDistrOptObjs.WaterMinDistAndHBondsFilterObj(distFilterRange=[0,2.1])
		self.filterObjB = clsDistrOptObjs.WaterMinDistAndHBondsFilterObj(distFilterRange=[2.1,10])
		self.filterObjC = clsDistrOptObjs.WaterMinDistAndHBondsFilterObj(distFilterRange=[2.1,10])
		self.filterObjs  = [self.filterObjA, self.filterObjB, self.filterObjC]

		#Min dist options (Note we only need a binResObj for this one)
		self.binResObjA = binResHelp.BinnedResultsStandard.fromBinEdges([0,5,10]) #I COULD probably use a dud but...
		self.minDistType = "O"
		self.minDistVal = 0.1

		#Options for the filtered atom combo opts object
		self.toIdxType = ["O"]
		self.useGroups = [ [0,1] ]

		self.createTestObjs()


	def createTestObjs(self):
		#Create the geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Create the classifier options obj
		currArgs = [ [self.binResObjA], self.oxyIndices, self.hyIndices, self.distFilterIndices, self.filterObjs] #The binResObjs should be irrelevant really
		self.classifierOptsObj = clsDistrOptObjs.WaterCountTypesMinDistAndHBondSimpleOpts.fromFilterObjs(*currArgs)

		#Create the min-dist distribution opts object [Note that we NEED all indices as toIndices to handle certain types of populators]

		currArgs = [self.binResObjA, self.oxyIndices, self.hyIndices, self.oxyIndices]
		currKwargs = {"primaryIdxType":"O", "minDistType":self.minDistType, "minVal":self.minDistVal}
#		self.minDistOptsObj = distrOptObjHelp.WaterMinDistOptions( *currArgs, **currKwargs )

		self.minDistOptsObjs = [ distrOptObjHelp.WaterMinDistOptions( *currArgs, **currKwargs ) for x in self.useGroups ]


		#Create the main options object
		currArgs = [self.oxyIndices, self.hyIndices, self.toIdxType, self.classifierOptsObj, self.minDistOptsObjs, self.useGroups] 
		self.optsObjFilteredCombo = filteredComboOptObjHelp.WaterToWaterFilteredAtomComboOptsObjGeneric(*currArgs)

		#create the sparse matrix calculator + populate it
		self.sparseMatrixCalculator = optsObjMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.optsObjFilteredCombo])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

		#Create the binval getter
		self.binValGetter = optsObjMapHelp.getMultiDimBinValGetterFromOptsObjs([self.optsObjFilteredCombo])

		#Create a second 2-d binval getter

	def testExpectedValsA(self):
		#NOTE: Each value is an iter since we're using the multi-dim binval getter. This iter structure makes more sense when using >1-d bins as examples
		expVals = [ (2.7248486196484385,),(1.9961963831246665,) ]
		actVals = self.binValGetter.getValsToBin(self.sparseMatrixCalculator)

		for expIter, actIter in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(exp,act,places=6) for exp,act in it.zip_longest(expIter,actIter)]

	def testExpected_selfDistrib_withMinDist(self):
		""" Test we get expected values for distribution between the same groups; want to avoid getting zero for all """
		self.minDistVal = 0.1
		self.useGroups = [ [0,0] ]
		self.createTestObjs()
		expVals = [ (1,), (1,) ]
		actVals = self.binValGetter.getValsToBin(self.sparseMatrixCalculator)

		for expIter, actIter in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(exp,act,places=6) for exp,act in it.zip_longest(expIter, actIter)]


	def testExpected_selfDistrib_withNegativeMinDist(self):
		self.minDistVal = -0.1
		self.useGroups = [ [0,0] ]
		self.createTestObjs()
		expVals = [ (0,), (0,) ]
		actVals = self.binValGetter.getValsToBin(self.sparseMatrixCalculator)

		for expIter, actIter in it.zip_longest(expVals, actVals):
			[self.assertAlmostEqual(exp,act,places=6) for exp,act in it.zip_longest(expIter, actIter)]


	def testExpected_toAll(self):
		self.useGroups = [ [1,0] ] #Means O-H is the closest contact between groups
		self.toIdxType = ["all"]
		self.createTestObjs()

		#Note that this also happens to be the O-H bondlength in all; which is a bit annoying but....
		expVals = [ (0.9980981915623334,), (0.9980981915623334,) ]
		actVals = self.binValGetter.getValsToBin(self.sparseMatrixCalculator)

		for expIter, actIter in it.zip_longest(expVals, actVals):
			[self.assertAlmostEqual(exp,act,places=6) for exp,act in it.zip_longest(expIter,actIter)]

	def testExpected_toHydrogen(self):
		self.useGroups = [ [0,1] ]
		self.toIdxType = ["h"]
		self.createTestObjs()
		expVals = [(2.938196725884773,),(1.9932385707686877,)]
		actVals = self.binValGetter.getValsToBin(self.sparseMatrixCalculator)

		for expIter, actIter in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expIter,actIter)]

	def testExpected_twoDistrBins(self):
		self.useGroups = [ [0,1], [0,1] ]
		self.toIdxType = [ "O", "H" ] 
		self.createTestObjs()

		expVals = [ (2.7248486196484385,2.938196725884773),
		            (1.9961963831246665,1.9932385707686877) ] #Each element is the 2-d value to bin for EACH ATOM
		actVals = self.binValGetter.getValsToBin(self.sparseMatrixCalculator)

		for expIter, actIter in it.zip_longest(expVals, actVals):
			[self.assertAlmostEqual(exp,act,places=6) for exp,act in it.zip_longest(expIter,actIter)]

	def testRaisesForInconsistentUseGroups(self):
		self.useGroups = [ [0,1], [1,0] ] #Makes no sense for binning purposes
		self.toIdxType = [ "O", "H" ]
		with self.assertRaises(ValueError):
			self.createTestObjs()



class TestGetBinvalsWaterWaterFilteredRdf(unittest.TestCase):

	#Basically stolen from above
	def setUp(self):
		#Define the geometry
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		#roll 90,pitch=45, OH len~1, HOH angle 104.5. Then just added translation vectors
		self.waterACoords = [ [1,0,0,"O"], [-1.13,0,0.99,"H"], [1.99,0,-0.13,"H"] ]
		#translation = [2,0,0]; no rotations
		self.waterBCoords = [ [2,0,0,"O"], [2.61, 0.79, 0, "H"], [2.61,-0.79,0,"H"] ]
		#translation = [2+(2*0.61), 2*0.79, 0]
		self.waterCCoords = [ [3.22, 1.58, 0, "O"],  [3.83, 2.37, 0, "H"], [3.83, 0.79, 0, "H"] ]
		#translation = [2+(2*0.61), -2*0.79,0]
		self.waterDCoords = [ [3.22, -1.58, 0, "O"], [3.83, -0.79, 0, "H"], [3.83, -2.37, 0, "H"] ]

		self.xCoord = [[0,0,0,"X"]]
		self.coords = self.xCoord + self.waterACoords + self.waterBCoords + self.waterCCoords + self.waterDCoords 

		#General options
		self.oxyIndices = [1,4,7,10] 
		self.hyIndices = [ [2,3], [5,6], [8,9], [11,12] ]

		#Filter/classification options
		self.distFilterIndices = [0]
		self.filterObjA = clsDistrOptObjs.WaterMinDistAndHBondsFilterObj(distFilterRange=[0,2.1])
		self.filterObjB = clsDistrOptObjs.WaterMinDistAndHBondsFilterObj(distFilterRange=[2.1,10])
		self.filterObjC = clsDistrOptObjs.WaterMinDistAndHBondsFilterObj(distFilterRange=[2.1,10])
		self.filterObjs  = [self.filterObjA, self.filterObjB, self.filterObjC]

		#Min dist options (Note we only need a binResObj for this one)
		self.binResObjA = binResHelp.BinnedResultsStandard.fromBinEdges([0.1,5,10]) #Needs to be sensible for default, since we filter vals to bin by whether dists are in these bins
		self.minDistType = "O"

		#Options for the filtered atom combo opts object
		self.toIdxType = ["O"]
		self.useGroups = [ [0,1] ]

		self.createTestObjs()

	def createTestObjs(self):
		#Create the geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Create the classifier options obj
		currArgs = [ [self.binResObjA], self.oxyIndices, self.hyIndices, self.distFilterIndices, self.filterObjs] #The binResObjs should be irrelevant really
		self.classifierOptsObj = clsDistrOptObjs.WaterCountTypesMinDistAndHBondSimpleOpts.fromFilterObjs(*currArgs)

		#Create the distribution opts object
		#Note that the actual indices we pass here are really just placeholders; the actual used indices are determined by the opts-object values
		#(Though oxyIndices may need to be correct.....)
		currArgs = [self.binResObjA, self.oxyIndices, self.oxyIndices]
		currKwargs = {"minDistAToB":False}
		self.rdfDistOptsObjs = [ distrOptObjHelp.CalcRdfOptions( *currArgs, **currKwargs ) for x in self.useGroups ] #Multiple lets us get rdfs between various groups

		#Create the main options object
		currArgs = [self.oxyIndices, self.hyIndices, self.toIdxType, self.classifierOptsObj, self.rdfDistOptsObjs, self.useGroups] 
		self.optsObjFilteredCombo = filteredComboOptObjHelp.WaterToWaterFilteredAtomComboOptsObjGeneric(*currArgs)

		#create the sparse matrix calculator + populate it
		self.sparseMatrixCalculator = optsObjMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.optsObjFilteredCombo])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

		#Create the binval getter
		self.binValGetter = optsObjMapHelp.getMultiDimBinValGetterFromOptsObjs([self.optsObjFilteredCombo])

	def _runTestFunct(self):
		return self.binValGetter.getValsToBin(self.sparseMatrixCalculator)

	def testExpected_oneDimToOxy(self):
		#Dists are only for oxygen-oxygen
		distAC, distAD = 2.7248486196484385, 2.7248486196484385
		distBC, distBD = 1.9961963831246665, 1.9961963831246665

		expVals =  [ (distAC,), (distAD,), (distBC,), (distBD,) ]
		actVals = self._runTestFunct()

		for expIter,actIter in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(exp,act, places=6) for exp,act in it.zip_longest(expIter,actIter)]


	def testExpected_oneDimToHy(self):
		self.toIdxType = [ "H" ] 
		self.createTestObjs()

		#Distances from oxygen to hydrogen
		distAC1, distAC2 = 3.6913141291415448 , 2.938196725884773 
		distAD1, distAD2 = 2.938196725884773, 3.6913141291415448 
		distBC1, distBC2 = 2.9942945746869998, 1.9932385707686877
		distBD1, distBD2 = 1.9932385707686877, 2.9942945746869998 

		expVals = [ (distAC1,), (distAC2,), (distAD1,), (distAD2,),
		            (distBC1,), (distBC2,), (distBD1,), (distBD2,) ]
		actVals = self._runTestFunct()

		for expIter, actIter in it.zip_longest(expVals, actVals):
			[self.assertAlmostEqual(exp,act, places=6) for exp,act in it.zip_longest(expIter,actIter)]

	def testExpected_twoDim(self):
		self.toIdxType = ["O","H"]
		self.useGroups = [ [0,1], [0,1] ]
		self.createTestObjs()

		#Distances from oxygen to oxygen
		distAC_Oxy, distAD_Oxy = 2.7248486196484385, 2.7248486196484385
		distBC_Oxy, distBD_Oxy = 1.9961963831246665, 1.9961963831246665


		#Distances from oxygen to hydrogen
		distAC1, distAC2 = 3.6913141291415448 , 2.938196725884773 
		distAD1, distAD2 = 2.938196725884773, 3.6913141291415448 
		distBC1, distBC2 = 2.9942945746869998, 1.9932385707686877
		distBD1, distBD2 = 1.9932385707686877, 2.9942945746869998 

		expVals = [ (distAC_Oxy, distAC1), (distAC_Oxy, distAC2), (distAC_Oxy, distAD1), (distAC_Oxy, distAD2),
		            (distAD_Oxy, distAC1), (distAD_Oxy, distAC2), (distAD_Oxy, distAD1), (distAD_Oxy, distAD2),
		            (distBC_Oxy, distBC1), (distBC_Oxy, distBC2), (distBC_Oxy, distBD1), (distBC_Oxy, distBD2),
		            (distBD_Oxy, distBC1), (distBD_Oxy, distBC2), (distBD_Oxy, distBD1), (distBD_Oxy, distBD2) ]

		actVals = self._runTestFunct()

		for expIter, actIter in it.zip_longest(expVals, actVals):
			[self.assertAlmostEqual(exp,act, places=6) for exp,act in it.zip_longest(expIter,actIter)]

class TestGetBinValsWaterWaterHozRdf(unittest.TestCase):

	#Basically stolen from above
	def setUp(self):
		#Define the geometry
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		#roll 90,pitch=45, OH len~1, HOH angle 104.5. Then just added translation vectors
		self.waterACoords = [ [1,0,1,"O"], [-1.13,0,0.99,"H"], [1.99,0,-0.13,"H"] ]
		#translation = [2,0,0]; no rotations
		self.waterBCoords = [ [2,0,1,"O"], [2.61, 0.79, 0, "H"], [2.61,-0.79,0,"H"] ]
		#translation = [2+(2*0.61), 2*0.79, 0]
		self.waterCCoords = [ [3.22, 1.58, 0, "O"],  [3.83, 2.37, 0, "H"], [3.83, 0.79, 0, "H"] ]
		#translation = [2+(2*0.61), -2*0.79,0]
		self.waterDCoords = [ [3.22, -1.58, 0, "O"], [3.83, -0.79, 0, "H"], [3.83, -2.37, 0, "H"] ]

		self.xCoord = [[0,0,0,"X"]]
		self.coords = self.xCoord + self.waterACoords + self.waterBCoords + self.waterCCoords + self.waterDCoords 

		#General options
		self.oxyIndices = [1,4,7,10] 
		self.hyIndices = [ [2,3], [5,6], [8,9], [11,12] ]

		#Filter/classification options
		self.distFilterIndices = [0]
		self.filterObjA = clsDistrOptObjs.WaterMinDistAndHBondsFilterObj(distFilterRange=[0,3.1])
		self.filterObjB = clsDistrOptObjs.WaterMinDistAndHBondsFilterObj(distFilterRange=[3.1,10])
		self.filterObjC = clsDistrOptObjs.WaterMinDistAndHBondsFilterObj(distFilterRange=[3.1,10])
		self.filterObjs  = [self.filterObjA, self.filterObjB, self.filterObjC]

		#Min dist options (Note we only need a binResObj for this one)
		self.binResObjA = binResHelp.BinnedResultsStandard.fromBinEdges([0.1,5,10]) #Needs to be sensible for default, since we filter vals to bin by whether dists are in these bins
		self.minDistCase = False
		#Options for the filtered atom combo opts object
		self.toIdxType = ["O"]
		self.useGroups = [ [0,1] ]

		self.createTestObjs()

	def createTestObjs(self):
		#Create the geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Create the classifier options obj
		currArgs = [ [self.binResObjA], self.oxyIndices, self.hyIndices, self.distFilterIndices, self.filterObjs] #The binResObjs should be irrelevant really
		self.classifierOptsObj = clsDistrOptObjs.WaterCountTypesMinDistAndHBondSimpleOpts.fromFilterObjs(*currArgs)

		#Create the distribution opts object; indices are just placeholders
		currArgs = [self.binResObjA, self.oxyIndices, self.oxyIndices]
		self.hozDistRdfOptsObjs = [distrOptObjHelp.CalcHozDistOptions(*currArgs, minDistAToB=self.minDistCase) for x in self.useGroups]

		#Create the main options object
		currArgs = [self.oxyIndices, self.hyIndices, self.toIdxType, self.classifierOptsObj, self.hozDistRdfOptsObjs, self.useGroups] 
		self.optsObjFilteredCombo = filteredComboOptObjHelp.WaterToWaterFilteredAtomComboOptsObjGeneric(*currArgs)

		#create the sparse matrix calculator + populate it
		self.sparseMatrixCalculator = optsObjMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.optsObjFilteredCombo])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

		#Create the binval getter
		self.binValGetter = optsObjMapHelp.getMultiDimBinValGetterFromOptsObjs([self.optsObjFilteredCombo])

	def _runTestFunct(self):
		return self.binValGetter.getValsToBin(self.sparseMatrixCalculator)


	def testExpectedA(self):
		#For oxy only
		distAC, distAD = 2.7248486196484385, 2.7248486196484385
		distBC, distBD = 1.9961963831246665, 1.9961963831246665

		expVals =  [ (distAC,), (distAD,), (distBC,), (distBD,) ]
		actVals = self._runTestFunct()

		for expIter,actIter in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(exp,act, places=6) for exp,act in it.zip_longest(expIter,actIter)]

	def testMinDist(self):
		self.waterDCoords[0][0] -= 0.1 #So the two arent equidistant
		self.minDistCase = True
		self.createTestObjs()

		minDistA = 2.6440121028467325
		minDistB = 1.9366982212001953

		expVals = [ (minDistA,), (minDistB,) ]
		actVals = self._runTestFunct()

		for expIter,actIter in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(exp,act, places=6) for exp,act in it.zip_longest(expIter,actIter)]

	def testMinDists_overlappingGroups(self):
		self.useGroups = [ [0,0] ]
		self.minDistCase = True
		self.createTestObjs()

		minDistAB = 1
		expVals = [ (minDistAB,), (minDistAB,) ]
		actVals = self._runTestFunct()

		for expIter,actIter in it.zip_longest(expVals,actVals):
			[self.assertAlmostEqual(exp,act, places=6) for exp,act in it.zip_longest(expIter,actIter)]

	def testTwoDimVersionA(self):
		self.toIdxType = ["O","H"]
		self.useGroups = [ [0,1], [0,1] ]
		self.createTestObjs()

		#Distances from oxygen to oxygen
		distAC_Oxy, distAD_Oxy = 2.7248486196484385, 2.7248486196484385
		distBC_Oxy, distBD_Oxy = 1.9961963831246665, 1.9961963831246665


		#Distances from oxygen to hydrogen
		distAC1, distAC2 = 3.6913141291415448 , 2.938196725884773 
		distAD1, distAD2 = 2.938196725884773, 3.6913141291415448 
		distBC1, distBC2 = 2.9942945746869998, 1.9932385707686877
		distBD1, distBD2 = 1.9932385707686877, 2.9942945746869998 

		expVals = [ (distAC_Oxy, distAC1), (distAC_Oxy, distAC2), (distAC_Oxy, distAD1), (distAC_Oxy, distAD2),
		            (distAD_Oxy, distAC1), (distAD_Oxy, distAC2), (distAD_Oxy, distAD1), (distAD_Oxy, distAD2),
		            (distBC_Oxy, distBC1), (distBC_Oxy, distBC2), (distBC_Oxy, distBD1), (distBC_Oxy, distBD2),
		            (distBD_Oxy, distBC1), (distBD_Oxy, distBC2), (distBD_Oxy, distBD1), (distBD_Oxy, distBD2) ]

		actVals = self._runTestFunct()

		for expIter, actIter in it.zip_longest(expVals, actVals):
			[self.assertAlmostEqual(exp,act, places=6) for exp,act in it.zip_longest(expIter,actIter)]





class TestGetBinvalsWaterWaterFilteredVariousDistribs(unittest.TestCase):

	#Basically stolen from above
	def setUp(self):
		#Define the geometry
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		#roll 90,pitch=45, OH len~1, HOH angle 104.5. Then just added translation vectors
		self.waterACoords = [ [1,0,0,"O"], [-1.13,0,0.99,"H"], [1.99,0,-0.13,"H"] ]
		#translation = [2,0,0]; no rotations
		self.waterBCoords = [ [2,0,0,"O"], [2.61, 0.79, 0, "H"], [2.61,-0.79,0,"H"] ]
		#translation = [2+(2*0.61), 2*0.79, 0]
		self.waterCCoords = [ [3.22, 1.58, 0, "O"],  [3.83, 2.37, 0, "H"], [3.83, 0.79, 0, "H"] ]
		#translation = [2+(2*0.61), -2*0.79,0]
		self.waterDCoords = [ [3.22, -1.58, 0, "O"], [3.83, -0.79, 0, "H"], [3.83, -2.37, 0, "H"] ]

		self.xCoord = [[0,0,0,"X"]]
		self.coords = self.xCoord + self.waterACoords + self.waterBCoords + self.waterCCoords + self.waterDCoords 

		#General options
		self.oxyIndices = [1,4,7,10] 
		self.hyIndices = [ [2,3], [5,6], [8,9], [11,12] ]

		#Filter/classification options
		self.distFilterIndices = [0]
		self.filterObjA = clsDistrOptObjs.WaterMinDistAndHBondsFilterObj(distFilterRange=[0,2.1])
		self.filterObjB = clsDistrOptObjs.WaterMinDistAndHBondsFilterObj(distFilterRange=[2.1,10])
		self.filterObjC = clsDistrOptObjs.WaterMinDistAndHBondsFilterObj(distFilterRange=[2.1,10])
		self.filterObjs  = [self.filterObjA, self.filterObjB, self.filterObjC]

		#Min dist options (Note we only need a binResObj for this one)
		self.binResObjA = binResHelp.BinnedResultsStandard.fromBinEdges([0.1,5,10]) #Needs to be sensible for default, since we filter vals to bin by whether dists are in these bins
		self.minDistType = "O"

		#Options for the filtered atom combo opts object
		self.toIdxType = ["O"]
		self.useGroups = [ [0,1] ]

		#Create the distribution opts object
		#Note that the actual indices we pass here are really just placeholders; the actual used indices are determined by the opts-object values
		#(Though oxyIndices may need to be correct.....)
		currArgs = [self.binResObjA, self.oxyIndices, self.oxyIndices]
		currKwargs = {"minDistAToB":False}
		self.distrOptObjs = [ distrOptObjHelp.CalcRdfOptions( *currArgs, **currKwargs ) for x in self.useGroups ] #Multiple lets us get rdfs between various groups

		self.createTestObjs()

	def createTestObjs(self):
		#Create the geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Create the classifier options obj
		currArgs = [ [self.binResObjA], self.oxyIndices, self.hyIndices, self.distFilterIndices, self.filterObjs] #The binResObjs should be irrelevant really
		self.classifierOptsObj = clsDistrOptObjs.WaterCountTypesMinDistAndHBondSimpleOpts.fromFilterObjs(*currArgs)


		#Create the main options object
		currArgs = [self.oxyIndices, self.hyIndices, self.toIdxType, self.classifierOptsObj, self.distrOptObjs, self.useGroups] 
		self.optsObjFilteredCombo = filteredComboOptObjHelp.WaterToWaterFilteredAtomComboOptsObjGeneric(*currArgs)

		#create the sparse matrix calculator + populate it
		self.sparseMatrixCalculator = optsObjMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.optsObjFilteredCombo])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

		#Create the binval getter
		self.binValGetter = optsObjMapHelp.getMultiDimBinValGetterFromOptsObjs([self.optsObjFilteredCombo])

	def _runTestFunct(self):
		return self.binValGetter.getValsToBin(self.sparseMatrixCalculator)

	def testExpectedOrientations(self):
		#Alter water orientations a bit
		self.waterACoords = [ [1,0,0,"O"], [1.4329029902719427, 0.7906895737438433, 0.4329029902719426,"H"],
		                      [1.4329029902719427, -0.7906895737438433, 0.4329029902719426, 'H'] ] #Pitch is 45
		self.waterBCoords = [ [2,0,0,"O"], [2.01, 0.79, 0.61, "H"], [2.01,-0.79,0.61,"H"] ] #Pitch is 89 ish
		self.coords = self.xCoord + self.waterACoords + self.waterBCoords + self.waterCCoords + self.waterDCoords 

		#Create the distribution object
		self.useGroups = [ [0] ]
		currArgs = [self.binResObjA, self.oxyIndices, self.hyIndices]
		self.distrOptObjs = [ distrOptObjHelp.WaterOrientationOptions(*currArgs, angleType="pitch") for x in self.useGroups ]
		self.createTestObjs()

		#Figure out expected values
		expVals = [(45,), (89.0608090542646,)]
		actVals = self._runTestFunct()

		#Run and check expected/actual are equal
		for expIter, actIter in it.zip_longest(expVals, actVals):
			[self.assertAlmostEqual(exp,act, places=6) for exp,act in it.zip_longest(expIter,actIter)]


	def testExpected_totalNumberHBonds(self):
		#Create the distribution object
		self.useGroups = [ [0,1] ]
		currArgs = [self.binResObjA, self.oxyIndices, self.hyIndices, self.oxyIndices, self.hyIndices]
		currKwargs = {"acceptor":True, "donor":True}
		self.distrOptObjs = [distrOptObjHelp.CountHBondsBetweenWaterGroupsOptions(*currArgs, **currKwargs) for x in self.useGroups]
		self.createTestObjs()

		#Figure out what values we expect
		expVals = [ (0,), (2,) ]
		actVals = self._runTestFunct()

		#Run + check expected/actual are equal
		for expIter, actIter in it.zip_longest(expVals, actVals):
			[self.assertAlmostEqual(exp,act, places=6) for exp,act in it.zip_longest(expIter,actIter)]


	def testExpectedPlanarDists_toOxyI(self):
		#Mod coords; A/B should stay within 2.1 (just the oxy)
		self.waterACoords = [ [0,0,1.1,"O"], [-1.13,0,0.99,"H"], [1.99,0,-0.13,"H"] ]
		#translation = [2,0,0]; no rotations
		self.waterBCoords = [ [1,0,1.5,"O"], [2.61, 0.79, 0, "H"], [2.61,-0.79,0,"H"] ]
		self.coords = self.xCoord + self.waterACoords + self.waterBCoords + self.waterCCoords + self.waterDCoords 

		#Mod distr
		self.useGroups = [ [0] ]
		currArgs = [ self.binResObjA, self.oxyIndices ]
		self.distrOptObjs = [distrOptObjHelp.CalcPlanarDistOptions(*currArgs)]
		self.createTestObjs()

		#Figure out what values we expect
		expVals = [ (1.1,), (1.5,) ]
		actVals = self._runTestFunct()

		#Run + check expected/actual are equal
		for expIter, actIter in it.zip_longest(expVals, actVals):
			[self.assertAlmostEqual(exp,act, places=6) for exp,act in it.zip_longest(expIter, actIter)]


	def testPlanarDists_hyIdxType(self):
		#Mod coords; A/B should stay within 2.1 (just the oxy)
		self.waterACoords = [ [0,0,1.1,"O"], [-1.13,0,0.99,"H"], [1.99,0,-0.13,"H"] ]
		#translation = [2,0,0]; no rotations
		self.waterBCoords = [ [1,0,1.5,"O"], [2.61, 0.79, 0, "H"], [2.61,-0.79,0,"H"] ]
		self.coords = self.xCoord + self.waterACoords + self.waterBCoords + self.waterCCoords + self.waterDCoords 

		#Mod distr
		self.useGroups = [ [0] ]
		self.toIdxType = ["H"]
		currArgs = [ self.binResObjA, self.oxyIndices ]
		self.distrOptObjs = [distrOptObjHelp.CalcPlanarDistOptions(*currArgs)]
		self.createTestObjs()

		#Figure out what values we expect
		expVals = [ (0.99,), (0.13,), (0,), (0,) ]
		actVals = self._runTestFunct()

		#Run + check expected/actual are equal
		for expIter, actIter in it.zip_longest(expVals, actVals):
			[self.assertAlmostEqual(exp,act, places=6) for exp,act in it.zip_longest(expIter, actIter)]


class TestGenericNonHyAndHyFilteredAtomComboDistribs_varyDistrOpts(unittest.TestCase):

	def setUp(self):
		#Define the geometry; 2 "free" water/hydroxyl and 2 hydrogen bonded to each other
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		self.hydroxylA = [ [0,0,0,"O"], [0,1,0.2,"H"] ]
		self.hydroxylB = [ [5,0,0.1,"O"], [5,1,0.3,"H"] ]
		self.waterA = [ [1,2,0,"O"], [0.5,3,0,"H"], [1.5,3,0,"H"] ]
		self.waterB = [ [6,6,0,"O"], [5.5,7,0,"H"], [6.5,7,0,"H"] ]

		self.cartCoords = self.hydroxylA + self.hydroxylB + self.waterA + self.waterB

		#Options for the classifier
		self.binResObjA = binResHelp.BinnedResultsStandard.fromBinEdges([-0.5,0.5,1.5,2.5,3.5,4.5])
		self.fromNonHy = [ [0], [2] ]
		self.fromHy    = [ [1], [3] ]
		self.toNonHy = [ [4], [7] ]
		self.toHy    = [ [5,6],[8,9] ]

		#Only 1st hydroxyl in 1st group; 2nd hydroxyl in the other
		self.nTotalFilterRanges = [ [0.5,100], [-0.5,0.5] ] #>0 h-bonds and 0 hbonds are the criteria
		self.maxOOHBond = 3
		self.maxAngleHBond = 50

		#Create the distribution object; which will change in each test but...
		currArgs = [ self.binResObjA, [x[0] for x in self.toNonHy] ] #Neither arg should actually matter; hence I've set the wrong value for the 2nd arg on purpose
		self.distrOptObjs = [distrOptObjHelp.CalcPlanarDistOptions(*currArgs, planeEqn=planeEqnHelp.ThreeDimPlaneEquation(0,0,1,4))]
		self.useGroups = [ [0] ]
		self.useNonHyIdx, self.useIdxEach = True, 0

		self.createTestObjs()

	def createTestObjs(self):
		#Create the geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

		#Create the classifier options object
		currArgs = [ [self.binResObjA, self.binResObjA], self.fromNonHy, self.fromHy, self.toNonHy, self.toHy ]
		currKwargs = {"nTotalFilterRanges":self.nTotalFilterRanges, "maxOOHBond":self.maxOOHBond, "maxAngleHBond":self.maxAngleHBond}
		self.classifierOptsObj = clsDistrOptObjs.ClassifyBasedOnHBondingToGroup_simple(*currArgs, **currKwargs)

		#Create the main options object
		currArgs = [self.fromNonHy, self.fromHy, self.classifierOptsObj, self.distrOptObjs, self.useGroups]
		currKwargs = {"useNonHyIdx":self.useNonHyIdx, "useIdxEach":self.useIdxEach}
		self.optsObjFilteredCombo = filteredComboOptObjHelp.GenericNonHyAndHyFilteredOptsObj_simple(*currArgs, **currKwargs)

		#Create the sparse matrix calculator and populate it
		self.sparseMatrixCalculator = optsObjMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.optsObjFilteredCombo])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

		#Create the binval getter
		self.binValGetter = optsObjMapHelp.getMultiDimBinValGetterFromOptsObjs([self.optsObjFilteredCombo])

	def _runTestFunct(self):
		return self.binValGetter.getValsToBin(self.sparseMatrixCalculator)

	def testExpectedPlanarDistribs_toHydroxylOxygen(self):
		expVals = [ (4,) ]
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)

	def testExpectedPlanarDistribs_toBothHydroxylHydrogen(self):
		#Change options
		self.useNonHyIdx = False
		self.nTotalFilterRanges = [ [-0.5,100], [-0.5,0.5] ]
		self.createTestObjs()

		expVals = [ (3.8,), (3.7,) ]
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)

	def testExpectedWaterOrientation_singleWater(self):
		#Create the distribution
		self.fromNonHy, self.toNonHy = self.toNonHy, self.fromNonHy
		self.fromHy, self.toHy = self.toHy, self.fromHy

		currBinResObj = binResHelp.getEmptyBinResultsFromMinMaxAndWidthStandard(-90,90,10,extremesAtCentre=False)
		fakeOxyIndices, fakeHyIndices = [ [1], [2,3] ] #Settings wrong indices here SHOULDNT matter
		currArgs = [currBinResObj, fakeOxyIndices, fakeHyIndices]
		currKwargs = {"angleType":"pitch"}
		self.distrOptObjs = [ distrOptObjHelp.WaterOrientationOptions(*currArgs,**currKwargs) ]
		self.createTestObjs()

		#
		expVals = [ (0,) ]
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)

	def testWaterOrientation_raisesForHydroxylGroups(self):
		""" Want to raise an Assertion error if we try to get water orientation data for hydroxyl groups """

		currBinResObj = binResHelp.getEmptyBinResultsFromMinMaxAndWidthStandard(-90,90,10,extremesAtCentre=False)
		fakeOxyIndices, fakeHyIndices = [ [1], [2,3] ] #Settings wrong indices here SHOULDNT matter
		currArgs = [currBinResObj, fakeOxyIndices, fakeHyIndices]
		currKwargs = {"angleType":"pitch"}
		self.distrOptObjs = [ distrOptObjHelp.WaterOrientationOptions(*currArgs,**currKwargs) ]

		with self.assertRaises(AssertionError):
			self.createTestObjs()

	def testMinHozDistFixedGroup(self):

		#Create the distribution object; only the toIndices should matter
		currBinResObj = binResHelp.getEmptyBinResultsFromMinMaxAndWidthStandard(0,10,1) #Shouldnt REALLY matter much
		dudIndicesFrom, indicesTo = [0], [4,5,6,7,8,9]
		self.useGroups = [ [1] ] #This has a slight non-zero z-value so...

		currArgs = [currBinResObj, dudIndicesFrom, indicesTo]
		self.distrOptObjs = [distrOptObjHelp.CalcHozDistOptions(*currArgs,minDistAToB=True)]
		self.createTestObjs()

		#Run + check
		expVals = [ (math.sqrt( (3**2) + (0.5**2) ),) ]
		actVals = self._runTestFunct()

		self.assertTrue( np.allclose( np.array(expVals), np.array(actVals) ) )

	def testHozDistFixedGroup(self):
		#Create the distribution object; only the toIndices should matter
		currBinResObj = binResHelp.getEmptyBinResultsFromMinMaxAndWidthStandard(0,10,1) #Shouldnt REALLY matter much
		dudIndicesFrom, indicesTo = [0], [4,5,6]
		self.useGroups = [ [1] ] #This has a slight non-zero z-value so...

		currArgs = [currBinResObj, dudIndicesFrom, indicesTo]
		self.distrOptObjs = [distrOptObjHelp.CalcHozDistOptions(*currArgs)]
		self.createTestObjs()

		#Run + check
		expVals = [ (4.47213595499958,), (5.408326913195984,), (4.6097722286464435,) ]
		actVals = self._runTestFunct()


		self.assertTrue( np.allclose( np.array(expVals), np.array(actVals) ) )

	def testMinHozDistInterGroup(self):
		#Create the distribution object; only the toIndices should matter
		currBinResObj = binResHelp.getEmptyBinResultsFromMinMaxAndWidthStandard(0,10,1) #Shouldnt REALLY matter much
		dudIndicesFrom, indicesTo = [0], list()
		self.useGroups = [ [1,0] ] #This has a slight non-zero z-value so...

		currArgs = [currBinResObj, dudIndicesFrom, indicesTo]
		self.distrOptObjs = [distrOptObjHelp.CalcHozDistOptions(*currArgs)]
		self.createTestObjs()

		#Run + check
		expVals = [ (5,) ]
		actVals = self._runTestFunct()

		self.assertTrue( np.allclose( np.array(expVals), np.array(actVals) ) )


	def testCountHBondsFixedGroup(self):
		currBinResObj = self.binResObjA #Built for counting h-bonds regardless
		currArgs = [ currBinResObj, self.fromNonHy, self.fromHy, self.toNonHy, self.toHy ]
		currKwargs = {"maxOO":self.maxOOHBond,"maxAngle":self.maxAngleHBond}
		self.distrOptObjs = [ distrOptObjHelp.CountHBondsBetweenGenericGroupsOptions(*currArgs,**currKwargs) ]
		self.useGroups = [ [0] ]
		self.createTestObjs()

		#Run + check
		expVals = [ (1,) ]
		actVals = self._runTestFunct()

		self.assertTrue( np.allclose( np.array(expVals), np.array(actVals) ) )


	def testCountHBondsInterGroup(self):
		currBinResObj = self.binResObjA #Built for counting h-bonds regardless
		currArgs = [ currBinResObj, self.fromNonHy, self.fromHy, self.toNonHy, self.toHy ]
		currKwargs = {"maxOO":self.maxOOHBond,"maxAngle":self.maxAngleHBond}
		self.distrOptObjs = [ distrOptObjHelp.CountHBondsBetweenGenericGroupsOptions(*currArgs,**currKwargs) ]
		self.useGroups = [ [0,1] ]
		self.createTestObjs()

		#Run + check
		expVals = [ (0,) ]
		actVals = self._runTestFunct()

		self.assertTrue( np.allclose( np.array(expVals), np.array(actVals) ) )

	def testExpectedNumbNearestHozNebs_valsA(self):
		currBinResObj = binResHelp.getEmptyBinResultsFromMinMaxAndWidthStandard(0,10,1) #Shouldnt REALLY matter much
		dudIndicesFrom, indicesTo = [0], list() 
		self.useGroups = [ [0,1] ] #Hydroxyl to hydroxyl

		#
		currArgs = [currBinResObj, dudIndicesFrom, indicesTo]
		hozSepDistrOpts = distrOptObjHelp.CalcHozDistOptions(*currArgs)
		distRange = [0.1,3]
		useDistrOpts = distrOptObjHelp.CountNWithinDistOptions(currBinResObj, hozSepDistrOpts, distRange)
		self.distrOptObjs = [useDistrOpts]
		self.createTestObjs()

		#
		expVals = [(0,)]
		actVals = self._runTestFunct()

		self.assertTrue( np.allclose( np.array(expVals), np.array(actVals) ) )

	def testExpectedNumbNearestHozNebs_valsB(self):
		currBinResObj =  binResHelp.getEmptyBinResultsFromMinMaxAndWidthStandard(0,10,1) #Shouldnt REALLY matter much
		dudIndicesFrom, indicesTo = [0], [4,7] #Only second val matters; these are the non-hy indices we look towards 
		self.useGroups = [ [0] ]

		#
		currArgs = [currBinResObj, dudIndicesFrom, indicesTo]
		hozSepDistrOpts = distrOptObjHelp.CalcHozDistOptions(*currArgs)
		distRange = [0.1,3.1]
		useDistrOpts = distrOptObjHelp.CountNWithinDistOptions(currBinResObj, hozSepDistrOpts, distRange)
		self.distrOptObjs = [useDistrOpts]
		self.createTestObjs()

		#
		expVals = [(1,)]
		actVals = self._runTestFunct()

		self.assertTrue( np.allclose( np.array(expVals), np.array(actVals) ) )


class TestGenericNonHyAndHyFilteredAtomComboDistribs_varyDistrOpts(unittest.TestCase):

	def setUp(self):
		#Define the geometry; 2 "free" water/hydroxyl and 2 hydrogen bonded to each other
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		self.hydroxylA = [ [0,0,0,"O"], [0,1,0.2,"H"] ]
		self.hydroxylB = [ [5,0,0.1,"O"], [5,1,0.3,"H"] ]
		self.waterA = [ [1,2,0,"O"], [0.5,3,0,"H"], [1.5,3,0,"H"] ]
		self.waterB = [ [6,6,0,"O"], [5.5,7,0,"H"], [6.5,7,0,"H"] ]

		self.cartCoords = self.hydroxylA + self.hydroxylB + self.waterA + self.waterB

		#Options for the classifier
		self.binResObjA = binResHelp.BinnedResultsStandard.fromBinEdges([-0.5,0.5,1.5,2.5,3.5,4.5])
		self.fromNonHy = [ [0], [2] ]
		self.fromHy    = [ [1], [3] ]
		self.toNonHy = [ [4], [7] ]
		self.toHy    = [ [5,6],[8,9] ]

		#Sort out filter criterion
		self.nTotalFilterRanges = [ [-0.5,100], [10,100] ] #Any number of h-bonds allowed for first; insane number for the 2nd (so all hydroxyl are in group 0)
		self.maxOOHBond = 3
		self.maxAngleHBond = 50

		#Create the distribution object; which will change in each test but...
		currArgs = [ self.binResObjA, [x[0] for x in self.toNonHy] ] #Neither arg should actually matter; hence I've set the wrong value for the 2nd arg on purpose
		self.distrOptObjs = distrOptObjHelp.CalcPlanarDistOptions(*currArgs, planeEqn=planeEqnHelp.ThreeDimPlaneEquation(0,0,1,4))
		self.useGroups = [ [0] ]
		self.useNonHyIdx, self.useIdxEach = True, 0

		self.createTestObjs()

	def createTestObjs(self):
		#Create the geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

		#Create the classifier options object
		currArgs = [ [self.binResObjA, self.binResObjA], self.fromNonHy, self.fromHy, self.toNonHy, self.toHy ]
		currKwargs = {"nTotalFilterRanges":self.nTotalFilterRanges, "maxOOHBond":self.maxOOHBond, "maxAngleHBond":self.maxAngleHBond}
		self.classifierOptsObj = clsDistrOptObjs.ClassifyBasedOnHBondingToGroup_simple(*currArgs, **currKwargs)

		#Create the main options object
		currArgs = [self.fromNonHy, self.fromHy, self.classifierOptsObj, self.distrOptObjs, self.useGroups]
		currKwargs = {"useNonHyIdx":self.useNonHyIdx, "useIdxEach":self.useIdxEach}
		self.optsObjFilteredCombo = filteredComboOptObjHelp.GenericNonHyAndHyFilteredOptsObj_getAverageVal(*currArgs, **currKwargs)

		#Create the sparse matrix calculator and populate it
		self.sparseMatrixCalculator = optsObjMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.optsObjFilteredCombo])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

		#Create the binval getter
		self.binValGetter = optsObjMapHelp.getMultiDimBinValGetterFromOptsObjs([self.optsObjFilteredCombo])

	def _runTestFunct(self):
		return self.binValGetter.getValsToBin(self.sparseMatrixCalculator)

	def _checkExpAndActValsEqual(self, expVals, actVals):
		#Check we really have 1 value in each
		self.assertTrue(len(expVals)==1)
		self.assertTrue(len(expVals[0])==1)
		self.assertEqual(len(expVals),len(actVals))
		self.assertEqual(len(expVals[0]),len(actVals[0]))

		#Then check we get the expected
		self.assertAlmostEqual(expVals[0][0],actVals[0][0])

	def testExpectedPlanarDistrA_singleValue(self):
		#Only 1st hydroxyl in 1st group; 2nd hydroxyl in the other
		self.nTotalFilterRanges = [ [0.5,100], [-0.5,0.5] ] #>0 h-bonds and 0 hbonds are the criteria
		self.createTestObjs()
		expVals = [ (4,) ]
		actVals = self._runTestFunct()
		self._checkExpAndActValsEqual(expVals,actVals)

	def testExpectedPlanarDistrA_twoValues(self):
		expVals = [ (3.95,) ]
		actVals = self._runTestFunct()
		self._checkExpAndActValsEqual(expVals,actVals)

	def testExpectedTotalHBonds(self):
		#Setup the h-bond counting distr opts
		dudIndices = [ [1] ], [ [1] ] #Use these for values that SHOULD get overwritten
		currArgs = [ self.binResObjA, dudIndices, dudIndices, self.toNonHy, self.toHy]
		currKwargs = {"acceptor":True, "donor":True}
		self.distrOptObjs = distrOptObjHelp.CountHBondsBetweenGenericGroupsOptions(*currArgs, **currKwargs)
		self.createTestObjs()

		#Sort the test out
		expVals = [ (0.5,) ]
		actVals = self._runTestFunct()
		self._checkExpAndActValsEqual(expVals,actVals)

class TestHBondedToDynamicGroup_varyDistrOpts(unittest.TestCase):

	def setUp(self):
		#The geometry(left OH donates 1-hbond to middle water; middle water donates one to right-hand water)
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.hydroxylA =  [ [2,2,2,"O"], [3,2,2,"H"] ]
		self.waterA    = [ [4,2,2.1,"O"], [5,2,2,"H"], [3.5,1.5,2,"H"] ]
		self.waterB    = [ [6,2,2.2,"O"], [7,2,2,"H"], [5.5,1.5,2,"H"] ]

		self.cartCoords = self.hydroxylA + self.waterA + self.waterB

		#Options for the individual classifiers [some options shared]
		self.binResObjA = binResHelp.BinnedResultsStandard.fromBinEdges([-0.5,0.5,1.5,2.5,3.5,4.5])
		self.fromNonHyIndices, self.fromHyIndices = [ [2], [5] ], [ [3,4], [6,7]  ] #Used in both the objects
		self.toNonHyIndices, self.toHyIndices = [ [0] ], [ [1] ]
		self.maxOOHBond = 2.4
		self.maxAngleHBond = 30 #Angle doesnt matter in original 3-molecule geom
		self.nTotalFilterRangesStatic  = [ [-1,0.1], [0.1,1.1] ] #Second one filters to 1 h-bond to hydroxyl 
		self.nTotalFilterRangesDynamic = [ [-1,0.1], [0.1,1.1] ] #Second group filters to 1 h-bond to atoms with 1 h-bond to hydroxyl

		self.dynToNonHyIndices = copy.deepcopy(self.fromNonHyIndices)
		self.dynToHyIndices = copy.deepcopy(self.fromHyIndices)

		#Any extra options for the combined classifier options object
		self.mutuallyExclusive = True
		self.firstClassifierObjs = None

		#Create a distribution object; will likely change each test but...
		currArgs = [ self.binResObjA, [x[0] for x in self.toNonHyIndices] ] #Neither arg should actually matter; hence I've set the wrong value for the 2nd arg on purpose
		self.distrOptObjs = [distrOptObjHelp.CalcPlanarDistOptions(*currArgs, planeEqn=planeEqnHelp.ThreeDimPlaneEquation(0,0,1,4))]

		#Some final options
		self.useGroups = [ [1] ]
		self.useNonHyIdx, self.useIdxEach = True, 0

		self.createTestObjs()

	def createTestObjs(self):
		#Create the geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

		#Create the first classifier options object [For the group with static indices]
		currArgs = [ [self.binResObjA, self.binResObjA], self.fromNonHyIndices, self.fromHyIndices, self.toNonHyIndices, self.toHyIndices]
		currKwargs = {"nTotalFilterRanges":self.nTotalFilterRangesStatic, "maxOOHBond":self.maxOOHBond,
		              "maxAngleHBond":self.maxAngleHBond}
		self.staticGroupOptObj = clsDistrOptObjs.ClassifyBasedOnHBondingToGroup_simple(*currArgs, **currKwargs)

		#Classifier options for the second group (which is based on h-bond relationships to the first)
		#Note we pass the .fromIndices of classifier 1 as "toIndices" here; since they ARE our target
		currArgs = [ [self.binResObjA, self.binResObjA], self.fromNonHyIndices, self.fromHyIndices, self.dynToNonHyIndices, self.dynToHyIndices]
		currKwargs = {"nTotalFilterRanges":self.nTotalFilterRangesDynamic, "maxOOHBond":self.maxOOHBond,
		              "maxAngleHBond":self.maxAngleHBond}
		self.dynamicGroupOptObj = clsDistrOptObjs.ClassifyBasedOnHBondingToGroup_simple(*currArgs, **currKwargs)

		#Overall classifier options object
		currArgs = [self.staticGroupOptObj, self.dynamicGroupOptObj]
		currKwargs = {"mutuallyExclusive":self.mutuallyExclusive, "firstClassifierObjs":self.firstClassifierObjs}
		self.classifierOptObj = clsDistrOptObjs.ClassifyBasedOnHBondingToDynamicGroup(*currArgs,**currKwargs)

		#Create the main options object
		currArgs = [self.fromNonHyIndices, self.fromHyIndices, self.classifierOptObj, self.distrOptObjs, self.useGroups]
		currKwargs = {"useNonHyIdx":self.useNonHyIdx, "useIdxEach":self.useIdxEach}
		self.optsObjFilteredCombo = filteredComboOptObjHelp.GenericNonHyAndHyFilteredOptsObj_simple(*currArgs, **currKwargs)

		#Create the sparse matrix calculator and populate it
		self.sparseMatrixCalculator = optsObjMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.optsObjFilteredCombo])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

		#Create the binval getter
		self.binValGetter = optsObjMapHelp.getMultiDimBinValGetterFromOptsObjs([self.optsObjFilteredCombo])

	def _runTestFunct(self):
		return self.binValGetter.getValsToBin(self.sparseMatrixCalculator)

	def testExpectedPlanarA(self):
		expVals = [ (1.8,) ]
		actVals = self._runTestFunct()
		self.assertTrue( np.allclose( np.array(expVals), np.array(actVals) ) )



class TestHydroxylDiatomFromNonHyAndHyFilteredVaryDistr(unittest.TestCase):

	def setUp(self):
		#All geometric params
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.hydroxylA = [ [0,0,0,"O"], [0,0,1,"H"] ] #1 length; 0 degree angle
		self.hydroxylB = [ [3,0,0,"O"], [3,1,1,"H"] ] #sqrt(2) length; 45 degree angle
		self.hydroxylC = [ [5,0,0,"O"], [7,0,0,"H"] ] #2 length;90 degrees angle

		self.coords = self.hydroxylA + self.hydroxylB + self.hydroxylC

		#Sort some classifiers using min-hoz distance [since its sort of simplest
		self.binResObjA = binResHelp.BinnedResultsStandard.fromBinEdges([0,4,10])
		self.fromNonHyIndices, self.fromHyIndices = [ [0], [2], [4] ], [ [1], [3], [5] ]
		self.toNonHyIndices, self.toHyIndices = [ [0] ], [ [1] ]
		self.distFilterRanges = [ [3.5,10] ] #Should leave only ONE hydroxyl to look at 
		self.useIndicesFrom, self.useIndicesTo = "nonHy", "nonHy"


		#Create a distribution object; will likely change each test but...
		currArgs = [ self.binResObjA, [ [0,1] ] ] #Neither arg should actually matter; hence I've set the wrong value for the 2nd arg on purpose
		self.distrOptObjs = [distrOptObjHelp.GetDistsForDiatomOpts(*currArgs)]

		#Some final options
		self.useGroups = [ [0] ]
		self.useNonHyIdx = True
		self.createTestObjs()

		self.createTestObjs()

	def createTestObjs(self):
		#Create the geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

		#Create the Classifier opts obj
		currArgs = [self.binResObjA, self.fromNonHyIndices, self.fromHyIndices, self.toNonHyIndices, self.toHyIndices, self.distFilterRanges]
		currKwargs = {"useIndicesFrom":self.useIndicesFrom, "useIndicesTo":self.useIndicesTo}
		self.classifierOptsA = clsDistrOptObjs.ClassifyNonHyAndHyBasedOnMinHozDistsToAtomGroups(*currArgs, **currKwargs)

		#Create the filtered options obj
		currArgs = [self.fromNonHyIndices, self.fromHyIndices, self.classifierOptsA, self.distrOptObjs, self.useGroups]
		currKwargs = {"useNonHyIdx":self.useNonHyIdx}
		self.filteredOptsObj = filteredComboOptObjHelp.HydroxylDiatomFromNonHyAndHyFilteredOptsObj_simple(*currArgs, **currKwargs)

		#Create the sparse matrix calculator and populate it
		self.sparseMatrixCalculator = optsObjMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.filteredOptsObj])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

		#Create the binval getter
		self.binValGetter = optsObjMapHelp.getMultiDimBinValGetterFromOptsObjs([self.filteredOptsObj])

	def _runTestFunct(self):
		return self.binValGetter.getValsToBin(self.sparseMatrixCalculator)

	#TODO: Use the "HydroxylDiatomFromNonHyAndHyFilteredOptsObj_simple" object in filtered_atom_combo_opt_objs
	def testExpectedBondLengthsA(self):
		expVals = [ (2,) ]
		actVals = self._runTestFunct()
		self.assertTrue( np.allclose( np.array(expVals), np.array(actVals) ) )

	def testExpectedBonLengthsB(self):
		self.distFilterRanges = [ [-0.1,3.5] ]
		self.createTestObjs()
		expVals = [(1,),(math.sqrt(2),)]
		actVals = self._runTestFunct()
		self.assertTrue( np.allclose( np.array(expVals), np.array(actVals) ) )

	def testExpectedAnglesA(self):
		currArgs = [ self.binResObjA, [ [0,1] ] ] #Neither arg should actually matter; hence I've set the wrong value for the 2nd arg on purpose
		inpVector = [0,0,1]
		self.distrOptObjs = [distrOptObjHelp.GetAngleWithGenericVectorForDiatomOpts(*currArgs,inpVector)]
		self.createTestObjs()

		expVals = [(90,)]
		actVals = self._runTestFunct()
		self.assertTrue( np.allclose( np.array(expVals), np.array(actVals) ) )


	def testExpectedHozDistsA(self):
		self.coords[-1][-2] += 1 #Should only change the absolute(not horizontal) distance
		self.distrOptObjs = [distrOptObjHelp.GetHozDistsForDiatomOpts( self.binResObjA, [ [0,1] ] )] #Neither arg should actually matter; hence I've set the wrong value for the 2nd arg on purpose
		self.createTestObjs()

		expVals = [ (2,) ]
		actVals = self._runTestFunct()
		self.assertTrue( np.allclose( np.array(expVals), np.array(actVals) ) )


class TestWaterDerivativeFilteredOptsDistr(unittest.TestCase):

	def setUp(self):
		#Define the geometry; 2 "free" water/hydroxyl and 2 hydrogen bonded to each other
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.waterA = [ [0,0,0,"O"], [0,0.5,0,"H"], [0,-0.5,0,"H"] ]
		self.hydroxylA = [ [2,2,4,"O"], [2,2,5,"H"] ]
		self.hydroxylB = [ [5,5,6,"O"], [5,5,7,"H"] ]

		self.cartCoords = self.waterA + self.hydroxylA + self.hydroxylB

		#Options for the classifier
		self.binResObjA = binResHelp.BinnedResultsStandard.fromBinEdges([-0.5,0.5,1.5,2.5,3.5,4.5])
		self.oxyIndices = [0,3,5]
		self.hyIndices = [1,2,4,6]
		self.maxOHDist = 1.5
		self.nNebs = [1,2]

		#Options for the planar distribution
		currArgs = [ self.binResObjA, self.hyIndices ] #Neither arg should actually matter; hence I've set the wrong value for the 2nd arg on purpose
		self.distrOptObjs = [distrOptObjHelp.CalcPlanarDistOptions(*currArgs, planeEqn=planeEqnHelp.ThreeDimPlaneEquation(0,0,1,1))]

		#Extra Options for the filtered Distr options object
		self.useGroups = [ [0] ]
		self.useNonHyIdx, self.useIdxEach = True, 0

		self.createTestObjs()
		

	def createTestObjs(self):
		#Create the geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

		#Create the classifier options object
		currArgs = [ [self.binResObjA, self.binResObjA], self.oxyIndices, self.hyIndices ]
		currKwargs = {"maxOHDist":self.maxOHDist, "nNebs":self.nNebs}
		self.classifierOptsObj = clsDistrOptObjs.WaterDerivativeBasedOnDistanceClassifierOptsObj(*currArgs, **currKwargs)

		#Create the main options object
		currArgs = [ self.oxyIndices, self.hyIndices, self.classifierOptsObj, self.distrOptObjs, self.useGroups]
		currKwargs = {"useNonHyIdx":self.useNonHyIdx, "useIdxEach":self.useIdxEach}
		self.optsObjFilteredCombo = filteredComboOptObjHelp.WaterDerivativeFilteredOptsObj_simple(*currArgs, **currKwargs)

		#Create the sparse matrix calculator and populate it
		self.sparseMatrixCalculator = optsObjMapHelp.getSparseMatrixCalculatorFromOptsObjIter([self.optsObjFilteredCombo])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

		#Create the binval getter
		self.binValGetter = optsObjMapHelp.getMultiDimBinValGetterFromOptsObjs([self.optsObjFilteredCombo])

	def _runTestFunct(self):
		return self.binValGetter.getValsToBin(self.sparseMatrixCalculator)

	def testExpectedHydroxylPlanarDistr(self):
		expVals = [ (3,), (5,) ]
		actVals = self._runTestFunct()
		self.assertTrue( np.allclose( np.array(expVals), np.array(actVals) ) )

