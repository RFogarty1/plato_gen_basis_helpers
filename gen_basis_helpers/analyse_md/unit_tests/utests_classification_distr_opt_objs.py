
import copy
import unittest

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.binned_res as binResHelp

import gen_basis_helpers.analyse_md.classification_distr_opt_objs as tCode


#As below; mainly testing alternative constructor
#This class should mostly be covered by the parent class which is extends
class TestWaterClassifierIncludingNearestFilterAtomMinHozDistsFilter(unittest.TestCase):

	#Just c+p from parent class mainly
	def setUp(self):
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.cartCoords = [ [0,0,0,"O"], [0,0,0,"H"], [0,0,0,"H"],
		                    [0,0,0,"O"], [0,0,0,"H"], [0,0,0,"H"],
		                    [0,0,0,"O"], [0,0,0,"H"], [0,0,0,"H"] ]

		self.binResObjA = binResHelp.BinnedResultsStandard.fromBinEdges([1,2])
		self.binResObjB = binResHelp.BinnedResultsStandard.fromBinEdges([3,4,5])

		self.oxyIndices = [0,3,6]
		self.hyIndices  = [ [1,2], [4,5], [7,8] ]

		self.distFilterIndices = [10,110,12]

		self.distFilterRanges = [ [0,3], [3,6.5] ]
		self.nDonorFilterRanges = [ [0.5,1.5], [1.5,2.5] ]
		self.nAcceptorFilterRanges = [ [1.5,2.5], [2.5,3.5] ]
		self.nTotalFilterRanges = [ [3.5,4.5], [5.5,6.5] ]
		self.adsSiteMinHozToOtherAdsSiteRanges = [ [0,4.0], [4.0,6.5] ]

		self.checkInputConsistent = True
		self.maxOOHBond = 3.2
		self.maxAngleHBond = 40.4

		self.createTestObjs()

	def createTestObjs(self):
		#Sort the cell
		self.cellA = uCellHelp.UnitCell(lattParams=[10,10,10],lattAngles=[90,90,90])
		self.cellA.cartCoords = self.cartCoords

		#
		self.binResObjs = [self.binResObjA, self.binResObjB]
		currArgs = [ self.binResObjs, self.oxyIndices, self.hyIndices, self.distFilterIndices, self.distFilterRanges ]
		currKwargs = { "nDonorFilterRanges":self.nDonorFilterRanges, "nAcceptorFilterRanges":self.nAcceptorFilterRanges,
		               "nTotalFilterRanges": self.nTotalFilterRanges,
		               "maxOOHBond":self.maxOOHBond, "maxAngleHBond": self.maxAngleHBond, "checkInputConsistent":self.checkInputConsistent,
		               "adsSiteMinHozToOtherAdsSiteRanges": self.adsSiteMinHozToOtherAdsSiteRanges}

		self.testObj = tCode.WaterAdsorbedClassifier_usingMinHozDistsBetweenAdsorptionSitesOptsObj(*currArgs, **currKwargs)

	def testCmpEqual(self):
		objA = copy.deepcopy(self.testObj)
		self.createTestObjs()
		objB = self.testObj
		self.assertEqual(objA, objB)

	def testCmpUnequal_diffAdsSiteAdsSiteMinHozDistRanges(self):
		objA = copy.deepcopy(self.testObj)
		self.adsSiteMinHozToOtherAdsSiteRanges[-1][-1] += 2
		self.createTestObjs()
		objB = self.testObj
		self.assertNotEqual(objA,objB)

	def testRaises_adsSiteMinHozRangeWrongLen(self):
		self.adsSiteMinHozToOtherAdsSiteRanges.append( [3.4,4.6] )
		with self.assertRaises(ValueError):
			self.createTestObjs()

	def testFromFilterObjs(self):
		self.nTotalFilterRanges = None #Want to test when i dont both setting explicitly
		self.createTestObjs()

		#Create filter objects
		kwargsA = {"distFilterRange": self.distFilterRanges[0], "nDonorFilterRange": self.nDonorFilterRanges[0],
		           "nAcceptorFilterRange": self.nAcceptorFilterRanges[0], "nTotalFilterRange": None,
		           "minHozDistAdsSitesRange": self.adsSiteMinHozToOtherAdsSiteRanges[0]}

		kwargsB = {"distFilterRange": self.distFilterRanges[1], "nDonorFilterRange": self.nDonorFilterRanges[1],
		           "nAcceptorFilterRange": self.nAcceptorFilterRanges[1], "nTotalFilterRange": None,
		           "minHozDistAdsSitesRange": self.adsSiteMinHozToOtherAdsSiteRanges[1]}

		filterObjA = tCode.WaterMinDistHbondsAndMinHozDistAdsSiteFilterObj(**kwargsA)
		filterObjB = tCode.WaterMinDistHbondsAndMinHozDistAdsSiteFilterObj(**kwargsB)

		currArgs = [ self.binResObjs, self.oxyIndices, self.hyIndices, self.distFilterIndices, [filterObjA, filterObjB] ]
		currKwargs = {"maxOOHBond":self.maxOOHBond, "maxAngleHBond": self.maxAngleHBond, "checkInputConsistent":self.checkInputConsistent}
		actObj = tCode.WaterAdsorbedClassifier_usingMinHozDistsBetweenAdsorptionSitesOptsObj.fromFilterObjs(*currArgs, **currKwargs)
		expObj = self.testObj

		self.assertEqual(expObj, actObj)


#Mainly testing the alternative constructor
class TestWaterCountsOptions(unittest.TestCase):

	def setUp(self):

		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.cartCoords = [ [0,0,0,"O"], [0,0,0,"H"], [0,0,0,"H"],
		                    [0,0,0,"O"], [0,0,0,"H"], [0,0,0,"H"],
		                    [0,0,0,"O"], [0,0,0,"H"], [0,0,0,"H"] ]

		self.binResObjA = binResHelp.BinnedResultsStandard.fromBinEdges([1,2])
		self.binResObjB = binResHelp.BinnedResultsStandard.fromBinEdges([3,4,5])

		self.oxyIndices = [0,3,6]
		self.hyIndices  = [ [1,2], [4,5], [7,8] ]

		self.distFilterIndices = [10,110,12]

		self.distFilterRanges = [ [0,3], [3,6.5] ]
		self.nDonorFilterRanges = [ [0.5,1.5], [1.5,2.5] ]
		self.nAcceptorFilterRanges = [ [1.5,2.5], [2.5,3.5] ]
		self.nTotalFilterRanges = [ [3.5,4.5], [5.5,6.5] ]

		self.checkInputConsistent = True
		self.maxOOHBond = 3.2
		self.maxAngleHBond = 40.4

		self.createTestObjs()

	def createTestObjs(self):
		#Sort the cell
		self.cellA = uCellHelp.UnitCell(lattParams=[10,10,10],lattAngles=[90,90,90])
		self.cellA.cartCoords = self.cartCoords

		#
		self.binResObjs = [self.binResObjA, self.binResObjB]
		currArgs = [ self.binResObjs, self.oxyIndices, self.hyIndices, self.distFilterIndices, self.distFilterRanges ]
		currKwargs = { "nDonorFilterRanges":self.nDonorFilterRanges, "nAcceptorFilterRanges":self.nAcceptorFilterRanges,
		               "nTotalFilterRanges": self.nTotalFilterRanges,
		               "maxOOHBond":self.maxOOHBond, "maxAngleHBond": self.maxAngleHBond, "checkInputConsistent":self.checkInputConsistent}
		self.testObj = tCode.WaterCountTypesMinDistAndHBondSimpleOpts(*currArgs, **currKwargs)

	def testCmpEqual(self):
		objA = copy.deepcopy(self.testObj)
		self.createTestObjs()
		objB = self.testObj
		self.assertEqual(objA,objB)

	def testCmpUnequal_diffFilterIndices(self):
		objA = copy.deepcopy(self.testObj)
		self.distFilterIndices[1] += 3
		self.createTestObjs()
		objB = self.testObj
		self.assertNotEqual(objA,objB)

	def testCmpUnequal_diffHBondParams(self):
		objA = copy.deepcopy(self.testObj)
		self.maxOOHBond += 0.4
		self.createTestObjs()
		objB = self.testObj
		self.assertNotEqual(objA,objB)

	def testCmpUnequal_diffLenFilterRanges(self):
		objA = copy.deepcopy(self.testObj)
		self.nDonorFilterRanges.append( [4.5,6.5] )
		self.checkInputConsistent = False
		self.createTestObjs()
		objB = self.testObj
		self.assertNotEqual(objA, objB)

	def testCmpUnequal_diffFilterRanges(self):
		objA = copy.deepcopy(self.testObj)
		self.nDonorFilterRanges[1][0] += 2
		self.createTestObjs()
		objB = self.testObj
		self.assertNotEqual(objA,objB)

	def testAlternativeConstructor_waterIndices(self):
		waterIndices = [ [0,1,2], [3,4,5], [6,7,8] ]
		currArgs = [self.binResObjs, waterIndices, self.cellA, self.distFilterIndices, self.distFilterRanges]
		currKwargs = { "nDonorFilterRanges":self.nDonorFilterRanges, "nAcceptorFilterRanges":self.nAcceptorFilterRanges,
		               "nTotalFilterRanges": self.nTotalFilterRanges, "maxOOHBond":self.maxOOHBond,
		               "maxAngleHBond": self.maxAngleHBond, "checkInputConsistent":self.checkInputConsistent }

		expObj = self.testObj
		actObj = tCode.WaterCountTypesMinDistAndHBondSimpleOpts.fromWaterIndicesAndGeom(*currArgs, **currKwargs)
		self.assertEqual(expObj, actObj)

	def testAlternativeConstructor_fromFilterObjs(self):


		kwargsA = {"distFilterRange":self.distFilterRanges[0], "nDonorFilterRange":self.nDonorFilterRanges[0],
		 "nAcceptorFilterRange":self.nAcceptorFilterRanges[0], "nTotalFilterRange":self.nTotalFilterRanges[0]}

		kwargsB = {"distFilterRange":self.distFilterRanges[1], "nDonorFilterRange":self.nDonorFilterRanges[1],
		 "nAcceptorFilterRange":self.nAcceptorFilterRanges[1], "nTotalFilterRange":self.nTotalFilterRanges[1]}


		filterObjA = tCode.WaterMinDistAndHBondsFilterObj(**kwargsA)
		filterObjB = tCode.WaterMinDistAndHBondsFilterObj(**kwargsB)

		expObj = self.testObj

		currArgs = [self.binResObjs, self.oxyIndices, self.hyIndices, self.distFilterIndices, [filterObjA, filterObjB]]
		currKwargs = {"maxOOHBond":self.maxOOHBond, "maxAngleHBond":self.maxAngleHBond, "checkInputConsistent":self.checkInputConsistent}
		actObj = tCode.WaterCountTypesMinDistAndHBondSimpleOpts.fromFilterObjs(*currArgs, **currKwargs)
		self.assertEqual(expObj, actObj)


	def testAlternativeConstructor_fromFilterObjs_someNone(self):

		kwargsA = {"distFilterRange":self.distFilterRanges[0], "nDonorFilterRange":None,
		 "nAcceptorFilterRange":self.nAcceptorFilterRanges[0], "nTotalFilterRange":None}

		kwargsB = {"distFilterRange":self.distFilterRanges[1], "nDonorFilterRange":self.nDonorFilterRanges[1],
		 "nAcceptorFilterRange":None, "nTotalFilterRange":None}

		filterObjA = tCode.WaterMinDistAndHBondsFilterObj(**kwargsA)
		filterObjB = tCode.WaterMinDistAndHBondsFilterObj(**kwargsB)

		#Create the expected object
		currArgs = [ self.binResObjs, self.oxyIndices, self.hyIndices, self.distFilterIndices, self.distFilterRanges ]
		currKwargs = { "nDonorFilterRanges":[ [-1,1000] , self.nDonorFilterRanges[1]],
		               "nAcceptorFilterRanges":[self.nAcceptorFilterRanges[0],[-1,1000]],
		               "nTotalFilterRanges": None,
		               "maxOOHBond":self.maxOOHBond, "maxAngleHBond": self.maxAngleHBond, "checkInputConsistent":self.checkInputConsistent}
		expObj = tCode.WaterCountTypesMinDistAndHBondSimpleOpts(*currArgs, **currKwargs)

		#Create actual object and test
		currArgs = [self.binResObjs, self.oxyIndices, self.hyIndices, self.distFilterIndices, [filterObjA, filterObjB]]
		currKwargs = {"maxOOHBond":self.maxOOHBond, "maxAngleHBond":self.maxAngleHBond, "checkInputConsistent":self.checkInputConsistent}
		actObj = tCode.WaterCountTypesMinDistAndHBondSimpleOpts.fromFilterObjs(*currArgs, **currKwargs)

		self.assertEqual(expObj, actObj)


	def testCheckInpConsistent_raisesForVaryLenRanges(self):
		self.nDonorFilterRanges.append( [4,1] )
		with self.assertRaises(ValueError):
			self.createTestObjs()


