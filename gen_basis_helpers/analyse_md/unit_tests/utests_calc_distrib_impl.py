
import copy
import itertools as it
import unittest
import unittest.mock as mock

import gen_basis_helpers.analyse_md.binned_res as binResHelp
import gen_basis_helpers.analyse_md.calc_angular_distrib_impl as calcAnglDistribImplHelp
import gen_basis_helpers.analyse_md.calc_distrib_core as calcDistribCoreHelp
import gen_basis_helpers.analyse_md.calc_radial_distrib_impl as calcRadDistribImplHelp
import gen_basis_helpers.analyse_md.traj_core as trajHelp

import gen_basis_helpers.analyse_md.calc_distrib_impl as tCode

class TestAverageDistribFunctForEachInpTraj_rdf(unittest.TestCase):

	def setUp(self):

		#rdf opts/bin options
		self.binEdges = [0,2,4]
		self.indicesA_objA = [1,2]
		self.indicesB_objA = [5,6]

		self.indicesA_objB = [3,4]
		self.indicesB_objB = [7,8]
		self.volume_objA = 50
		self.volume_objB = 100

		#Options for the values INSIDE the bins
		self.rdfValsTrajA_one = [5,10]
		self.rdfValsTrajA_two = [15,25]

		self.rdfValsTrajB_one = [20, 30]
		self.rdfValsTrajB_two = [40, 50]

		#Trajectory input
		self.trajA = mock.Mock()
		self.trajB = mock.Mock()

		self.createTestObjs()

	def createTestObjs(self):
		self.binObjA = binResHelp.BinnedResultsStandard.fromBinEdges(self.binEdges)
		self.binObjB = binResHelp.BinnedResultsStandard.fromBinEdges(self.binEdges)

		self.optObjA = calcDistribCoreHelp.CalcRdfOptions(self.binObjA, self.indicesA_objA, self.indicesB_objA, volume=self.volume_objA)
		self.optObjB = calcDistribCoreHelp.CalcRdfOptions(self.binObjB, self.indicesA_objB, self.indicesB_objB, volume=self.volume_objB)

		self.inpOptions = [self.optObjA, self.optObjB]
		self.inpTrajs = [self.trajA, self.trajB]

	def testCreateOptsToUse(self):
		expOpts = [ [copy.deepcopy(self.optObjA), copy.deepcopy(self.optObjB)],
		            [copy.deepcopy(self.optObjA), copy.deepcopy(self.optObjB)] ]

		actOpts = tCode._getInpOptionsForEachTraj(self.inpOptions, self.inpTrajs)

		for expOptList,actOptList in it.zip_longest(expOpts,actOpts):
			for expObj,actObj in it.zip_longest(expOptList, actOptList):
				self._checkTwoOptObjsEqual(expObj, actObj)


	#All ints here
	def _checkTwoOptObjsEqual(self, objA, objB):
		directCmpKeys = ["binResObj", "distribKey", "indicesA", "indicesB", "volume"]
		for key in directCmpKeys:
			self.assertEqual( getattr(objA,key), getattr(objB,key) )


	@mock.patch("gen_basis_helpers.analyse_md.calc_distrib_impl.calcDistribCoreHelp.populateRdfValsOnOptionObjs")
	def testExpectedCallToPopulate_forRdfObj(self, mockPopulateOpts):
		""" Make sure that, using CalcRdfOptions, we call the expected function to populate them """
		tCode._populateInpOptsUsingTrajs([self.inpOptions, self.inpOptions], self.inpTrajs)
		mockPopulateOpts.assert_called()

	#Set the inpOptions objects to planar rdf; then just make sure the relevant function gets called (dont worry about args)
	@mock.patch("gen_basis_helpers.analyse_md.calc_distrib_impl.calcRadialDistribImplHelp.populatePlanarRdfsFromOptionsObjs")
	def testExpectedCallToPopulate_forPlanarRdfObj(self, mockPopulateOpts):
		#Only care about what type of class these options are
		optsA = calcRadDistribImplHelp.CalcPlanarRdfOptions(self.binObjA, [1])
		optsB = calcRadDistribImplHelp.CalcPlanarRdfOptions(self.binObjA, [1])
		tCode._populateInpOptsUsingTrajs([[optsA, optsB]], self.inpTrajs)
		mockPopulateOpts.assert_called()

	@mock.patch("gen_basis_helpers.analyse_md.calc_angular_distrib_impl.populateInteratomicAngularDistribsFromOptionsObjs")
	def testExpectedCallToPopulate_forInteratomicAngularObj(self, mockPopulateOpts):
		optsA = calcAnglDistribImplHelp.CalcInteratomicAngularDistribOptions(self.binObjA, [ [1,2,3] ])
		optsB = calcAnglDistribImplHelp.CalcInteratomicAngularDistribOptions(self.binObjA, [ [3,4,5] ])
		tCode._populateInpOptsUsingTrajs([[optsA, optsB]], self.inpTrajs)
		mockPopulateOpts.assert_called()

	@mock.patch("gen_basis_helpers.analyse_md.calc_distrib_impl.calcDistribCoreHelp.populateRdfValsOnOptionObjs")
	def testRaisesWhenDiffInstanceTypesPassedToPopulateFunct(self, mockPopulateOpts):
		optsA = self.inpOptions
		optsB = [mock.Mock()]
		with self.assertRaises(AssertionError):
			tCode._populateInpOptsUsingTrajs([optsA,optsB], self.inpTrajs)

	def testExpectedAveragedValsGiven(self):
		#Setup the populated bins
		inpOpts = tCode._getInpOptionsForEachTraj(self.inpOptions,self.inpTrajs) 
		inpOpts[0][0].binResObj.binVals["rdf"] = self.rdfValsTrajA_one
		inpOpts[0][1].binResObj.binVals["rdf"] = self.rdfValsTrajA_two
		inpOpts[1][0].binResObj.binVals["rdf"] = self.rdfValsTrajB_one
		inpOpts[1][1].binResObj.binVals["rdf"] = self.rdfValsTrajB_two

		#Figure out expected
		expValsTrajA = [ [1,(15+5)/2], [3,(10+25)/2] ]
		expValsTrajB = [ [1,(20+40)/2], [3,(30+50)/2] ]
		expVals = [expValsTrajA, expValsTrajB]

		actVals = tCode._getAveragedDistribsFromPopulatedBinsMultiTraj(inpOpts)

		for expValList, actValList in it.zip_longest(expVals, actVals):
			for exp,act in it.zip_longest(expValList, actValList):
				[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp,act)]
 

class TestGetCumulativeAveragesOfBinCentresVsVals(unittest.TestCase):

	def setUp(self):
		self.binsVsValsA = [ [1,3], [2,5], [3, 9] ]
		self.binsVsValsB = [ [1,4], [2,8], [3,12] ]
		self.binsVsValsC = [ [1,8], [2,4], [3,20] ]
		self.createTestObjs()

	def createTestObjs(self):
		self.bins = [self.binsVsValsA, self.binsVsValsB, self.binsVsValsC]

	def _runTestFunct(self):
		return tCode.getCumulativeAverageOfBinCentresVsVals(self.bins)

	def testExpectedValsA(self):
		expBinsA = [ [1,3], [2,5], [3, 9] ]
		expBinsB = [ [1,(3+4)/2], [2, (5+8)/2], [3, (9+12)/2] ]
		expBinsC = [ [1,(3+4+8)/3], [2,(5+8+4)/3], [3,(9+12+20)/3] ]

		expBins = [expBinsA, expBinsB, expBinsC]
		actBins = self._runTestFunct()

		for expBinRes, actBinRes in it.zip_longest(expBins, actBins):
			for expBin, actBin in it.zip_longest(expBinRes, actBinRes):
				[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(expBin,actBin)] #Each is [binCentre,binVal]


