
import copy
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.binned_res as binResHelp
import gen_basis_helpers.analyse_md.traj_core as trajCoreHelp

import gen_basis_helpers.analyse_md.calc_angular_distrib_impl as tCode


class TestPopulateInteratomicAngularDistribs(unittest.TestCase):

	def setUp(self):
		#The co-ordinates to use
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.cartCoords = [ [1,1,7,"X"],
		                    [1,1,9,"X"],
		                    [1,3,1,"Y"], #135 degrees once PBCs taken into acount
		                    [1,5,9,"Y"], #90 degrees with 0,1,this
		                    [1,3,9,"Y"] ] #90 degrees with 0,1,this


		#Indices for each binner
		self.indicesA = [ [0,1,3], [0,1,4], [0,1,2] ] #90,90,135 degrees
		self.indicesB = [ [0,1,3], [0,1,2] ] #90,135 degrees

		#Bin edges
		self.binEdgesA = [0,60,120,180]
		self.binEdgesB = [0,60,120,180]

		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords
		self.trajStepA = trajCoreHelp.TrajStepFlexible(unitCell=self.cellA)
		self.trajA = trajCoreHelp.TrajectoryInMemory([self.trajStepA])

		#sort the bins
		self.binsA = binResHelp.BinnedResultsStandard.fromBinEdges(self.binEdgesA)
		self.binsB = binResHelp.BinnedResultsStandard.fromBinEdges(self.binEdgesB)

		#Sort the options objects
		self.optObjsA = tCode.CalcInteratomicAngularDistribOptions(self.binsA, self.indicesA)
		self.optObjsB = tCode.CalcInteratomicAngularDistribOptions(self.binsB, self.indicesB)

	def _runTestFunct(self):
		tCode.populateInteratomicAngularDistribsFromOptionsObjs(self.trajA, [self.optObjsA, self.optObjsB])

	def testExpectedValsSingleTrajStep(self):
		#Figure out the expected numbers
		expCountsA, expCountsB = [0,2,1], [0,1,1]
		expPdfA, expPdfB = [0,2/3,1/3], [0,0.5,0.5] #The probability distribution function
		binWidth, domainWidth = 60,180
		pdfToAdf = domainWidth/binWidth
		expAdfA, expAdfB = [0, pdfToAdf*(2/3), pdfToAdf*(1/3)], [0, pdfToAdf*0.5, pdfToAdf*0.5]

		#Generate expected bin objects
		expBinsA, expBinsB = copy.deepcopy(self.binsA), copy.deepcopy(self.binsB)
		expBinsA.binVals["counts"], expBinsB.binVals["counts"] = expCountsA, expCountsB
		expBinsA.binVals["adf"], expBinsB.binVals["adf"] = expAdfA, expAdfB
		expBinsA.binVals["pdf"], expBinsB.binVals["pdf"] = expPdfA, expPdfB

		#Compare actual vs expected bin objects
		self._runTestFunct()
		
		actBinsA, actBinsB = self.binsA, self.binsB

		self.assertEqual(expBinsA, actBinsA)
		self.assertEqual(expBinsB, actBinsB)

	def testExpected_binsNotSpanningFullSampledDomain(self):
		#Change the bins
		self.binEdgesA = [60,120,180]
		self.binEdgesB = [60,120,180]
		self.createTestObjs()

		#Figure out the expected numbers
		expCountsA, expCountsB = [2,1], [1,1]
		expPdfA, expPdfB = [2/3,1/3], [0.5,0.5] #The probability distribution function
		binWidth, domainWidth = 60,180
		pdfToAdf = domainWidth/binWidth
		expAdfA, expAdfB = [pdfToAdf*(2/3), pdfToAdf*(1/3)], [pdfToAdf*0.5, pdfToAdf*0.5]

		#Generate expected bin objects
		expBinsA, expBinsB = copy.deepcopy(self.binsA), copy.deepcopy(self.binsB)
		expBinsA.binVals["counts"], expBinsB.binVals["counts"] = expCountsA, expCountsB
		expBinsA.binVals["adf"], expBinsB.binVals["adf"] = expAdfA, expAdfB
		expBinsA.binVals["pdf"], expBinsB.binVals["pdf"] = expPdfA, expPdfB

		#Compare actual vs expected bin objects
		self._runTestFunct()
		
		actBinsA, actBinsB = self.binsA, self.binsB

		self.assertEqual(expBinsA, actBinsA)
		self.assertEqual(expBinsB, actBinsB)


class TestInteratomicOptsInit(unittest.TestCase):

	def setUp(self):
		self.binEdges = [0,90,180]
		self.indices = [[1,2,3]]
		self.checkEdges = True
		self.createTestObjs()

	def createTestObjs(self):
		self.binResObj = binResHelp.BinnedResultsStandard.fromBinEdges(self.binEdges)

	def _createTestObj(self):
		self.testObj = tCode.CalcInteratomicAngularDistribOptions(self.binResObj, self.indices, checkEdges=self.checkEdges)

	def testRaisesWhenBinsOutsideLowerEndOfDomain(self):
		self.binEdges = [-1,90,180]
		self.createTestObjs()
		with self.assertRaises(ValueError):
			self._createTestObj()

	def testRaisesWhenBinsOutsideUpperEndOfDomain(self):
		self.binEdges = [0,90,181]
		self.createTestObjs()
		with self.assertRaises(ValueError):
			self._createTestObj()

	def testDoesntRaiseWhenCheckBinsFalse_binsOutsideDomain(self):
		self.binEdges = [-1,90,181]
		self.checkEdges = False
		self.createTestObjs()
		self._createTestObj() #Shouldnt raise an error; even though bin edges are outside our domain

class TestInteratomicMultiBinner(unittest.TestCase):

	#Note: I used basically the same values for testing the calcANgles function we use
	def setUp(self):

		#The co-ordinates to use
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.cartCoords = [ [1,1,7,"X"],
		                    [1,1,9,"X"],
		                    [1,3,1,"Y"], #135 degrees once PBCs taken into acount
		                    [1,5,9,"Y"], #90 degrees with 0,1,this
		                    [1,3,9,"Y"] ] #90 degrees with 0,1,this


		#Indices for each binner
		self.indicesA = [ [0,1,3], [0,1,4], [0,1,2] ] #90,90,135 degrees
		self.indicesB = [ [0,1,3], [0,1,2] ] #90,135 degrees

		#Bin edges
		self.binEdgesA = [0,60,120,180]
		self.binEdgesB = [0,60,120,180]

		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords
		self.trajStepA = trajCoreHelp.TrajStepFlexible(unitCell=self.cellA)

		#sort the bins
		self.binsA = binResHelp.BinnedResultsStandard.fromBinEdges(self.binEdgesA)
		self.binsB = binResHelp.BinnedResultsStandard.fromBinEdges(self.binEdgesB)

		#Sort the single binners
		self.binnerA = tCode._InteratomicAngularBinnerFixedIndices(resBins=self.binsA, indices=self.indicesA)
		self.binnerB = tCode._InteratomicAngularBinnerFixedIndices(resBins=self.binsB, indices=self.indicesB)

		#Create the test object
		self.testObj = tCode._InteratomicAngularMultiBinnerFixedIndices([self.binnerA, self.binnerB])

	def _runTestFunct(self):
		self.testObj.updateCountsFromTrajStep(self.trajStepA)

	def testDoubleBinnerCaseA(self):
		expBinA, expBinB = copy.deepcopy(self.binsA), copy.deepcopy(self.binsB)
		expBinA.binVals["counts"], expBinB.binVals["counts"] = [0,2,1], [0,1,1]
		self._runTestFunct()
		actBinA, actBinB = self.binsA, self.binsB

		self.assertEqual(expBinA, actBinA)
		self.assertEqual(expBinB, actBinB)


