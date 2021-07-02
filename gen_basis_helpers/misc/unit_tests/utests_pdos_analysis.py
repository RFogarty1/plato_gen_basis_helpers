
import copy
import itertools as it
import unittest
import unittest.mock as mock

import gen_basis_helpers.cp2k.parse_pdos_files as parsePdosHelp
import gen_basis_helpers.misc.pdos_analysis as tCode

class TestGetSpecFragments(unittest.TestCase):

	def setUp(self):
		#Options for parsed Pdos
		self.eigenValues = [2,3]
		self.occs = [1,2]
		self.breakdowns = [ [4,5,6],
		                    [7,8,9] ]
		self.breakdownHeaders = ["s","p","d"]

		#Option for calling the function
		self.fragName = "fragA"
		self.multByOcc = False
		self.separateShells = False
		self.createTestObjs()

	def createTestObjs(self):
		currKwargs = {"eigenValues":self.eigenValues, "occs":self.occs, "fragName":self.fragName,
		              "breakdowns": self.breakdowns, "breakdownHeaders": self.breakdownHeaders}
		self.testObj = parsePdosHelp.PdosFragmentStandard(**currKwargs)
		self.expLabelMergedShells = tCode.simXpsStdObjsHelp.MolFragLabel(fragKey=self.fragName, eleKey="",aoKey="")
		self.mockRetValA = mock.Mock()


	def _runTestFunct(self):
		currKwargs = {"fragName":self.fragName, "multByOcc":self.multByOcc, "separateShells":self.separateShells}
		return tCode.getSimXpsSpecFragmentsFromPdosFragmentStandard(self.testObj, **currKwargs)
#specCreatorHelp

	@mock.patch("gen_basis_helpers.misc.pdos_analysis.specCreatorHelp.SpectrumFragmentStandard")
	def testExpectedCallsMadeA(self, mockedSpecFragObj):
		expIntensities = [ sum(self.breakdowns[0]), sum(self.breakdowns[1]) ]
		expEnergies = self.eigenValues
		mockedSpecFragObj.side_effect = lambda *args,**kwargs: self.mockRetValA		

		actRetVals = self._runTestFunct()
		mockedSpecFragObj.assert_called_with(expEnergies, expIntensities, self.expLabelMergedShells)
		self.assertEqual( 1,len(actRetVals) )
		actRetVal = actRetVals[0]
		self.assertEqual(self.mockRetValA, actRetVal)

	@mock.patch("gen_basis_helpers.misc.pdos_analysis.specCreatorHelp.SpectrumFragmentStandard")
	def testExpectedCallsMade_multByOccs(self, mockedSpecFragObj):
		self.multByOcc = True
		self.createTestObjs()

		expEnergies = self.eigenValues
		unWeightedIntensities = [ sum(self.breakdowns[0]), sum(self.breakdowns[1]) ]
		expIntensities = [occ*i for occ,i in it.zip_longest(self.occs,unWeightedIntensities)]
		mockedSpecFragObj.side_effect = lambda *args,**kwargs: self.mockRetValA		

		actRetVals = self._runTestFunct()
		mockedSpecFragObj.assert_called_with(expEnergies, expIntensities, self.expLabelMergedShells)
		self.assertEqual( 1,len(actRetVals) )
		actRetVal = actRetVals[0]
		self.assertEqual(self.mockRetValA, actRetVal)

	@mock.patch("gen_basis_helpers.misc.pdos_analysis.specCreatorHelp.SpectrumFragmentStandard")
	def testExpectedCallsMade_separateShells(self, mockedSpecFragObj):
		self.separateShells = True
		self.createTestObjs()

		expEnergies = self.eigenValues
		expIntensities = [ [4,7], [5,8], [6,9] ]
		expLabels = [ tCode.simXpsStdObjsHelp.MolFragLabel(fragKey=self.fragName, eleKey="",aoKey=aoKey) for aoKey in self.breakdownHeaders ]

		retVals = self._runTestFunct()
		for intensities,label in it.zip_longest(expIntensities,expLabels):
			mockedSpecFragObj.assert_any_call(expEnergies, intensities, label)
		self.assertEqual( len(expIntensities), len(retVals) )



class TestGetSummedPdosFragments(unittest.TestCase):

	def setUp(self):
		#Setup fragment args
		self.eigenValsA = [-8, -6, 2]
		self.eigenValsB = copy.deepcopy(self.eigenValsA)

		self.occsA = [2,1,0]
		self.occsB = copy.deepcopy(self.occsA)

		self.fragNameA = "frag_a"
		self.fragNameB = "frag_b"

		self.breakdownHeadersA = ["s","p"]
		self.breakdownHeadersB = ["s","p"]

		#
		self.breakdownsA = [ [0.5, 0.1], [0.3,0.2], [0.1,0.3] ]
		self.breakdownsB = [ [0.3, 0.1], [0.2,0.5], [0.2,0.2] ]

		#Options to pass to the run function
		self.outFragName = None

		self.createTestObjs()

	def createTestObjs(self):

		kwargsA = {"eigenValues":self.eigenValsA, "occs":self.occsA, "fragName":self.fragNameA,
		           "breakdowns":self.breakdownsA, "breakdownHeaders":self.breakdownHeadersA}
		kwargsB = {"eigenValues":self.eigenValsB, "occs":self.occsB, "fragName":self.fragNameB,
		           "breakdowns":self.breakdownsB, "breakdownHeaders":self.breakdownHeadersB}

		self.testObjA = parsePdosHelp.PdosFragmentStandard(**kwargsA)
		self.testObjB = parsePdosHelp.PdosFragmentStandard(**kwargsB)

		self.pdosFragments = [self.testObjA, self.testObjB]

	def _runTestFunct(self):
		args = [self.pdosFragments]
		kwargs = {"outFragName": self.outFragName}
		return tCode.getSummedPdosFragments_singleElement(*args, **kwargs)

	def testExpected_validInput(self):
		#Create the expected object
		expBreakdowns = [ [0.8, 0.2], [0.5, 0.7], [0.3, 0.5] ]
		expKwargs = {"eigenValues":self.eigenValsA, "occs":self.occsA, "fragName":self.outFragName,
		             "breakdowns":expBreakdowns, "breakdownHeaders":self.breakdownHeadersA}
		expObj = parsePdosHelp.PdosFragmentStandard(**expKwargs)
		actObj = self._runTestFunct()
		self.assertEqual(expObj,actObj)

	def testRaises_diffHeaders(self):
		self.breakdownHeadersB = [x for x in reversed(self.breakdownHeadersA)]
		self.createTestObjs()
		with self.assertRaises(ValueError):
			self._runTestFunct()

	def testRaises_diffEigenVals(self):
		self.eigenValsA[-1] += 0.5
		self.createTestObjs()
		with self.assertRaises(ValueError):
			self._runTestFunct()


class TestApplyShiftToPdosFrags(unittest.TestCase):

	def setUp(self):
		self.eigenValsA = [2,3,4]
		self.eigenValsB = [3,6]
		self.shiftVal = 5
		self.createTestObjs()

	def _runTestFunct(self):
		args = [self.pdosFragsA, self.shiftVal]
		return tCode.applyEnergyShiftToPdosFragments(*args)

	def createTestObjs(self):
		self.testObjA = parsePdosHelp.PdosFragmentStandard(eigenValues=self.eigenValsA)
		self.testObjB = parsePdosHelp.PdosFragmentStandard(eigenValues=self.eigenValsB)
		self.pdosFragsA = [self.testObjA, self.testObjB]

	def testExpectedCaseA(self):
		eigensA, eigensB = [7,8,9], [8,11]
		expObjA = parsePdosHelp.PdosFragmentStandard(eigenValues=eigensA)
		expObjB = parsePdosHelp.PdosFragmentStandard(eigenValues=eigensB)
		self._runTestFunct()

		actObjA, actObjB = self.testObjA, self.testObjB
		self.assertEqual(expObjA, actObjA)
		self.assertEqual(expObjB, actObjB)


class TestNormalisePdosFrags(unittest.TestCase):

	def setUp(self):
		self.occsA = [5,10,15]
		self.targMaxOcc = 1
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = parsePdosHelp.PdosFragmentStandard(occs=self.occsA)

	def _runTestFunct(self):
		tCode.normalisePdosFragOccupancy([self.testObjA], targMaxOcc=self.targMaxOcc)

	def testExpectedTargetMaxOne(self):
		expOccs = [1/3, 2/3, 1]
		expObj = parsePdosHelp.PdosFragmentStandard(occs=expOccs)
		self._runTestFunct()
		actObj = self.testObjA
		self.assertEqual(expObj, actObj)

	def testExpectedTargetMaxThree(self):
		self.targMaxOcc = 3
		expOccs = [1,2,3]
		expObj = parsePdosHelp.PdosFragmentStandard(occs=expOccs)
		self._runTestFunct()
		actObj = self.testObjA
		self.assertEqual(expObj, actObj)


