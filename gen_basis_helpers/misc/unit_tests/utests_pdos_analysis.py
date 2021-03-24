
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


