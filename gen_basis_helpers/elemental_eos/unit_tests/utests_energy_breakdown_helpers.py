#!/usr/bin/python3


import types
import unittest
import unittest.mock as mock

import numpy as np



#platoOut.parsePlatoOutFile

import gen_basis_helpers.elemental_eos.eos_energy_breakdown_helpers as tCode

class TestGetE0EnergiesFromWorkFlow(unittest.TestCase):

	def setUp(self):
		self.testWorkFlow = createStubWorkFlowA()

	def runTestFunct(self):
		tCode.volVsE0EnergiesFromWorkFlow(self.testWorkFlow, perAtom=True, rydToEv=False)

	@mock.patch("gen_basis_helpers.elemental_eos.eos_energy_breakdown_helpers.platoOut.parsePlatoOutFile")
	def testVolVsE0EnergiesAttachedToOutput(self, mockedParser):
		fakeEnergies = types.SimpleNamespace( e0Coh=80 )
		fakeUCell = types.SimpleNamespace( volume = 200 )
		mockedParser.return_value = {"energies": fakeEnergies, "numbAtoms": 2, "unitCell":fakeUCell}
		self.runTestFunct()

		expOutArray = np.array(([100,40],[100,40]))
		expOutput = types.SimpleNamespace( volVsE0=expOutArray )

		mockedParser.assert_any_call( "fake_path_a.out" )
		mockedParser.assert_any_call( "fake_path_b.out" )
		actOutArray = self.testWorkFlow.output.volVsE0
		self.assertTrue( np.allclose(actOutArray,expOutArray) )
		



def createStubWorkFlowA():
	outDict = {"_inpFilePaths":["fake_path_a", "fake_path_b"]}
	return types.SimpleNamespace(**outDict)


if __name__ == '__main__':
	unittest.main()

