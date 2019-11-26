#!/usr/bin/python3


import os
import types
import unittest
import unittest.mock as mock

import gen_basis_helpers.job_utils.surface_energies_helpers as tCode


class TestCreateSurfaceRunner(unittest.TestCase):

	def setUp(self):
		self.energyCohesiveElectronic = 20
		self.energyTotalElectronic = 50
		self.nAtoms = 5
		self.surfArea = 7
		self.energyAttr = "electronicCohesiveE"
		self.createTestObj()
	
	def createTestObj(self):
		self.stubPlatoCalcObj = self._createStubPlatoCalcObj()
		self.testObj = tCode.PlatoSurfRunner(self.stubPlatoCalcObj, self.surfArea, self.energyAttr)

	def testExpectedSurfaceEnergyGiven(self):
		expSurfArea = self.surfArea
		actSurfArea = self.testObj.surfaceArea
		self.assertAlmostEqual(expSurfArea, actSurfArea)

	def testWriteFilesCalledWithoutArgs(self):
		self.testObj.writeFiles()
		self.stubPlatoCalcObj.writeFile.assert_called_once_with()

	def testCorrectWorkFolderObtained(self):
		expWorkFolder = os.path.join("fake","path","to")
		actWorkFolder = self.testObj.workFolder
		self.assertEqual(expWorkFolder, actWorkFolder)

	def testCorrectNAtoms(self):
		expNumbAtoms = self.nAtoms
		actNumbAtoms = self.testObj.nAtoms
		self.assertEqual(expNumbAtoms, actNumbAtoms)

	def testCorrectEnergyPerAtomCohesiveElectronic(self):
		self.energyAttr = "electronicCohesiveE"
		self.createTestObj()
		expEnergyPerAtom = self.energyCohesiveElectronic / self.nAtoms
		actEnergyPerAtom = self.testObj.energiesPerAtom
		self.assertAlmostEqual(expEnergyPerAtom, actEnergyPerAtom)

	def testCorrectEnergyPerAtomTotalElectronic(self):
		self.energyAttr = "electronicTotalE"
		self.createTestObj()
		expEnergyPerAtom = self.energyTotalElectronic / self.nAtoms
		actEnergyPerAtom = self.testObj.energiesPerAtom
		self.assertAlmostEqual(expEnergyPerAtom, actEnergyPerAtom)

	def testGetRunComm(self):
		expRunComm = ["hello"]
		actRunComm = self.testObj.runComm
		self.assertEqual(expRunComm, actRunComm)

	def _createStubPlatoCalcObj(self):
		stubEnergies = types.SimpleNamespace( electronicCohesiveE= self.energyCohesiveElectronic, #Stub for plato_pylib energies object 
		                                      electronicTotalE = self.energyTotalElectronic )
		stubObj = types.SimpleNamespace( getRunComm = mock.Mock(),
		                                 writeFile = mock.Mock(),
		                                 filePath = os.path.join("fake","path","to","file"),
		                                 parseOutFile = lambda: {"energies": stubEnergies, "numbAtoms":self.nAtoms} )

		stubObj.getRunComm.return_value = "hello"

		return stubObj


if __name__ == '__main__':
	unittest.main()

