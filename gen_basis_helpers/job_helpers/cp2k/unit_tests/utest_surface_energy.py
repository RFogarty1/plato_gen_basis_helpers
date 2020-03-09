
import os

import unittest
import unittest.mock as mock

import gen_basis_helpers.job_helpers.cp2k.surface_energy_help as tCode

class TestStandardInputCreator(unittest.TestCase):

	def setUp(self):
		self.absGridCutoff = mock.Mock()
		self.addedMOs = mock.Mock()
		self.basisObjs = mock.Mock()
		self.cp2kMethodStr = "cp2k_test_object" #Points to a method that should never change
		self.relGridCutoff = mock.Mock()
		self.workFolder = "fake_folder"
		self.eleKey = "fake_ele"
		self.methodKey = "fake_method"
		self.structKey = "fake_struct"
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"absGridCutoff":self.absGridCutoff, "addedMOs":self.addedMOs, "basisObjs":self.basisObjs,
		             "cp2kMethodStr":self.cp2kMethodStr, "relGridCutoff":self.relGridCutoff, "baseWorkFolder":self.workFolder,
		             "eleKey":self.eleKey, "methodKey":self.methodKey, "structKey":self.structKey} 
		self.testObjA = tCode.SurfaceEnergyStandardInputCreator(**kwargDict)

	def testFactoryCalledWithMethodStr(self):
		creator = self.testObjA._createCalcObjCreator()
		self.assertEqual(self.cp2kMethodStr,creator.methodStr) 

	def testVariousKwargsSet(self):
		""" Test various keywords are set on the creator object """
		expWorkFolder = os.path.join(self.workFolder, self.eleKey, self.methodKey, self.structKey)
		expKwargs = {"addedMOs":self.addedMOs,"basisObjs":self.basisObjs,
		             "methodStr":self.cp2kMethodStr, "absGridCutoff":self.absGridCutoff,
		             "relGridCutoff":self.relGridCutoff, "workFolder":expWorkFolder}
		outObj = self.testObjA._createCalcObjCreator()
		for key in expKwargs:
			objAttr = getattr(outObj,key)
			self.assertEqual( expKwargs[key], objAttr )


