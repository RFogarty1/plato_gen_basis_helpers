#!/usr/bin/python3

import os
import types
import unittest
import unittest.mock as mock


import plato_pylib.plato.mod_plato_inp_files as modInp
import gen_basis_helpers.convergers.integ_grid_conv as tCode



class TestVaryAngularSpacingOptDictMapper(unittest.TestCase):

	def setUp(self):
		self.testObj = tCode.createDft2AngularSpacingOptDictMapper()

	def runTestFunct(self, optDict, varyParam):
		self.testObj.modOptDictWithVariableParam( optDict, varyParam )

	def testExpectedGridOptPresent(self):
		inpParam = 40
		basicOptDict = modInp.getDefOptDict("dft2")
		expectedVal = [basicOptDict["integralmeshspacing"][0], inpParam, inpParam]
		self.assertTrue( expectedVal != basicOptDict["integralmeshspacing"] )
		self.runTestFunct(basicOptDict, inpParam)
		self.assertEqual( expectedVal, basicOptDict["integralmeshspacing"] )	

	def testRaisesErrorIfGridTypeWrong(self):
		inpParam = 40
		basicOptDict = modInp.getDefOptDict("dft2")
		basicOptDict["integralmeshtype"] = "uniform"
		with self.assertRaises(AssertionError):
			self.runTestFunct(basicOptDict, inpParam)

	def testRaisesErrorIfGridTypeMissing(self):
		inpParam = 40
		basicOptDict = dict()
		with self.assertRaises(AssertionError):
			self.runTestFunct(basicOptDict, inpParam)


class TestVaryRadialSpacingOptDictMapper(unittest.TestCase):

	def setUp(self):
		self.testObj = tCode.createDft2RadialSpacingOptDictMapper()

	def runTestFunct(self, optDict, varyParam):
		self.testObj.modOptDictWithVariableParam( optDict, varyParam )

	def testExpectedGridOptPresent(self):
		inpParam = 12
		basicOptDict = modInp.getDefOptDict("dft2")
		basicOptDict["integralmeshspacing"][2] = 8
		expectedVal = [inpParam, basicOptDict["integralmeshspacing"][1], basicOptDict["integralmeshspacing"][2]]
		self.assertTrue( expectedVal != basicOptDict["integralmeshspacing"] )
		self.runTestFunct(basicOptDict, inpParam)
		self.assertEqual( expectedVal, basicOptDict["integralmeshspacing"] )	


class TestVaryKPtOptDictMapper(unittest.TestCase):

	def setUp(self):
		self.testObj = tCode.createKPointsMPGridOptDictMapper()

	def runTestFunct(self, optDict, varyParam):
		self.testObj.modOptDictWithVariableParam( optDict, varyParam )

	def testExpectedValsPresent(self):
		inpParam = [5,5,4]
		basicOptDict = modInp.getDefOptDict("dft2")
		basicOptDict["BlochStates".lower()] = [3,2,1]
		self.runTestFunct(basicOptDict, inpParam)
		self.assertEqual( inpParam, basicOptDict["blochstates"] )

class TestVaryFFTSpacingOptDictMapper(unittest.TestCase):

	def setUp(self):
		self.testObj = tCode.createDftFFTSpacingOptDictMapper()

	def runTestFunct(self, optDict, varyParam):
		self.testObj.modOptDictWithVariableParam(optDict, varyParam)

	def testExpectedGridOptPresent(self):
		inpParam = 0.4
		basicOptDict = modInp.getDefOptDict("dft")
		self.assertNotAlmostEqual( inpParam, basicOptDict["fftGridSpacing".lower()] )
		self.runTestFunct( basicOptDict, inpParam )
		self.assertAlmostEqual( inpParam, basicOptDict["fftGridSpacing".lower()] )

#def createVariableToGridInfoObjStubA():
#	outDict = { "varyType": "angular" }
#	return types.SimpleNamespace( **outDict )
#
#def createMethodObjStubA():
#	outDict = {"runCommFunction": "fake_run_comm"}
#	return types.SimpleNamespace( **outDict )

if __name__ == '__main__':
	unittest.main()

