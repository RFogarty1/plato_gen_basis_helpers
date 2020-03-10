

import unittest
import unittest.mock as mock

import gen_basis_helpers.job_helpers.cp2k.self_defects_help as tCode


class TestSelfDefectStandardInputCreator(unittest.TestCase):

	def setUp(self):
		self.absGridCutoff = mock.Mock()
		self.addedMOsBulk = mock.Mock()
		self.addedMOsDefect = mock.Mock()
		self.basisObjs = mock.Mock()
		self.cp2kMethodStr = mock.Mock()
		self.relGridCutoff = mock.Mock()

		self.expKwargDictBoth = {"basisObjs":self.basisObjs, "methodStr":self.cp2kMethodStr,
		                        "absGridCutoff":self.absGridCutoff, "relGridCutoff":self.relGridCutoff}
		self.createTestObjs()


	def createTestObjs(self):
		kwargDict = {"absGridCutoff":self.absGridCutoff, "addedMOsBulk":self.addedMOsBulk,
		             "addedMOsDefect":self.addedMOsDefect, "basisObjs":self.basisObjs,
		             "cp2kMethodStr":self.cp2kMethodStr, "relGridCutoff":self.relGridCutoff}
		self.testObjA = tCode.SelfDefectStandardInputCreator(**kwargDict)


	def testExpectedKwargsPassedForBulk(self):
		expKwargDictToPass = self.expKwargDictBoth
		expKwargDictToPass.update({"addedMOs":self.addedMOsBulk})
		outCreator = self.testObjA._createCalcObjCreatorBulk()
		for key in expKwargDictToPass.keys():
			self.assertEqual( expKwargDictToPass[key], getattr(outCreator,key) ) 

	def testExpectedKwargsPassedForDefect(self):
		expKwargDictToPass = self.expKwargDictBoth
		expKwargDictToPass.update({"addedMOs":self.addedMOsDefect})
		outCreator = self.testObjA._createCalcObjCreatorDefect()
		for key in expKwargDictToPass.keys():
			self.assertEqual( expKwargDictToPass[key], getattr(outCreator,key) )


