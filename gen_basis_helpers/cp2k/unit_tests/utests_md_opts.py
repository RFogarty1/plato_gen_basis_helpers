

import unittest
import unittest.mock as mock

import gen_basis_helpers.cp2k.cp2k_md_options as tCode

class TestStandardMDOptsObject(unittest.TestCase):

	def setUp(self):
		self.timeStep = 0.5
		self.nSteps = 50
		self.ensemble = "NVE"
		self.temperature = 300.5
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"timeStep":self.timeStep, "nSteps":self.nSteps,
		             "ensemble":self.ensemble, "temperature":self.temperature}
		self.testObjA = tCode.MolDynamicsOptsCP2KStandard(**kwargDict)

	def testExpectedKeysPresentInOptDictA_relevantOptsSet(self):
		expDict = {"mdTimeStep":"{:.2f}".format(self.timeStep),
		           "mdSteps":self.nSteps, "mdEnsemble":self.ensemble,
		            "mdTemperature": "{:.2f}".format(self.temperature) } 
		actDict = self.testObjA.optDict

		for key in expDict.keys():
			self.assertEqual( expDict[key], actDict[key] )

	def testExpectedKeysPresentInOptDictB_someRelevantOptsNone(self):
		self.timeStep = None
		self.createTestObjs()

		expMissingKeys = ["mdTimeStep"]
		expDict =  {"mdSteps":self.nSteps, "mdEnsemble":self.ensemble,
		            "mdTemperature": "{:.2f}".format(self.temperature)} 
		actDict = self.testObjA.optDict

		for key in expDict.keys():
			self.assertEqual( expDict[key], actDict[key] )

		for key in expMissingKeys:
			self.assertTrue( key not in actDict.keys() )


