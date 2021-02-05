

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

class TestMetaDynOptsObject(unittest.TestCase):

	def setUp(self):
		self.metaVars = mock.Mock()
		self.ntHills = 2
		self.doHills = False
		self.hillHeight = 3
		self.heightNumbFmt = "{:.2f}".format(self.hillHeight)
		self.printColvarCommonIter = 4
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"ntHills":self.ntHills, "doHills":self.doHills, "hillHeight":self.hillHeight,
		             "heightNumbFmt":self.heightNumbFmt,
		             "printColvarCommonIter":self.printColvarCommonIter}
		self.testObjA = tCode.MetaDynamicsOptsCP2KStandard(self.metaVars, **kwargDict)

	def testExpectedKeysPresentInOptDictA_relevantOptsSet(self):
		expDict = {"metaDyn_doHills":self.doHills, "metaDyn_ntHills":self.ntHills,
		           "metaDyn_hillHeight": self.heightNumbFmt.format(self.hillHeight),
		           "metaDyn_printColvarCommonIterLevels": self.printColvarCommonIter,
		           "metaVars":self.metaVars}

		actDict = self.testObjA.optDict

		for key in expDict.keys():
			self.assertEqual( expDict[key], actDict[key] )


	def testExpectedKeysPresentInOptDictB_someRelevantOptsNone(self):
		self.ntHills = None
		self.createTestObjs()
		expMissingKeys = ["metaDyn_ntHills"]
		expDict = {"metaDyn_doHills":self.doHills, 
		           "metaDyn_hillHeight": self.heightNumbFmt.format(self.hillHeight),
		           "metaDyn_printColvarCommonIterLevels": self.printColvarCommonIter,
		           "metaVars":self.metaVars}

		actDict = self.testObjA.optDict

		for key in expDict.keys():
			self.assertEqual( expDict[key], actDict[key] )

		for key in expMissingKeys:
			self.assertTrue( key not in actDict.keys() )



