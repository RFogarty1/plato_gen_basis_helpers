

import unittest
import unittest.mock as mock


import gen_basis_helpers.cp2k.method_register as methReg
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



class TestThermostatOpts(unittest.TestCase):

	def setUp(self):
		self.length = 8
		self.mts = 4
		self.timeCon = 20
		self.timeConFmt = "{:.2f}"
		self.yoshida = None
		self.createTestObjs()

	def createTestObjs(self):
		self.pycp2kObj = methReg.createCP2KObjFromMethodStr("cp2k_test_object") 
		kwargs = {"length":self.length, "mts":self.mts, "yoshida":self.yoshida,
		          "timeCon":self.timeCon, "timeConFmt":self.timeConFmt}
		self.testObjA = tCode.NoseThermostatOpts(**kwargs)

	def testExpectedModsToPycp2kObj(self):
		thermostatSection = self.pycp2kObj.CP2K_INPUT.MOTION.MD.THERMOSTAT
		expType, expLength, expMts = "nose", self.length, self.mts
		expTimeCon = self.timeConFmt.format(self.timeCon)
		self.testObjA.addToPyCp2kObj(self.pycp2kObj)

		self.assertEqual(expType, thermostatSection.Type)
		self.assertEqual(expLength, thermostatSection.NOSE.Length)
		self.assertEqual(expMts, thermostatSection.NOSE.Mts)
		self.assertEqual(expTimeCon, thermostatSection.NOSE.Timecon)



