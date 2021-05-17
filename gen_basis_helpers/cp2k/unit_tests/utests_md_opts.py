

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
		self.printKindTemp = False
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"timeStep":self.timeStep, "nSteps":self.nSteps,
		             "ensemble":self.ensemble, "temperature":self.temperature,
		             "printKindTemp":self.printKindTemp}
		self.testObjA = tCode.MolDynamicsOptsCP2KStandard(**kwargDict)

	def testExpectedKeysPresentInOptDictA_relevantOptsSet(self):
		self.printKindTemp = True
		self.createTestObjs()

		expDict = {"mdTimeStep":"{:.2f}".format(self.timeStep),
		           "mdSteps":self.nSteps, "mdEnsemble":self.ensemble,
		            "mdTemperature": "{:.2f}".format(self.temperature),
		            "mdPrintKindTemps":self.printKindTemp } 
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
		self.printHills = True
		self.printColvarCommonIter = 4
		self.printHillsCommonIter = 5
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"ntHills":self.ntHills, "doHills":self.doHills, "hillHeight":self.hillHeight,
		             "heightNumbFmt":self.heightNumbFmt, "printHills":self.printHills,
		             "printColvarCommonIter":self.printColvarCommonIter, "printHillsCommonIter":self.printHillsCommonIter}
		self.testObjA = tCode.MetaDynamicsOptsCP2KStandard(self.metaVars, **kwargDict)

	def testExpectedKeysPresentInOptDictA_relevantOptsSet(self):
		expDict = {"metaDyn_doHills":self.doHills, "metaDyn_ntHills":self.ntHills,
		           "metaDyn_hillHeight": self.heightNumbFmt.format(self.hillHeight),
		           "metaDyn_printColvarCommonIterLevels": self.printColvarCommonIter,
		           "metaVars":self.metaVars, "metaDyn_printHills":self.printHills,
		           "metaDyn_printHillsCommonIterLevels":self.printHillsCommonIter}

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


class TestLangevinThermostatOpts(unittest.TestCase):

	def setUp(self):
		self.gamma = 2
		self.noisyGamma = 3
		self.createTestObjs()

	def createTestObjs(self):
		self.pycp2kObj = methReg.createCP2KObjFromMethodStr("cp2k_test_object")
		kwargs = {"gamma":self.gamma, "noisyGamma":self.noisyGamma}
		self.testObjA = tCode.LangevinThermostatOpts(**kwargs)

	def testExpectedCaseA(self):
		langevinSection = self.pycp2kObj.CP2K_INPUT.MOTION.MD.LANGEVIN
		expGamma, expNoisyGamma = self.gamma, self.noisyGamma
		self.testObjA.addToPyCp2kObj(self.pycp2kObj)

		self.assertEqual(expGamma, langevinSection.Gamma)
		self.assertEqual(expNoisyGamma, langevinSection.Noisygamma)

class TestAdaptiveLangevin(unittest.TestCase):

	def setUp(self):
		self.timeConLangevin = 2
		self.timeConNose = 3
		self.createTestObjs()

	def createTestObjs(self):
		self.pycp2kObj = methReg.createCP2KObjFromMethodStr("cp2k_test_object")
		kwargs = {"timeConLangevin":self.timeConLangevin, "timeConNose":self.timeConNose}
		self.testObjA = tCode.AdaptiveLangevinThermostatOpts(**kwargs)

	def testExpectedCaseA(self):
		expTimeConLangevin, expTimeConNose = self.timeConLangevin, self.timeConNose
		expType = "AD_LANGEVIN"
		thermostatSection = self.pycp2kObj.CP2K_INPUT.MOTION.MD.THERMOSTAT
		self.testObjA.addToPyCp2kObj(self.pycp2kObj)

		self.assertEqual(expType, thermostatSection.Type)
		self.assertEqual(expTimeConLangevin, thermostatSection.AD_LANGEVIN.Timecon_langevin)
		self.assertEqual(expTimeConNose, thermostatSection.AD_LANGEVIN.Timecon_nh)


class TestThermalRegion(unittest.TestCase):

	def setUp(self):
		self.doLangevin = True
		self.atomList = [0,2,4]
		self.expAtomList = [1,3,5] #We're using  zero based numbering; this is the base 1 numbering
		self.noisyGammaRegion = 0.2
		self.baseZeroAtomList = True
		self.temperature = 200
		self.createTestObjs()

	def createTestObjs(self):
		self.pycp2kObj = methReg.createCP2KObjFromMethodStr("cp2k_test_object")
		kwargs = {"atomList": self.atomList, "doLangevin":self.doLangevin, "noisyGamma":self.noisyGammaRegion,
		          "baseZeroAtomList":self.baseZeroAtomList, "temperature":200}
		self.testObj = tCode.ThermalRegion(**kwargs)

	def _runTestFunct(self):
		self.testObj.addToPyCp2kObj( self.pycp2kObj )

	def testExpectedCaseA_baseZeroAtomList(self):
		self._runTestFunct()
		self._testExpectedAllOptionsSet()

	def testExpected_baseOneAtomList(self):
		self.baseZeroAtomList = False
		self.createTestObjs()
		self.expAtomList = [x for x in self.atomList]
		self._runTestFunct()
		self._testExpectedAllOptionsSet()

	def _testExpectedAllOptionsSet(self):
		thermalSection = self.pycp2kObj.CP2K_INPUT.MOTION.MD.THERMAL_REGION
		self.assertEqual( 1,len(thermalSection.DEFINE_REGION_list) )
		self.assertEqual( self.doLangevin, thermalSection.DEFINE_REGION_list[-1].Do_langevin )
		self.assertEqual( self.expAtomList, thermalSection.DEFINE_REGION_list[-1].List )
		self.assertAlmostEqual( self.noisyGammaRegion, thermalSection.DEFINE_REGION_list[-1].Noisy_gamma_region)
		self.assertAlmostEqual( self.temperature, thermalSection.DEFINE_REGION_list[-1].Temperature)

class TestThermostatRegion(unittest.TestCase):

	def setUp(self):
		self.baseZeroAtomList = True
		self.atomList = [3,5,7]
		self.expAtomList = [4,6,8] #We're using zero based numbering; this is the base 1 numbering
		self.createTestObjs()

	def createTestObjs(self):
		currKwargs = {"atomList":self.atomList, "baseZeroAtomList":self.baseZeroAtomList}
		self.testObjA = tCode.ThermostatRegionInfo(**currKwargs)
		self.pycp2kObj = methReg.createCP2KObjFromMethodStr("cp2k_test_object")

	def _runTestFunct(self):
		self.testObjA.addToPyCp2kObj( self.pycp2kObj )

	#TODO: Create a mixin property (OutAtomListFromBaseZeroMixin)
	def testExpectedCalls_baseZeroAtomList(self):
		self._runTestFunct()
		self._testExpectedAllOptionsSet()

	def testExpectedCalls_baseOneAtomList(self):
		self.baseZeroAtomList = False
		self.atomList = self.expAtomList
		self.createTestObjs()
		self._runTestFunct()
		self._testExpectedAllOptionsSet()

	def _testExpectedAllOptionsSet(self):
		thermostatSection = self.pycp2kObj.CP2K_INPUT.MOTION.MD.THERMOSTAT
		self.assertEqual(thermostatSection.DEFINE_REGION_list[-1].List, self.expAtomList)


class TestThermostatMultiRegions(unittest.TestCase):

	def setUp(self):
		self.basicThermoOpts = mock.Mock()
		self.regions = [mock.Mock(), mock.Mock()]
		self.mockPyCp2kObj = mock.Mock()
		self.regionKwarg = "test_val"
		self.createTestObjs()

	def createTestObjs(self):
		args = [self.basicThermoOpts, self.regions]
		kwargs = {"regionKwarg":self.regionKwarg}
		self.testObjA = tCode.ThermostatOptsMultiRegions(*args, **kwargs)

	def _runTestFunct(self):
		self.testObjA.addToPyCp2kObj(self.mockPyCp2kObj)

	def testExpectedCallsA(self):

		self._runTestFunct()

		self.basicThermoOpts.addToPyCp2kObj.assert_called_with( self.mockPyCp2kObj )
		for reg in self.regions:
			reg.addToPyCp2kObj.assert_called_with(self.mockPyCp2kObj)

		self.assertEqual( self.regionKwarg , self.mockPyCp2kObj.CP2K_INPUT.MOTION.MD.THERMOSTAT.Region )


