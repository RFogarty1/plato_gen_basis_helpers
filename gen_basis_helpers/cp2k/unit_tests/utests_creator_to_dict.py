

import types
import unittest
import unittest.mock as mock


import gen_basis_helpers.analyse_md.analyse_metadyn_hills as aMetadynHillsHelp
import gen_basis_helpers.cp2k.cp2k_basis_obj as basisObjHelp
import gen_basis_helpers.cp2k.cp2k_creator as creatorHelp
import gen_basis_helpers.cp2k.cp2k_md_options as mdOptsHelp

import gen_basis_helpers.cp2k.cp2k_creator_to_dict as tCode


class TestGetMetadynHillsObjSimple(unittest.TestCase):

	def setUp(self):
		#Init hills
		self.initScales = [ [1,2], [3,4] ]
		self.initPositions = [ [5,6], [7,8] ]
		self.initHeights = [3,5]
		self.initTimes = [5,10] #Possible for initial to have times saved from a previous MD

		#
		self.spawnedScales = [ [2,3] ]
		self.spawnedPositions = [ [8,9] ]
		self.spawnedHeights = [ 6 ]
		self.spawnedTimes = [20]

		#Run opts
		self.raiseIfNoInitHills = False
		self.createTestObjs()

	def createTestObjs(self):
		#Get the initial hills
		currKwargs = {"scales":self.initScales, "heights":self.initHeights,
		              "positions":self.initPositions}
		self.spawnOptsA = mdOptsHelp.MetadynamicsSpawnHillsOptions.fromIters(**currKwargs)

		#Get the creator obj used
		self.inpMetadynOpts = mdOptsHelp.MetaDynamicsOptsCP2KStandard(None,spawnHillsOpts=self.spawnOptsA)
		self.creatorObjA = creatorHelp.CP2KCalcObjFactoryStandard(metaDynOpts=self.inpMetadynOpts)

		#Get an initial object representing hills spawned 
		currKwargs = {"times":self.initTimes, "scales":self.initScales, "heights": [ [h,h] for h in self.initHeights],
		              "positions":self.initPositions}
		self.initHillInfoDict = aMetadynHillsHelp.MetadynHillsInfo(**currKwargs)

		#Create a mock standard output obj; need this to get the data
		self.spawnedDictsA = {"time": self.spawnedTimes, "position": [x[0] for x in self.spawnedPositions],
		                      "scale": [x[0] for x in self.spawnedScales], "height": self.spawnedHeights}
		self.spawnedDictsB = {"time": self.spawnedTimes, "position": [x[1] for x in self.spawnedPositions],
		                      "scale": [x[1] for x in self.spawnedScales],"height": self.spawnedHeights}

		spawnedMetadynHillsDicts = [ self.spawnedDictsA, self.spawnedDictsB ] 

		parsedFile = types.SimpleNamespace(metadyn_hills=spawnedMetadynHillsDicts)
		outData = [ [types.SimpleNamespace(parsedFile=parsedFile)] ]

		self.stdOutObj = types.SimpleNamespace(data=outData)

	def _loadExpectedCaseA_fromCreator(self):
		expTimes = [0, 0, 20]
		expPositions = self.initPositions + self.spawnedPositions
		expHeights = [[h,h] for h in self.initHeights + self.spawnedHeights] #2 values since we have 2 dimensions
		expScales = self.initScales + self.spawnedScales

		currKwargs = {"times": expTimes, "positions": expPositions,
		              "scales": expScales, "heights": expHeights}
		expOutObj = aMetadynHillsHelp.MetadynHillsInfo(**currKwargs)
		return expOutObj


	def _runTestFunct(self):
		args = [self.stdOutObj]
		kwargs = {"creator": self.creatorObjA, "startHillsInfo": self.initHillInfoDict}
		return tCode.getMetadynHillsInfoFromStdOutObj(*args, **kwargs)

	def testRaisesIfBothCreatorAndStartHillsPassed(self):
		with self.assertRaises(ValueError):
			self._runTestFunct()

	def testExpectedCaseA_fromCreator(self):
		self.initHillInfoDict = None
		expObj = self._loadExpectedCaseA_fromCreator()
		actObj = self._runTestFunct()
		self.assertEqual(expObj,actObj)

	def testExpected_fromInpHillsInfo(self):
		self.initTimes = [0,0]
		self.createTestObjs()
		self.creatorObjA = None
		expObj = self._loadExpectedCaseA_fromCreator()
		actObj = self._runTestFunct()
		self.assertEqual(expObj,actObj)

	def testExpected_noInitiallySpawnedHills(self):
		#get expected
		expTimes = [20]
		expPositions, expScales = self.spawnedPositions, self.spawnedScales
		expHeights = [ [h,h] for h in self.spawnedHeights ]
		currKwargs = {"times":expTimes, "positions":expPositions, "scales":expScales, "heights":expHeights}
		expObj = aMetadynHillsHelp.MetadynHillsInfo(**currKwargs)

		#setup
		self.initHillInfoDict, self.creatorObjA = None, None

		#run and test
		actObj = self._runTestFunct()
		self.assertEqual(expObj, actObj)


	@unittest.skip("")
	def testRaisesWhenNoHillsInOutput(self):
		self.assertTrue(False)

	@unittest.skip("")
	def testDoesntRaiseWhenNoHillsInOutput(self):
		self.assertTrue(False)


	@unittest.skip("")
	def testExpected_onlyInitiallySpawnedHills(self):
		self.assertTrue(False)




class TestGetDictFromCreatorFuncts(unittest.TestCase):

	def setUp(self):
		self.addedMOs = 4
		self.absGridCutoff = 600
		self.relGridCutoff = 400
		self.charge = 0
		self.kPts = [12,12,1]
		self.basisObjs = [mock.Mock()]
		self.colVars = [mock.Mock(), mock.Mock()]
		self.createTestObjs()

	def createTestObjs(self):
		self.kwargDictA = {"absGridCutoff":self.absGridCutoff, "addedMOs":self.addedMOs,
		                   "basisObjs":self.basisObjs,"charge": self.charge,
		                   "colVars":self.colVars, "kPts":self.kPts} 
		self.creatorObjA = creatorHelp.CP2KCalcObjFactoryStandard(**self.kwargDictA)

	def testGetBasicSettingsFunct(self):
		expDict = {"absGridCutoff":self.absGridCutoff, "addedMOs":self.addedMOs,
		                   "charge": self.charge, "kPts":self.kPts} 
		actDict = tCode.getSelectedBasicInfoDictFromCreatorObj(self.creatorObjA)
		self.assertEqual(expDict,actDict)

	@mock.patch("gen_basis_helpers.cp2k.cp2k_creator_to_dict.getDictFromBasisObjList")
	def testGetBasisObjs(self, mockedDictFromBasisObjs):
		expObj = mock.Mock()
		expDict = {"basisObjs":expObj}
		mockedDictFromBasisObjs.side_effect = lambda *args: expObj
		actDict = tCode._getBasisObjsInfoFromCreatorObj(self.creatorObjA)
		mockedDictFromBasisObjs.assert_called_with(self.basisObjs)
		self.assertEqual(expDict,actDict)

	def testGetColvarsInfo_toDictPresent(self):
		expDictA, expDictB = {"colVarA_attr": "colVarA_val"}, {"colVarB_attr":"colVarB_val"}
		self.colVars[0].toDict.side_effect = lambda : expDictA
		self.colVars[1].toDict.side_effect = lambda : expDictB
		expDict = {"colVars": [expDictA, expDictB]}
		actDict = tCode._getColvarsInfoFromCreatorObj(self.creatorObjA)
		self.assertEqual(expDict, actDict)

	def testGetColVarsInfo_attrErrorForToDict(self):
		expDictA = {"colVarA_attr":"colVarA_val"}
		self.colVars[0].toDict.side_effect = lambda : expDictA
		def _raiseAttribError():
			raise AttributeError("")
		self.colVars[1].toDict.side_effect = _raiseAttribError
		expDict = {"colVars": [expDictA,None]}
		actDict = tCode._getColvarsInfoFromCreatorObj(self.creatorObjA)
		self.assertEqual(expDict, actDict)


class TestGetBasisDictFromBasisObjsList(unittest.TestCase):

	def setUp(self):
		self.eleA, self.eleB = "Mg", "O"
		self.basisA, self.basisB = "basis_set_a", "basis_set_b"
		self.potA, self.potB = "pot_a", "pot_b"
		self.basisFileA, self.basisFileB = "basis_file_a", "basis_file_b"
		self.ghostA, self.ghostB = False, True
		self.kindA, self.kindB = "Mg", "O_other" 
		self.createTestObjs()

	def createTestObjs(self):
		currKwargs = {"element":self.eleA, "basis":self.basisA, "potential":self.potA,
		              "basisFile":self.basisFileA, "potFile":self.potA, "ghost":self.ghostA,
		              "kind":self.kindA}
		self.testObjA = basisObjHelp.CP2KBasisObjStandard(**currKwargs)
		currKwargs = {"element":self.eleB, "basis":self.basisB, "potential":self.potB,
		              "basisFile":self.basisFileB, "potFile":self.potB, "ghost":self.ghostB,
		              "kind":self.kindB}
		self.testObjB = basisObjHelp.CP2KBasisObjStandard(**currKwargs)
		self.basisObjs = [self.testObjA, self.testObjB]

	def testToAndFromDictConsistentForCaseA(self):
		basisObjDict = tCode.getDictFromBasisObjList(self.basisObjs)
		expBasisObjsList = self.basisObjs
		actBasisObjsList = tCode.getBasisObjListFromDict(basisObjDict)

		#compare actual/expected
		sortedExpList = sorted(expBasisObjsList, key=lambda x:x.basis)
		sortedActList = sorted(actBasisObjsList, key=lambda x:x.basis)
		self.assertEqual(sortedExpList, sortedActList)




