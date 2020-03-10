
import copy
import collections
import os
import itertools as it
import types

import numpy as np

import unittest
import unittest.mock as mock

import gen_basis_helpers.workflows.elastic_workflows as tCode



class TestHcpElasticWorkflowFactory(unittest.TestCase):

	def setUp(self):
		self.baseGeom = mock.Mock()
		self.strainValues = [-2,-1,0,1,2]
		self.creator = mock.Mock()
		self.eType = None
		self.workFolder = "fake/folder"
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.HcpElasticWorkflowCreator(baseGeom=self.baseGeom, strainValues=self.strainValues,
		                                                creator=self.creator, workFolder=self.workFolder, eType=self.eType)

	@mock.patch("gen_basis_helpers.workflows.elastic_workflows.HcpElasticConstantsWorkflow")
	@mock.patch("gen_basis_helpers.workflows.elastic_workflows.getRequiredStrainObjsForStructType")
	def testCorrectCallForStrainMatrices(self, mockedStrainMatrixGetter, mockedHcpFlow):
		self.testObjA.create()
		mockedStrainMatrixGetter.assert_called_once_with("hcp")


	@mock.patch("gen_basis_helpers.workflows.elastic_workflows.StressStrainWorkflowCreator")
	def testBaseFactoryObject(self, mockedStressStrainFactory):
		fakeFactory = "fake_factory"
		mockedStressStrainFactory.side_effect = lambda *args,**kwargs: "fake_factory"
		actFactory = self.testObjA._stressStrainBaseFactory
		mockedStressStrainFactory.assert_called_once_with(baseGeom=self.baseGeom, creator=self.creator, strainValues=self.strainValues, eType=self.eType)
		self.assertEqual(fakeFactory,actFactory)


	@mock.patch("gen_basis_helpers.workflows.elastic_workflows.HcpElasticConstantsWorkflow")
	@mock.patch("gen_basis_helpers.workflows.elastic_workflows.HcpElasticWorkflowCreator._getUnitStrainMatrices")
	@mock.patch("gen_basis_helpers.workflows.elastic_workflows.HcpElasticWorkflowCreator._stressStrainBaseFactory",new_callable=mock.PropertyMock)
	def testExpectedCallsToStressStrainFactory(self, mockedBaseFactory, mockedUStrain, mockedHcpFlow):
		baseFactory = mock.Mock()
		fakeStrains = [mock.Mock(), mock.Mock()]
		mockedBaseFactory.return_value = baseFactory
		mockedUStrain.side_effect = lambda: fakeStrains

		expWorkFolders = [os.path.join(self.workFolder, "strain_{}".format(x)) for x in range(len(fakeStrains))]
		expStrains = fakeStrains
		self.testObjA.create()

		for strain,workFolder in it.zip_longest(fakeStrains,expWorkFolders):
			baseFactory.create.assert_any_call(strain=strain, workFolder=workFolder)



class TestStressStrainWorkflowFactory(unittest.TestCase):

	def setUp(self):
		self.baseGeom = mock.Mock()
		self.strainValues = [-2,-1,0,1,2]
		self.strain = mock.Mock()
		self.creator = mock.Mock()
		self.eType = None
		self.workFolder = "fake/folder/path"
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.StressStrainWorkflowCreator(baseGeom=self.baseGeom, strainValues=self.strainValues,
		                                                  workFolder=self.workFolder, creator=self.creator, strain=self.strain, eType=self.eType)

	@mock.patch("gen_basis_helpers.workflows.elastic_workflows.StressStrainWorkflowCreator._getGeomList")
	def testCreatorCalledWithCorrectArgs(self, mockedGeomGetter):
		self.createTestObjs()
		expGeoms = self.strainValues
		expFolders = [self.workFolder for x in self.strainValues]
		expFileNames = ["strain_{:.3f}".format(x).replace(".","pt").replace("-","m") for x in self.strainValues]
		mockedGeomGetter.side_effect = lambda: expGeoms
		self.testObjA.create()
		for folderName, fileName, geom in it.zip_longest(expFolders, expFileNames, expGeoms):
			self.creator.create.assert_any_call(workFolder=folderName, fileName=fileName, geom=geom)

	@mock.patch("gen_basis_helpers.workflows.elastic_workflows.elasticHelp.getStrainedUnitCellStructsForUnitStrainVects")
	def testCorrectStrainedGeomsWithMock(self, mockedGeomStrainGetter):
		fakeGetterFunct = lambda uCell, sParams, uStrainMatices: [ sParams ]
		expReturnVal = self.strainValues
		mockedGeomStrainGetter.side_effect = fakeGetterFunct
		actReturnVal = self.testObjA._getGeomList()
		mockedGeomStrainGetter.assert_called_once_with(self.baseGeom, self.strainValues, [self.strain.strainMatrix])
		self.assertEqual(expReturnVal,actReturnVal)	

	@mock.patch("gen_basis_helpers.workflows.elastic_workflows.StressStrainWorkflow")
	@mock.patch("gen_basis_helpers.workflows.elastic_workflows.StressStrainWorkflowCreator._getGeomList")
	def testExpectedArgsPassedToStressStrainFlows(self, mockedGetGeoms, mockedFlow):
		expWorkflow = "fake_workflow"
		mockedFlow.side_effect = lambda *args, **kwargs: expWorkflow
		fakeGeoms = [2*x for x in self.strainValues]
		self.creator.create.side_effect = lambda geom=None,**kwargs: geom #Bounce back geom so we can test we pass the calcObjs it creates
		mockedGetGeoms.side_effect = lambda: fakeGeoms
		expCalcObjs = fakeGeoms
		actWorkflow = self.testObjA.create()
		mockedFlow.assert_called_once_with(expCalcObjs, self.strainValues, self.strain, eType=self.eType)
		self.assertEqual(expWorkflow, actWorkflow)	



class TestHcpElasticWorkflow(unittest.TestCase):

	def setUp(self):
		self.strainE3   = tCode.CrystalStrain([0,0,1,0,0,0])
		self.strainE12  = tCode.CrystalStrain([1,1,0,0,0,0])
		self.strainE123 = tCode.CrystalStrain([1,1,1,0,0,0])
		self.strainE45  = tCode.CrystalStrain([0,0,0,2,2,0])
		self.strainE6   = tCode.CrystalStrain([0,0,0,0,0,2])
		self.secondDerivVals = [1,2,4,4,5]

		self.allStrainsA = [ self.strainE3, self.strainE12, self.strainE123, self.strainE45, self.strainE6 ]
		self.allRunCommsListA = [ list(), ["runCommB","runCommBB"], list(), ["runCommC"], ["runCommD"] ]
		self.createTestObjs()

	def createTestObjs(self):
		self.mockFlowsA = [mock.Mock() for x in self.allStrainsA]

		for wflow,strain,runComms,secondDeriv in it.zip_longest(self.mockFlowsA,self.allStrainsA,self.allRunCommsListA,self.secondDerivVals):
			wflow.strain = strain
			wflow.output = [types.SimpleNamespace(secondDeriv=secondDeriv)]
			wflow.preRunShellComms = runComms
		self.testObjA = tCode.HcpElasticConstantsWorkflow( self.mockFlowsA )


	def testInitRaisesIfInputListTooLong(self):
		self.allStrainsA.append(None)
		with self.assertRaises(ValueError):
			self.createTestObjs()

	def testInitRaisesIfOneStrainIsWrong(self):
		self.allStrainsA[2] = tCode.CrystalStrain([1 for x in range(6)])
		with self.assertRaises(ValueError):
			self.createTestObjs()

	def testInitRaisesIfTwoStrainsAreTheSame(self):
		self.allStrainsA[2] = self.allStrainsA[3]
		with self.assertRaises(ValueError):
			self.createTestObjs()

	def testInitAlwaysLeavesStressStrainInCorrectOrder(self):
		orderA = copy.deepcopy( self.allStrainsA )
		origObj = copy.deepcopy( self.testObjA )

		self.allStrainsA = [ self.strainE123, self.strainE12, self.strainE45, self.strainE3, self.strainE6 ]
		self.assertNotEqual( orderA, self.allStrainsA )
		self.createTestObjs() 

		oldStrainsOnObj = [x.strain for x in origObj.stressStrainFlows]
		newStrainsOnObj = [x.strain for x in self.testObjA.stressStrainFlows]
		self.assertEqual(newStrainsOnObj, oldStrainsOnObj)

	def testExpPreRunCommsGenerated(self):
		expPreRunComms = list()
		for x in self.allRunCommsListA:
			expPreRunComms.extend(x)
		actRunComms = self.testObjA.preRunShellComms
		self.assertEqual(expPreRunComms, actRunComms)

	def testExpectedElasticConstantsGenerated(self):
		expConstants = collections.OrderedDict([ ["11",1.75], ["12",-0.75], ["13",0.25],
		                                         ["33",1], ["44",0.5] ])
		self.testObjA.run()
		actConstants = self.testObjA.output[0].elasticConsts


		for key in expConstants.keys():
			self.assertAlmostEqual( expConstants[key], actConstants[key] )


class TestStressStrainWorkflow(unittest.TestCase):

	def setUp(self):
		self.strainValsA = [-2,-1,0,1,2]
		self.calcEnergiesA = [x for x in range(len(self.strainValsA))]
		self.volumesA = [x for x in range(1,len(self.strainValsA)+1)]
		self.calcObjsA = [mock.Mock() for x in self.strainValsA]
		self.strainObjA = mock.Mock()
		self.eTypeA = "electronicTotalE" 
		self.createTestObjects()

	def createTestObjects(self):
		self.runComms = ["runComm{}".format(x) for x in range(len(self.strainValsA))]
		for obj,energy,volume in it.zip_longest(self.calcObjsA, self.calcEnergiesA,self.volumesA):
			obj.parsedFile.energies.electronicTotalE = energy #Linked to value of self.eTypeA
			obj.parsedFile.unitCell.volume = volume

		self.testObjA = tCode.StressStrainWorkflow(self.calcObjsA, self.strainValsA, self.strainObjA, eType=self.eTypeA)

	def testWriteFilesCalled(self):
		for x in self.calcObjsA:
			x.writeFile.assert_called_once_with()

	def testRunGivesExpectedStrainVsEnergy(self):
		expectedVals = [ [x,y] for x,y in it.zip_longest(self.strainValsA,self.calcEnergiesA) ]
		self.testObjA.run()
		actualVals = self.testObjA.output[0].strainVsEnergy
		for expVals, actVals in it.zip_longest(expectedVals, actualVals):
			self.assertAlmostEqual(expVals[0],actVals[0])
			self.assertAlmostEqual(expVals[1],actVals[1])
			self.assertTrue( len(expVals)==2 ), self.assertTrue( len(actVals)==2 )

	@mock.patch("gen_basis_helpers.workflows.elastic_workflows.StressStrainWorkflow._applyInpUnitsToGPaConversionFactor")
	def testRunGivesExpectedStrainVsStress(self, applyConvFactorMock):
		applyConvFactorMock.side_effect = lambda x: x
		self.testObjA.run()
		energiesPerVolume = [e/v for e,v in it.zip_longest(self.calcEnergiesA,self.volumesA)]
		expectedVals = [ [x,y] for x,y in it.zip_longest(self.strainValsA,energiesPerVolume) ]
		actualVals = self.testObjA.output[0].actVals
		for expVals, actVals in it.zip_longest(expectedVals, actualVals):
			self.assertAlmostEqual(expVals[0],actVals[0])
			self.assertAlmostEqual(expVals[1],actVals[1])
			self.assertTrue( len(expVals)==2 ), self.assertTrue( len(actVals)==2 )
		self.assertTrue(applyConvFactorMock.called) 

	@mock.patch("gen_basis_helpers.workflows.elastic_workflows.elasticHelp.polyFitAndGetSecondDeriv")
	def testCorrectArgsPassedToGetFitFunction(self, mockedFitter):
		expFitFunct = "fake_function"
		expSecondDeriv = "fake_second_deriv"

		mockedFitter.side_effect = lambda *args: types.SimpleNamespace(getFittedValuesForXVals=expFitFunct,
		                                                               secondDeriv=expSecondDeriv)
		self.testObjA.run()
		mockedFitter.assert_called_once_with( self.testObjA.output[0].actVals )
		actFitFunct = self.testObjA.output[0].fitFunct
		actSecondDeriv = self.testObjA.output[0].secondDeriv
		self.assertEqual(expFitFunct, actFitFunct)
		self.assertEqual(expSecondDeriv, actSecondDeriv)


class TestStrainObject(unittest.TestCase):

	def setUp(self):
		self.strainValsA = [1,2,3,4,5,6]
		self.strainValsB = [2,2,2,3,3,3]
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.CrystalStrain(self.strainValsA)
		self.testObjB = tCode.CrystalStrain(self.strainValsB)

	@mock.patch("gen_basis_helpers.workflows.elastic_workflows._getUnitStrainMatrix")
	def testExpectedStrainMatrixWithMocking(self, mockedGetter):
		mockedGetter.side_effect = lambda x: x if (x>0 and x<7) else None
		expStrainMatrix = sum([x*x for x in self.strainValsA]) #In reality its the sum of some matrices
		actStrainMatrix = self.testObjA.strainMatrix
		self.assertEqual(expStrainMatrix, actStrainMatrix)

	def testExpectedStrainMatrixWithoutMocking(self):
		expMatrix = np.array( [ [1,3,2.5], [3,2,2], [2.5,2,3] ] )
		actMatrix = self.testObjA.strainMatrix
		self.assertTrue(np.allclose(expMatrix,actMatrix))

	def testEquality_equalObjsCompareEqual(self):
		self.strainValsA = [x*0.1 for x in self.strainValsA] #Want non-integers
		self.strainValsB = list(self.strainValsA)
		self.createTestObjs()
		self.assertEqual( self.testObjA, self.testObjB )

	def testEquality_unequalObjsCompareUnequal(self):
		self.assertNotEqual( self.testObjA, self.testObjB )

	def testToStrAsExpectedWhenAllSet(self):
		expStr = "+".join(["{}eps{}".format(x,idx) for idx,x in enumerate(self.strainValsB,1)])
		actStr = self.testObjB.toStr()
		self.assertEqual(expStr,actStr)

	def testToStrAsExpectedWhenTwoSet(self):
		self.strainValsA = [0,1,0,0,2,0]
		self.createTestObjs()
		expStr = "1eps2+2eps5"
		actStr = self.testObjA.toStr()
		self.assertEqual(expStr,actStr)





