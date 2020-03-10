
import os
import types

import unittest
import unittest.mock as mock

import gen_basis_helpers.job_helpers.surface_energies as tCode
import gen_basis_helpers.shared.label_objs as labelHelp
import gen_basis_helpers.shared.calc_runners as calcRunners



class TestSurfaceEnergiesCreatorTemplate(unittest.TestCase):

	def setUp(self):
		self.eleKey, self.methodKey, self.structKey = "eKey","mKey","sKey"
		self.baseGeom = mock.Mock()
		self.cellDims = [1,1,1]
		self.kPoints = [12,12,12]
		self.expSurfKPoints = [12,12,1]
		self.surfType = "hcp0001"
		self.lenVac = 10
		self.nLayers = 4
		self.applyNLayersToBulk = None
		self.baseWorkfolder = os.path.join("fake","folder")
		self.createTestObjs()

	def createTestObjs(self):
		self.expSurfFileName = "surface_n{}_vac_{:.2f}".format(self.nLayers,self.lenVac).replace(".","pt")
		self.testObjA = tCode.CodeSpecificStandardInputCreatorTemplate(baseGeom=self.baseGeom, cellDims=self.cellDims,
		                                                               surfType=self.surfType, nLayers=self.nLayers,
		                                                               lenVac=self.lenVac, kPts=self.kPoints, eleKey=self.eleKey,
		                                                               methodKey=self.methodKey, structKey=self.structKey,
		                                                               applyNLayersToBulk=self.applyNLayersToBulk)

	@mock.patch("gen_basis_helpers.job_helpers.surface_energies.CodeSpecificStandardInputCreatorTemplate._getSurfaceObjClass")
	@mock.patch("gen_basis_helpers.job_helpers.surface_energies.supCell")
	def testBulkCellCreationNotApplyingNLayers(self, supCellMock, surfObjClassMock):
		bulkCell = self.testObjA._bulkCell
		supCellMock.superCellFromUCell.assert_called_once_with(self.baseGeom,self.cellDims)
		surfObjClassMock.assert_not_called() #If we're using the proper bulk cell we dont need to call the surface object

	@mock.patch("gen_basis_helpers.job_helpers.surface_energies.CodeSpecificStandardInputCreatorTemplate._getSurfaceObjClass")
	@mock.patch("gen_basis_helpers.job_helpers.surface_energies.supCell")
	def testApplyNLayersToBulk(self, supCellMock, surfObjClassMock):
		self.assertTrue( self.testObjA.applyNLayersToBulk is False ) #Check default argument set correctly
		self.applyNLayersToBulk = True
		surfObjClass = mock.Mock()
		surfObjClassMock.side_effect = lambda *args,**kwargs: surfObjClass

		self.createTestObjs()
		bulkCell = self.testObjA._bulkCell
		callArgs, callKwargs = surfObjClass.call_args
		self.assertEqual( callArgs[1:], (self.nLayers,0) ) #We want 0 vaccum


	def testRaisesForIncorrectSurfaceType(self):
		self.surfType = "fake_surfaceType"
		self.createTestObjs()
		with self.assertRaises(ValueError):
			self.testObjA._getSurfaceObjClass()

	@mock.patch("gen_basis_helpers.job_helpers.surface_energies.CodeSpecificStandardInputCreatorTemplate._getSurfaceObjClass")
	@mock.patch("gen_basis_helpers.job_helpers.surface_energies.CodeSpecificStandardInputCreatorTemplate._bulkCellNoSurfaceLayers",new_callable=mock.PropertyMock)
	def testCreatesSurfaceUnitCell(self, mockedBulkGetter, mockedSurfClassGetter):
		mockSurfClass = mock.Mock()
		expBaseCell = mock.Mock()
		mockedBulkGetter.return_value = expBaseCell
		mockedSurfClassGetter.side_effect = lambda *args,**kwargs : mockSurfClass
		outCell = self.testObjA._surfaceCell
		mockSurfClass.assert_called_once_with(expBaseCell, self.nLayers, self.lenVac)

	def testSurfaceFileName(self):
		actSurfFileName = self.testObjA._surfFileName
		self.assertEqual(self.expSurfFileName,actSurfFileName)

	@mock.patch("gen_basis_helpers.job_helpers.surface_energies.CodeSpecificStandardInputCreatorTemplate._surfaceCell",new_callable=mock.PropertyMock)
	@mock.patch("gen_basis_helpers.job_helpers.surface_energies.CodeSpecificStandardInputCreatorTemplate._createCalcObjCreator")
	def testGetSurfaceCreator(self, mockedObjCreator,mockedSurfGeom):
		expSurfGeom = mock.Mock()
		mockBaseCreatorObj = mock.Mock()
		mockedSurfGeom.return_value = expSurfGeom
		mockedObjCreator.side_effect = lambda *args,**kwargs: mockBaseCreatorObj
		self.testObjA._getSurfaceCalcObj()

		self.assertEqual(expSurfGeom, mockBaseCreatorObj.geom)
		self.assertEqual(self.expSurfKPoints, mockBaseCreatorObj.kPts)
		self.assertEqual(self.expSurfFileName, mockBaseCreatorObj.fileName)

	def testLabelProp(self):
		expLabel = labelHelp.StandardLabel(eleKey=self.eleKey, structKey=self.structKey, methodKey=self.methodKey)
		actLabel = self.testObjA.label
		self.assertEqual(expLabel,actLabel)




class TestMapFunction(unittest.TestCase):

	def setUp(self):
		self.surfEnergy = 20
		self.ePerAtomSurf = 12
		self.ePerAtomBulk = 13
		self.methodStr = "methA"
		self.eleKey = "ele"
		self.structKey = "hcp0001"
		self.createTestObjs()

	def createTestObjs(self):
		self.testFunctA = tCode.MapSurfaceEnergiesToStandardFormat()
		self.testWorkflowA = mock.Mock()
		self.testWorkflowA.output = [types.SimpleNamespace(surfaceEnergy=self.surfEnergy,surfEPerAtom=self.ePerAtomSurf,
		                                                   bulkEPerAtom=self.ePerAtomBulk)]
		testLabelA = labelHelp.StandardLabel( eleKey=self.eleKey, methodKey=self.methodStr, structKey=self.structKey )
		self.standardInpObjA = calcRunners.StandardInputObj(self.testWorkflowA, testLabelA) #Dont need to set map funct

	def runTestFunct(self):
		return self.testFunctA(self.standardInpObjA)

	def testRunMethodCalled(self):
		self.runTestFunct()
		self.testWorkflowA.run.assert_called_once_with()

	def testTableDataOutput(self):
		self.runTestFunct()
		expTableData = [self.methodStr, "{:.4f}".format(self.surfEnergy)]
		actOutput = self.testFunctA(self.standardInpObjA)
		actTableData = actOutput.tableData
		self.assertEqual(expTableData,actTableData)

	def testTableWithEPerAtomOutput(self):
		self.runTestFunct()
		expTableData = [self.methodStr, "{:.4f}".format(self.surfEnergy),
		                "{:.3g}".format(self.ePerAtomBulk), "{:.3g}".format(self.ePerAtomSurf)]
		actOutput = self.testFunctA(self.standardInpObjA)
		actTableData = actOutput.tableWithEPerAtomVals
		self.assertEqual(expTableData, actTableData)



