
import os
import unittest
import unittest.mock as mock

import gen_basis_helpers.job_helpers.plato.elastic_help as tCode

class TestHcpInputCreator(unittest.TestCase):

	def setUp(self):
		self.eleStr = "Mg"
		self.methodStr = "dft2_some_method" #The methodStr actually tells me what grid to take so....
		self.structStr = "some_struct"
		self.convDatabase = mock.Mock()
		self.eleDatabase = mock.Mock()
		self.strainValues = [-1,0,1,2]
		self.baseWorkFolder = "fake/folder/path"
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = { "eleKey":self.eleStr, "methodKey":self.methodStr, "structKey":self.structStr,
		              "convDatabase":self.convDatabase, "eleDatabase":self.eleDatabase,
		              "strainValues":self.strainValues, "baseWorkFolder":self.baseWorkFolder }
		self.testObjA = tCode.PlatoHcpElasticStandardInputCreator(**kwargDict)

	@mock.patch("gen_basis_helpers.job_helpers.plato.elastic_help.platoCreator")
	@mock.patch("gen_basis_helpers.job_helpers.plato.elastic_help.elasHelp.HcpElasticStandardInputCreator")
	def testOutFolderNamesCorrectlyPassedToCreator(self, mockedCreator, mockedPlatoCreator):
		self.testObjA.create()
		expKwargs = {"baseWorkFolder": self.baseWorkFolder,
		             "extToWorkFolder": os.path.join("elastic",self.eleStr,"hcp",self.methodStr)}
		actArgs, actKwargs = mockedCreator.call_args
		for key in expKwargs:
			self.assertEqual( expKwargs[key], actKwargs[key] )

	@mock.patch("gen_basis_helpers.job_helpers.plato.elastic_help.platoCreator")
	@mock.patch("gen_basis_helpers.job_helpers.plato.elastic_help.elasHelp.HcpElasticStandardInputCreator")
	def testExpectedBaseGeomPassedOn(self, mockedCreator, mockedPlatoCreator):
		expGeom = "fake_geom"
		ourEleDatabase = getattr(self.eleDatabase, self.eleStr.capitalize())
		ourEleDatabase.getPlaneWaveGeom.side_effect = lambda key: expGeom
		self.testObjA.create()
		expKwargs = {"baseGeom" : expGeom}
		actArgs, actKwargs = mockedCreator.call_args
		ourEleDatabase.getPlaneWaveGeom.assert_called_once_with(self.structStr) #Could hard-code hcp tbh
		for key in expKwargs:
			self.assertEqual( expKwargs[key], actKwargs[key] )


	@mock.patch("gen_basis_helpers.job_helpers.plato.elastic_help.platoCreator.PlatoCalcObjFactoryStandard")
	@mock.patch("gen_basis_helpers.job_helpers.plato.elastic_help.elasHelp.HcpElasticStandardInputCreator")
	def testExpectedCreatorObj(self, mockedStandardInpCreator, mockedPlatoCreator):
		expKPts, expDataset, expGridVals, expMethodKey = "fake_k", "fake_dataset", "fake_grid_vals", self.methodStr
		ourEleStructDb, ourEleConvDb = [getattr(x,self.eleStr.capitalize()) for x in [self.eleDatabase, self.convDatabase]]
		ourEleStructDb.modelFiles.dft2PlatoPath = expDataset #any of the platoPath vars should work really
		ourEleConvDb.kptGridVals.getKptsPrimCell.side_effect = lambda key: expKPts
		ourEleConvDb.integGridVals.getPrimCellDft2AngularGrid.side_effect = lambda key: expGridVals

		self.testObjA.create()
		expKwargs = {"dataSet":expDataset, "gridVals":expGridVals, "kPts":expKPts, "methodStr":expMethodKey}
		actArgs,actKwargs = mockedPlatoCreator.call_args
		ourEleConvDb.kptGridVals.getKptsPrimCell.assert_called_once_with(self.structStr)
		ourEleConvDb.integGridVals.getPrimCellDft2AngularGrid.assert_called_once_with(self.structStr)
		for key in expKwargs:
			self.assertEqual( expKwargs[key], actKwargs[key] )


#TODO: Check expected structure and creator passed on (thats the main difficult part)

#	@mock.patch("gen_basis_helpers.job_helpers.plato.elastic_help.elasHelp.HcpElasticStandardInputCreator")
#	def test


