
import copy
import itertools as it
import types

import numpy as np

import unittest
import unittest.mock as mock

import gen_basis_helpers.workflows.elastic_workflows as tCode


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
		self.strainValsB = list(self.strainValsB)
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





