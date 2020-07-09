
import itertools as it
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp
import types

import gen_basis_helpers.shared.method_objs as methObjs
import gen_basis_helpers.workflows.stacking_fault_workflows as tCode


class TestStackingFaultWorkflow(unittest.TestCase):

	def setUp(self):
		self.displacements = [0,0.2,0.4,0.6,0.8,1.0]
		self.energies = [x+2 for x in range(len(self.displacements))]
		self.perfectEnergy = 1
		self.surfaceAreaAll = 0.5
		self.fitterObj = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		#Create the stub parsed files
		stubParsedFiles = list()
		for energy in self.energies:
			stubEnergiesObj = types.SimpleNamespace(electronicTotalE=energy)
			currStub = types.SimpleNamespace(energies=stubEnergiesObj, unitCell=mock.Mock())
			stubParsedFiles.append(currStub)

		stubCalcObjs = [methObjs.StubCalcMethodFromParsedFileObject(x) for x in stubParsedFiles]

		#Create a single perfect structure calc obj
		perfectEnergiesObj = types.SimpleNamespace(electronicTotalE=self.perfectEnergy)
		perfectParsedFile = types.SimpleNamespace(energies=perfectEnergiesObj, unitCell=mock.Mock())
		stubPerfectCalcObj = methObjs.StubCalcMethodFromParsedFileObject(perfectParsedFile)

		self.testWorkflowA = tCode.StackingFaultWorkflow(stubPerfectCalcObj, stubCalcObjs, self.displacements, self.fitterObj)

	@mock.patch("gen_basis_helpers.workflows.stacking_fault_workflows._getABSurfaceAreaFromParsedFile")
	def testExpectedRawEnergies(self, mockedSurfArea):
		mockedSurfArea.side_effect = lambda *args: self.surfaceAreaAll
		expVals = self.energies
		self.testWorkflowA.run()
		actVals = self.testWorkflowA.output.totalEnergiesRaw
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expVals,actVals)]

	@mock.patch("gen_basis_helpers.workflows.stacking_fault_workflows._getABSurfaceAreaFromParsedFile")
	def testExpectedEnergiesPerSurfaceAreaRaw(self, mockedSurfArea):
		mockedSurfArea.side_effect = lambda *args: self.surfaceAreaAll
		expVals = [x/self.surfaceAreaAll for x in self.energies]
		self.testWorkflowA.run()
		actVals = self.testWorkflowA.output.energiesPerSurfaceAreaRaw
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expVals,actVals)]

	@mock.patch("gen_basis_helpers.workflows.stacking_fault_workflows._getABSurfaceAreaFromParsedFile")
	def testExpectedDisplacements(self, mockedSurfArea):
		mockedSurfArea.side_effect = lambda *args: self.surfaceAreaAll
		self.testWorkflowA.run()
		expDisplacements = self.displacements
		actDisplacements = self.testWorkflowA.output.displacements
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expDisplacements,actDisplacements)]

	@mock.patch("gen_basis_helpers.workflows.stacking_fault_workflows._getABSurfaceAreaFromParsedFile")
	def testExpectedPerfectObjEnergy(self,mockedSurfArea):
		self.testWorkflowA.run()
		expEnergy = self.perfectEnergy
		actEnergy = self.testWorkflowA.output.totalEForPerfectStruct
		self.assertAlmostEqual(expEnergy,actEnergy)

	@mock.patch("gen_basis_helpers.workflows.stacking_fault_workflows._getABSurfaceAreaFromParsedFile")	
	def testExpectedPerfectObjEnergyPerArea(self,mockedSurfArea):
		mockedSurfArea.side_effect = lambda *args: self.surfaceAreaAll
		self.testWorkflowA.run()
		expEnergy = self.perfectEnergy / self.surfaceAreaAll
		actEnergy = self.testWorkflowA.output.ePerAreaForPerfectStruct
		self.assertAlmostEqual(expEnergy, actEnergy)

	@mock.patch("gen_basis_helpers.workflows.stacking_fault_workflows._getABSurfaceAreaFromParsedFile")	
	def testExpectedStackingFaultEnergiesForEachStruct(self, mockedSurfArea):
		mockedSurfArea.side_effect = lambda *args: self.surfaceAreaAll
		self.testWorkflowA.run()
		perfectVal = self.perfectEnergy / self.surfaceAreaAll #Define the zero
		ePerSurfArea = [e/self.surfaceAreaAll for e in self.energies]
		expVals = [energy-perfectVal for energy in ePerSurfArea]
		actVals = self.testWorkflowA.output.stackFaultEnergiesRaw
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expVals,actVals)]

	@mock.patch("gen_basis_helpers.workflows.stacking_fault_workflows.StackingFaultWorkflow._getStackFaultEnergies")
	@mock.patch("gen_basis_helpers.workflows.stacking_fault_workflows._getABSurfaceAreaFromParsedFile")	
	def testExpectedCallMadeToFitterObj(self, mockedSurfArea, mockStackFaultGetter):
		expStackFaultEnergies = mock.Mock()
		mockStackFaultGetter.side_effect = lambda *args:expStackFaultEnergies
		self.testWorkflowA.run()
		expCallVals = [self.displacements, expStackFaultEnergies]
		actCallVals,unused = self.fitterObj.fitStackingFaultEnergies.call_args
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expCallVals[0],actCallVals[0])]
		self.assertEqual(expCallVals[1],actCallVals[1])

	@mock.patch("gen_basis_helpers.workflows.stacking_fault_workflows._getABSurfaceAreaFromParsedFile")	
	def testExpectedFitResultsInOutput(self,mockedSurfArea):
		expFitResults = mock.Mock()
		self.fitterObj.fitStackingFaultEnergies.side_effect = lambda *args: expFitResults
		self.testWorkflowA.run()
		actFitResults = self.testWorkflowA.output.fitResult
		self.assertEqual(expFitResults,actFitResults)


class TestGetABSurfaceAreaFromParsedFile(unittest.TestCase):

	def setUp(self):
		self.lattParams = [2,2,2]
		self.lattAngles = [90,90,90]
		self.fractPositions = [ [0.0,0.0,0.0] ]
		self.eleList = ["X" for x in self.fractPositions]
		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"lattParams":self.lattParams, "lattAngles":self.lattAngles,
		             "fractCoords":self.fractPositions, "elementList":self.eleList}
		self.testCellA = uCellHelp.UnitCell(**kwargDict)
		self.parsedFileStubA = types.SimpleNamespace(unitCell=self.testCellA)

	def testOnCubicCellA(self):
		expSurfaceArea = self.lattParams[0]*self.lattParams[1]
		actSurfaceArea = tCode._getABSurfaceAreaFromParsedFile(self.parsedFileStubA)
		self.assertAlmostEqual(expSurfaceArea, actSurfaceArea)


class TestStackingFaultPolyFitter(unittest.TestCase):

	def setUp(self):
		self.dispVals = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
		fitFunctA = lambda x: (2*(x**2)) + (3*x) + 4
		self.stackFaultVals = [fitFunctA(x) for x in self.dispVals]
		self.polyOrder = 2
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.StackingFaultFitterPolyStandard( self.polyOrder )

	def _getActFitObj(self):
		return self.testObjA.fitStackingFaultEnergies(self.dispVals, self.stackFaultVals)

	def testQuadraticFitWorksForQuadraticFunction(self):
		actFitResult = self._getActFitObj()
		expFitParams = [2,3,4] #in this particular case, coefficients of a quadratic equation
		actFitParams = actFitResult.fitParams
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expFitParams, actFitParams)]

	@mock.patch("gen_basis_helpers.workflows.stacking_fault_workflows.StackingFaultFitterPolyStandard._getFitCoeffs")
	def testExpectedFitValsAtInputDisps(self, mockedCoeffGetter):
		mockedCoeffGetter.side_effect = lambda *args: [1,2,3] #corresponds to x**2 + 2x + 3
		expFitFunct = lambda x: (x**2) + (2*x) + 3
		expVals = [expFitFunct(x) for x in self.dispVals]
		actFitResult = self._getActFitObj()
		actVals = actFitResult.fitFaultValsAtInputDisps
		[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expVals, actVals)]

	@mock.patch("gen_basis_helpers.workflows.stacking_fault_workflows.StackingFaultFitterPolyStandard._getFitCoeffs")
	def testExpectedFaultEnergiesFromFit(self, mockedCoeffGetter):
		#First we set the polynomial such that the maximum is at an expected value
		aParam, maxPos = -2,0.7
		bParam = -1*2*aParam*maxPos

		mockedCoeffGetter.side_effect = lambda *args: [aParam,bParam,3] #corresponds to x**2 + 2x + 3
		expFitFunct = lambda x: (aParam*(x**2)) + (bParam*x) + 3
		expUnstableFaultEnergy = expFitFunct(maxPos)
		expIntrinsicFaultEnergy = expFitFunct(1)

		actFitRes = self._getActFitObj()
		actUnstableFaultEnergy = actFitRes.unstableFaultEnergy
		actIntrinsicFaultEnergy = actFitRes.intrinsicFaultEnergy

		self.assertAlmostEqual(expIntrinsicFaultEnergy, actIntrinsicFaultEnergy)
		self.assertAlmostEqual(expUnstableFaultEnergy, actUnstableFaultEnergy)

	@mock.patch("gen_basis_helpers.workflows.stacking_fault_workflows.StackingFaultFitterPolyStandard._getFitCoeffs")
	def testGetGoodnessOfFit(self,mockedCoeffGetter):
		mockedCoeffGetter.side_effect = lambda *args: [0,0,0.5]
		summedRelativeErrors = sum( [abs((x-0.5)/x) for x in self.stackFaultVals] )
		expError = summedRelativeErrors / len(self.stackFaultVals)
		actError = self._getActFitObj().goodnessOfFit
		self.assertAlmostEqual(expError,actError)

