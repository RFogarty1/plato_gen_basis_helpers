

import unittest

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.atom_combo_core as atomComboCoreHelp
import gen_basis_helpers.analyse_md.atom_combo_populators as atomComboPopulatorsHelp

import gen_basis_helpers.analyse_md.water_derivative_species_classifiers as tCode


class TestWaterDerivativeClassifiers_ohDistsOnly(unittest.TestCase):

	def setUp(self):
		#Geometry
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		self.dudA = [ [0,0,0,"X"] ]
		self.water = [ [0,0,0,"O"], [0.5,0.5,0,"H"], [0.5,-0.5,0,"H"] ]
		self.freeOxy = [ [4,0,0,"O"] ]
		self.dudB = [ [0,0,0,"Y"] ]
		self.hydroxyl = [ [0,7,0,"O"], [0,8,0,"H"] ] #Only 2 from the first water oxygen; so tests min-dist really
		self.hydronium = [ [5,5,0,"O"], [5.5,5.5,0,"H"], [4.5,5,0,"H"], [5.5,4.5,0,"H"] ]
		self.cartCoords = self.dudA + self.water + self.freeOxy + self.dudB + self.hydroxyl + self.hydronium

		#Will change in each test
		self.oxyIndices = [1,4,6,8]
		self.hyIndices = [2,3,7,9,10,11]
		self.maxOHDist = 2
		self.nNebs = 2

		self.createTestObjs()

	def createTestObjs(self):
		#Create geometry
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

		#Create + populate the sparseMatrix object; we only need a distmatrix populator for this
		allIndices = [x for x in reversed(range(len(self.cartCoords)))] #Reverse so i dont assume there in order implicitly
		distMatrixPopulator = atomComboPopulatorsHelp._DistMatrixPopulator(allIndices, allIndices)
		self.sparseMatrixCalculator = atomComboCoreHelp._SparseMatrixCalculatorStandard([distMatrixPopulator])
		self.sparseMatrixCalculator.calcMatricesForGeom(self.cellA)

		#Create classifier object
		currArgs = [self.oxyIndices, self.hyIndices]
		self.classifierObj = tCode._WaterDerivativeDistanceOnlyClassifierGeneric(*currArgs, maxOHDist=self.maxOHDist, nNebs=self.nNebs)


	def _runTestFunct(self):
		return self.classifierObj.classify(self.sparseMatrixCalculator)

	def testWaterIdentified(self):
		self.nNebs = 2
		self.createTestObjs()
		expVals = [ [ [1] ], [ [2,3] ] ]
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)

	def testHydroxylIdentified(self):
		self.nNebs = 1
		self.createTestObjs()
		expVals = [ [[6]], [[7]] ]
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)

	def testFreeOxygenIdentified(self):
		self.nNebs = 0
		self.createTestObjs()
		expVals = [ [[4]], [list()] ]
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)

	def testHydroniumIdentified(self):
		self.nNebs = 3
		self.createTestObjs()
		expVals = [ [ [8] ], [ [9,10,11] ] ]
		actVals = self._runTestFunct()
		self.assertEqual(expVals, actVals)



