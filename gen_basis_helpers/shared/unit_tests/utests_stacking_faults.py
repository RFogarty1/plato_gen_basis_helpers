


import copy
import itertools as it
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCell

import gen_basis_helpers.shared.geom_constraints as geomConstrHelp
import gen_basis_helpers.shared.simple_vector_maths as vectHelp
import gen_basis_helpers.shared.stacking_fault_geoms as tCode


class TestHcpI2StackingFaultGeomGenerator(unittest.TestCase):

	def setUp(self):
		self.lattParams = [2,2,3]
		self.lattAngles = [90,90,120]
		#Default to a 1x1x3 hcp Cell
		self.fractPositions = [[0.0, 0.0, 0.0],
		                       [1/3, 2/3, 1/6],
		                       [0.0, 0.0, 2/6],
		                       [0.0, 0.0, 4/6],
		                       [1/3, 2/3, 3/6],
		                       [1/3, 2/3, 5/6]]

		self.dispFactor = 1.0
		self.centralAtomIdx = 3
		self.createTestObjs()

	def createTestObjs(self):
		atomList = ["X" for x in self.fractPositions]
		self.testCellA = uCell.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles,
		                                fractCoords=self.fractPositions, elementList=atomList)
		self.testObjA = tCode.HcpI2StackingFaultGeomGenerator(centralIdx=self.centralAtomIdx)

	def _runTestFunct(self):
		return self.testObjA.getGeomForGivenDisplacement(self.testCellA, self.dispFactor, centralIdx=self.centralAtomIdx)

	def testExpectedGeomForFullDisp_testA(self):
		expCartCoords = self._getExpCartCoordsForFullDispA()
		expCell = copy.deepcopy(self.testCellA)
		expCell.cartCoords = expCartCoords
		actCell = self._runTestFunct()
		self.assertEqual(expCell,actCell)

	def testUsesSensibleCentralAtomIdx_1x1x3Cell(self):
		""" Make sure we use the near-central surface planes for the dislocation by default """
		self.centralAtomIdx = None
		self.createTestObjs()
		expIdx = 4 #second A in ABABAB, listing high->low
		actIdx = self.testObjA._getCentralAtomIdx(self.testCellA)
		self.assertEqual(expIdx,actIdx)

	def testUsesSensibleCentralAtomIdx_1x1x4Cell(self):
		""" Even number of layers (1 layer=2atoms) makes it easier to acidentally split along an A-B plane """
		self.fractPositions = [[0.0, 0.0, 0.0],
		                       [1/3, 2/3, 1/8],
		                       [0.0, 0.0, 2/8],
		                       [1/3, 2/3, 3/8],
		                       [0.0, 0.0, 4/8],
		                       [1/3, 2/3, 5/8],
		                       [0.0, 0.0, 6/8],
		                       [1/3, 2/3, 7/8]]
		self.centralAtomIdx = None
		self.createTestObjs()
		expIdx = 3 #Third A in ABABABAB assuming high to low
		actIdx = self.testObjA._getCentralAtomIdx(self.testCellA)
		self.assertEqual(expIdx,actIdx)

	#Algorithm for finding the central index works on the assumption that each plane has the same number of atoms in it
	def testRaisesErrorGettingDefaultCentralIdxForDifferentNumbersOfAtomsInDiffPlanes(self):
		self.fractPositions.append(self.fractPositions[-1])
		self.centralAtomIdx = None
		self.createTestObjs()
		with self.assertRaises(AssertionError):
			self._runTestFunct()

	def testRaisesIfInputCellAnglesIncorrect(self):
		self.lattAngles = [90,90,90]
		self.createTestObjs()
		dispVal = 1.0
		with self.assertRaises(ValueError):
			self.testObjA.getGeomForGivenDisplacement(self.testCellA,dispVal)

	def testExpectedGeomConstriantsReturned(self):
		expAtomicConstraints = geomConstrHelp.AtomicPositionConstraints(atomicCartConstraints=self._getAllExpectedAtomicGeomConstraints())
		expCellConstraints = geomConstrHelp.CellConstraints([True,True,True],[True,True,True])
		expConstraintObj = geomConstrHelp.GeomConstraints(expAtomicConstraints,expCellConstraints)
		actConstraintObj = self.testObjA.getGeomConstraints(self.testCellA)
		self.assertEqual(expConstraintObj,actConstraintObj)

	def _getAllExpectedAtomicGeomConstraints(self):
		outConstraints = list()
		for idx,unused in enumerate(self.fractPositions):
			currConstraint = geomConstrHelp.AtomicCartesianConstraint(idx, fixX=True,fixY=True)
			outConstraints.append(currConstraint)
		return outConstraints

	#Relies on specific properties of the default lattParams/fractPositions
	def _getExpCartCoordsForFullDispA(self):
		expCartCoords = copy.deepcopy(self.testCellA.cartCoords)

		#Figure out the displacement vector
		aVect, bVect, unused = self.testCellA.lattVects
		uVectA, uVectB = [vectHelp.getUnitVectorFromInpVector(x) for x in [aVect,bVect]]
		unitVectorDisplacement = uVectA
		dispOneDisplacementVector = self.lattParams[0] * (1/3) 
		dispVector = [x*dispOneDisplacementVector for x in unitVectorDisplacement]

		#Apply the displacement to the bottom half of the crystal
		maxZVal = expCartCoords[self.centralAtomIdx][2]
		zTol = 1e-2
		for currCoords in expCartCoords:
			if (currCoords[-2] <= maxZVal+zTol):
				currCoords[:3] = [x+d for x,d in it.zip_longest(currCoords[:3],dispVector)]

		return expCartCoords
