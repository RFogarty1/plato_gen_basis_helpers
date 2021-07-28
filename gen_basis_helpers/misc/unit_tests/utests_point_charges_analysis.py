

import unittest
import unittest.mock as mock

import numpy as np

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.misc.point_charges_analysis as tCode


class TestEstimateCoulombEnergyBetweenIndicesInGeom(unittest.TestCase):

	def setUp(self):
		self.lenConv = 1

		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]

		self.charges = [0.5, -0.3, 0.4, 0.2]

		self.coordsA = [ [0,0,8,"X"],
		                 [0,0,1,"Y"],
		                 [0,0,5,"Y"],
		                 [0,0,3,"Y"] ]
		self.indicesA = [0,1,2,3]
		self.indicesB = [0,1,2,3]

		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams,lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coordsA

	def _runTestFunct(self):
		args = [self.cellA, self.charges, self.indicesA, self.indicesB]
		currKwargs = {"lenConv":self.lenConv}
		return tCode.getCoulombEnergyBetweenIndicesForPointCharges(*args, **currKwargs)

	def _runGetCoulombMatrix(self):
		args = [self.cellA, self.charges, self.indicesA, self.indicesB]
		currKwargs = {"lenConv":self.lenConv}
		return tCode.getCoulombEnergyInteractionMatrix(*args, **currKwargs)

	def testGetCoulombMatrixAllIndices(self):
		intAB, intAC, intAD = -0.71998, 0.959973333333333, 0.287992
		intBC, intBD, intCD = -0.431988, -0.431988, 0.575984
		
		expMatrix = [ [np.inf, intAB, intAC, intAD],
		              [intAB, np.inf, intBC, intBD],
		              [intAC, intBC, np.inf, intCD],
		              [intAD, intBD, intCD, np.inf] ]

		actMatrix = self._runGetCoulombMatrix()
		self.assertTrue( np.allclose( np.array(expMatrix), np.array(actMatrix) ) )

	def testGetCoulombMatrixSelectedIndices(self):
		self.indicesA = [1,2] #BC
		self.indicesB = [0,2,3] #ACD

		intBA, intBC, intBD = -0.71998, -0.431988, -0.431988
		intCA, intCC, intCD = 0.959973333333333, np.inf, 0.575984

		expMatrix = [ [intBA, intBC, intBD],
		              [intCA, intCC, intCD] ]

		actMatrix = self._runGetCoulombMatrix()

		self.assertTrue( np.allclose( np.array(expMatrix), np.array(actMatrix) ) )

	def testGetCoulombInteractionEnergyA(self):
		self.indicesA = [2,3]
		self.indicesB = [0]

		int20, int30 = 0.959973333333333, 0.287992
		expInt = int20 + int30
		actInt = self._runTestFunct()
		self.assertAlmostEqual(expInt, actInt)

	def testExpectedCaseA_lenConvThree(self):
		self.lenConv = 3

		self.indicesA = [2,3]
		self.indicesB = [0]

		int20, int30 = 0.959973333333333, 0.287992
		expInt = (int20 + int30)*self.lenConv
		actInt = self._runTestFunct()
		self.assertAlmostEqual(expInt, actInt)

	def testExpectedIntEnergyAll(self):
		""" Total Coulomb energy between all charges in the system """
		intAB, intAC, intAD = -0.71998, 0.959973333333333, 0.287992
		intBC, intBD, intCD = -0.431988, -0.431988, 0.575984
		expInt = intAB + intAC + intAD + intBC + intBD + intCD
		actInt = self._runTestFunct()
		self.assertAlmostEqual(expInt,actInt)

class TestCoulombEnergyFromDistsAndCharges(unittest.TestCase):

	def setUp(self):
		self.qA, self.qB = -0.2, 0.4
		self.dist = 3
		self.lenConv, self.energyConv = 1, 1

	def _runTestFunct(self):
		currArgs = [self.qA, self.qB, self.dist] 
		currKwargs = {"lenConv":self.lenConv, "energyConv":self.energyConv}
		return tCode.getCoulombEnergyTwoPointsStandard(*currArgs, **currKwargs)

	def testExpectedChargesAndDistsA(self):
		expVal = -0.383989333333333
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal)

	def testExpectedChargesAndDistsB(self):
		self.qA, self.qB = 0.3, 0.6
		self.dist = 4.5
		expVal = 0.575984
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal)

	def testExpectedVals_energyConv(self):
		self.energyConv = 3
		expVal = -1.151968
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal)

	def testExpectedVals_distConv(self):
		self.lenConv = 4
		expVal = -1.53595733333333
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal)





