



import unittest
import unittest.mock as mock

import numpy as np

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.misc.smeared_charges_analysis as tCode


class TestGetInteractionEnergy(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.cartCoords = [ [0,0,2,"X"],
		                    [0,0,9,"Y"],
		                    [0,0,6,"Z"] ]

		self.charges = [0.2,-0.5,0.7]
		self.exponentDict = {"X":2,"Y":4,"Z":1}
		self.indicesA = [0,1,2]
		self.indicesB = [0,1,2]
		self.distLenConv = 1/3
		self.eConv = 2
		self._createTestObjs()

	def _createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

	def _runTestFunct(self):
		args = [self.cellA, self.exponentDict, self.charges, self.indicesA, self.indicesB]
		kwargs = {"distLenConv":self.distLenConv,"eConv":self.eConv}
		return tCode.getCoulombInteractionEnergyStandard(*args,**kwargs)

	def testInteractionEnergyFullMatrix(self):
		abInt, acInt, bcInt = -0.089752956514025, 0.092015940407703, -0.277933876243776
		expVal = self.eConv*(abInt + acInt + bcInt)
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal)

	def testInteractionEnergyPartialMatrix(self):
		self.indicesA = [0,1]
		self.indicesB = [2]
		abInt, acInt, bcInt = -0.089752956514025, 0.092015940407703, -0.277933876243776
		expVal = self.eConv*(acInt + bcInt)
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal)

class TestGetCoulombInteractionMatrix(unittest.TestCase):

	def setUp(self):
		self.lattParams, self.lattAngles = [10,10,10], [90,90,90]
		self.cartCoords = [ [0,0,2,"X"],
		                    [0,0,9,"Y"],
		                    [0,0,6,"Z"] ]

		self.charges = [0.2,-0.5,0.7]
		self.exponents = [2, 4, 1]
		self.indicesA = [0,1,2]
		self.indicesB = [0,1,2]
		self.distLenConv = 1
		self._createTestObjs()

	def _createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.cartCoords

	def _runTestFunct(self):
		args = [self.cellA, self.exponents, self.charges, self.indicesA, self.indicesB]
		kwargs = {"distLenConv":self.distLenConv}
		return tCode.getCoulombInteractionMatrix(*args, **kwargs)

	def testExpectedCaseA_allIndices(self):
		abInt, acInt, bcInt = -0.033333301221433, 0.034999864913425, -0.116649423063794

		selfInt = np.nan
		expMatrix = np.array( [ [selfInt, abInt  , acInt  ],
		                        [abInt  , selfInt, bcInt  ],
		                        [acInt  , bcInt  , selfInt] ] )
		actMatrix = self._runTestFunct()
		self.assertTrue( np.allclose(expMatrix,actMatrix, equal_nan=True) )
 
	def testExpectedCaseB_someIndices(self):
		self.indicesA = [1]
		self.indicesB = [0,2]

		abInt, bcInt = -0.033333301221433, -0.116649423063794

		expMatrix = np.array( [[abInt, bcInt]] )
		actMatrix = self._runTestFunct()
		self.assertTrue( np.allclose(expMatrix,actMatrix) )

	def testExpectedCaseC_distConv(self):
		self.distLenConv = 1/3
		abInt, acInt, bcInt = -0.089752956514025, 0.092015940407703, -0.277933876243776

		selfInt = np.nan
		expMatrix = np.array( [ [selfInt, abInt  , acInt  ],
		                        [abInt  , selfInt, bcInt  ],
		                        [acInt  , bcInt  , selfInt] ] )
		actMatrix = self._runTestFunct()
		self.assertTrue( np.allclose(expMatrix,actMatrix, equal_nan=True) )


class TestGetEnergyForTwoSmearedCharges(unittest.TestCase):

	def setUp(self):
		self.alphaA = 4
		self.alphaB = 6
		self.dist = 3

	def _runTestFunct(self):
		args = [self.dist, self.alphaA, self.alphaB]
		return tCode.getIntegralTwoSphericalGaussianSmearedDensities(*args)

	def testExpectedCaseA(self):
		""" Basically equivalent to me coding the algebra twice so...... """
		expVal = 0.33333333331686 #Just figured out in excel
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal)

	def testExpectedCaseB_erfNearZero(self):
		self.dist = 0.5
		expVal = 1.4533566434154
		actVal = self._runTestFunct()
		self.assertAlmostEqual(expVal, actVal)

