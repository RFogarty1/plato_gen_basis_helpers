
import itertools as it
import unittest
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp
import gen_basis_helpers.shared.surfaces as surfHelp
import gen_basis_helpers.adsorption.surf_to_sites_shared as tCode


class TestOrderSitesBasedOnABCentrality(unittest.TestCase):

	def setUp(self):
		self.lattParamsA = [2,2,2]
		self.lattAnglesA = [90,90,90]
		self.siteC = [1,1,1]
		self.siteB = [0.7,0.7,1]
		self.siteA = [0.3,0.3,1]
		self.createTestObjs()

	def createTestObjs(self):
		self.testSitesA = [self.siteA, self.siteB, self.siteC]
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParamsA, lattAngles=self.lattAnglesA,
		                                fractCoords=[[0,0,0]], elementList=["Mg"])
		self.surfA = surfHelp.GenericSurface(self.cellA, 1, lenVac=2)
		self.testFunctA = lambda inpSurf:self.testSitesA 

	def testExpectedForSimpleCaseA(self):
		testFunct = tCode.getWrappedSurfaceToSitesToReorderBasedOnABCentrality(self.testFunctA)
		expPositions = [self.siteC, self.siteB, self.siteA]
		actPositions = testFunct(self.surfA)
		for exp,act in it.zip_longest(expPositions,actPositions):
			[self.assertAlmostEqual(e,a) for e,a in it.zip_longest(exp,act)]



