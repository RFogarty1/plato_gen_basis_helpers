

import unittest
import unittest.mock

import plato_pylib.shared.ucell_class as uCellHelp

import gen_basis_helpers.analyse_md.get_neb_lists as tCode

class TestGetUniqueNeighbourLists(unittest.TestCase):

	def setUp(self):
		self.cutoff = 6
		self.lattParams = [10,10,10]
		self.lattAngles = [90,90,90]
		self.coords = [ [5,5,5,"X"],
		                [5,5,6,"Y"],
		                [5,5,7,"Z"] ]
		self.createTestObjs()

	def createTestObjs(self):
		self.cellA = uCellHelp.UnitCell(lattParams=self.lattParams, lattAngles=self.lattAngles)
		self.cellA.cartCoords = self.coords

	def testExpected_allNebs_PBCsIrrelevant(self):
		expNebLists = [ [1,2], [0,2], [0,1] ]
		actNebLists = tCode.getNeighbourListsForInpCell_imagesMappedToCentral(self.cellA, self.cutoff)
		self.assertEqual(expNebLists, actNebLists)

	def testExpected_notAllNebs_PBCsIrrelevant(self):
		self.cutoff = 1.5
		expNebLists = [ [1], [0,2], [1] ]
		actNebLists = tCode.getNeighbourListsForInpCell_imagesMappedToCentral(self.cellA, self.cutoff)
		self.assertEqual(expNebLists, actNebLists)

	def testExpected_allNebs_withPBCs(self):
		self.coords[0][-2] = 9
		self.coords[1][-2] = 9
		self.coords[2][-2] = 1
		self.createTestObjs()

		expNebLists = [ [1,2], [0,2], [0,1] ]
		actNebLists = tCode.getNeighbourListsForInpCell_imagesMappedToCentral(self.cellA, self.cutoff)
		self.assertEqual(expNebLists,actNebLists)

