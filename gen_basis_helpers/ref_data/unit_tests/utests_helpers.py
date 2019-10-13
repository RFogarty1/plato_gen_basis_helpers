#!/usr/bin/python3

import unittest
import unittest.mock as mock

import gen_basis_helpers.ref_data.helpers_ref_data as tCode

import plato_pylib.shared.ucell_class as UCell

class TestGetDimerSepFromUCell(unittest.TestCase):

	def setUp(self):
		self.edgeLength = 20
		self.atomSep = 5.0
		self.createUCellObj()

	def createUCellObj(self):
		lattVects = [ [self.edgeLength, 0.0, 0.0],
		              [0.0, self.edgeLength, 0.0],
		              [0.0, 0.0, self.edgeLength] ]

		cartCoords = [ [0.0, 0.0, 1.0, "Zr"],
		                [0.0, 0.0, 1.0+self.atomSep, "Zr"] ]

		self.uCellObj = UCell.UnitCell.fromLattVects(lattVects)
		self.uCellObj.cartCoords = cartCoords

	def testAssertErrorForNonDimerCalc(self):
		endCart = self.uCellObj.cartCoords
		endCart.pop()
		self.uCellObj.cartCoords = endCart

		with self.assertRaises(AssertionError):
			tCode.getDimerSepFromUCellObj(self.uCellObj)

	def testGetExpectedValInSimpleCase(self):
		expSep = self.atomSep
		actualSep = tCode.getDimerSepFromUCellObj(self.uCellObj)
		self.assertAlmostEqual(expSep,actualSep)


class TestGetDimerSepFromCastepOutFile(unittest.TestCase):

	def setUp(self):
		self.testInpFile = "fake_inp_file"

	@mock.patch("gen_basis_helpers.ref_data.helpers_ref_data.getDimerSepFromUCellObj")
	@mock.patch("gen_basis_helpers.ref_data.helpers_ref_data.parseCastep.parseCastepOutfile")
	def testExpectedCallsSimpleCase(self,parseMock, sepFromUCellMock):
		expSep = 4.0
		fakeUCell = mock.Mock()
		parseMockRetVal = {"unitCell":fakeUCell}
		sepFromUCellMock.return_value = expSep
		parseMock.return_value = parseMockRetVal

		actSep = tCode.getDimerSepFromCastepOutFile(self.testInpFile)
		parseMock.assert_called_once_with(self.testInpFile)
		fakeUCell.convAngToBohr.assert_called_once_with()
		sepFromUCellMock.assert_called_once_with(fakeUCell)
		self.assertAlmostEqual(expSep,actSep)



if __name__ == '__main__':
	unittest.main()

