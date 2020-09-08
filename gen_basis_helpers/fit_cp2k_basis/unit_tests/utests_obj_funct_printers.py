
import unittest
import unittest.mock as mock

import gen_basis_helpers.fit_cp2k_basis.obj_funct_printers as tCode


class StubObjFunct():

	def __init__(self, observers=None):
		self.observers = list(observers) if observers is not None else list()

	def addObjValObserver(self, observer):
		self.observers.append(observer)

	def __call__(self,coeffs):
		objVal = sum(coeffs)
		for x in self.observers:
			x.updateObjVal(objVal)
		
class TestPrinterIterNumbVsObjVal(unittest.TestCase):

	def setUp(self):
		self.printEveryNIters = 1
		self.coeffsA = [ [1,2], [3,4], [5,6] ]
		self.expObjVals = [ sum(x) for x in self.coeffsA ] 
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjFunctA = StubObjFunct()
		tCode.makeObjFunctPrintObjValEveryNIters(self.printEveryNIters, self.testObjFunctA)

	@mock.patch("gen_basis_helpers.fit_cp2k_basis.obj_funct_printers.PrintIterNumberVsObjVal._printVals")
	def testFirstIterCallsPrintFunct(self, mockedPrintVals):
		expIterNumb, expObjVal = 1, self.expObjVals[0]
		self.testObjFunctA(self.coeffsA[0])
		mockedPrintVals.assert_called_with(expIterNumb, expObjVal)

	@mock.patch("gen_basis_helpers.fit_cp2k_basis.obj_funct_printers.PrintIterNumberVsObjVal._printVals")
	def testPrintOnlyCalledForNthIter(self, mockedPrintVals):
		self.printEveryNIters = 3
		self.createTestObjs()
		for coeffs in self.coeffsA:
			self.testObjFunctA(coeffs)
		mockedPrintVals.assert_called_once_with(3,self.expObjVals[2])

