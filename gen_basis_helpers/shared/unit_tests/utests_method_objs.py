

import unittest
import unittest.mock as mock

import gen_basis_helpers.shared.method_objs as tCode



class TestStandardParsedFileObj(unittest.TestCase):

	def setUp(self):
		self.energies = mock.Mock()
		self.numbAtoms = mock.Mock()
		self.unitCell = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		self.kwargDict = {"energies":self.energies, "numbAtoms":self.numbAtoms, "unitCell":self.unitCell}

	def testInitFromKwargDict(self):
		testObj = tCode.StandardParsedOutputFile(**self.kwargDict)
		for attr in self.kwargDict.keys():
			objAttr = getattr(testObj, attr)
			self.assertEqual( self.kwargDict[attr], objAttr )



