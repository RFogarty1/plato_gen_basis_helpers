

import os
import unittest
import unittest.mock as mock

import gen_basis_helpers.misc.shared_io as tCode

class TestJsonIoForObjWithToAndFromDict(unittest.TestCase):

	def setUp(self):
		self.outDictA = {"key_a":"val_a", "key_b":"val_b"}
		self.filenameA = "temp_file.json"
		self.createTestObjs()

	def tearDown(self):
		os.remove(self.filenameA)

	def createTestObjs(self):
		self.testObjA = _TempClassWithToAndFromDictInterface(self.outDictA)

	def testReadAndWriteConsistentA(self):
		tCode.dumpObjWithToDictToJson(self.testObjA, self.filenameA)
		objB = tCode.readObjWithFromDictFromJsonFile(_TempClassWithToAndFromDictInterface, self.filenameA)
		self.assertEqual( self.testObjA, objB )



class _TempClassWithToAndFromDictInterface():

	def __init__(self, inpDict):
		self.inpDict = inpDict

	@classmethod
	def fromDict(cls, inpDict):
		return cls(inpDict)

	def toDict(self):
		return self.inpDict

	def __eq__(self, other):
		return self.inpDict==other.inpDict

