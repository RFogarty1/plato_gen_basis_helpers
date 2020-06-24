

import unittest
import unittest.mock as mock

import gen_basis_helpers.shared.plato_calc_objs as tCode

class TestPlatoCalcObj(unittest.TestCase):

	def setUp(self):
		self.filePath = "fake_path"
		self.strDict = mock.Mock()
		self.strDictWriteFunction = mock.Mock()
		self.fileParser = mock.Mock()
		self.runCommFunction = mock.Mock()
		self.inpFileExt = mock.Mock()
		self.outFileExt = "fake_ext"
		self.raiseInParser = True

		self.createTestObjs()

	def createTestObjs(self):
		kwargDict = {"filePath":self.filePath, "strDict":self.strDict,
		             "strDictWriteFunction":self.strDictWriteFunction,
		             "fileParser":self.fileParser, "runCommFunction":self.runCommFunction,
		             "inpFileExt":self.inpFileExt, "outFileExt":self.outFileExt,
		             "raiseInParserIfScfNotConverged":self.raiseInParser}
		self.testObjA = tCode.CalcObj.fromEnforcedKwargs(**kwargDict)

	def testRaiseIfScfNotConverged(self):
		returnedDict = {"scf_is_converged":False}
		self.fileParser.side_effect = [returnedDict]
		with self.assertRaises(ValueError):
			self.testObjA.parseOutFile()


