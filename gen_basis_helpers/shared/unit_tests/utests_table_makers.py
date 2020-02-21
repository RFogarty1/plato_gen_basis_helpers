
import unittest
import unittest.mock as mock

import gen_basis_helpers.shared.table_maker_base as tCode

class TestStandardTableMakerBasic(unittest.TestCase):

	def setUp(self):
		self.testDataA = [ [1,1], [2,4], [3,9] ]
		self.fmt = None
		self.headers = None
		self.mapFunct = tCode.blankMapFunct #Means our interface is very much like tabulate, which makes testing easier
		self.createTestObjs()

	def createTestObjs(self):
		self.testTableMakerA = tCode.TableMakerStandard(fmt=self.fmt, headers=self.headers,
		                                                mapInputToRowsFunct=self.mapFunct)


	#The most simple test function possible in effect
	@mock.patch("gen_basis_helpers.shared.table_maker_base.tabulate.tabulate")
	def testPassesSimplestDataToTabulate(self, mockedTabulate):
		self.testTableMakerA.createTable( self.testDataA )
		mockedTabulate.assert_called_once_with(self.testDataA)

	#These functions essentially test the temporary setting of parameters upon createTable function call
	@mock.patch("gen_basis_helpers.shared.table_maker_base.tabulate.tabulate")
	def testFmtParameterPassedWhenSetOnCall(self, mockedTabulate):
		testFmt = "html"
		self.testTableMakerA.createTable( self.testDataA, fmt=testFmt )
		mockedTabulate.assert_called_once_with( self.testDataA, tablefmt=testFmt )

	def testFmtParamerResetsAfterCreateTableCall(self):
		startFmt = self.testTableMakerA.fmt
		testFmt = "github"
		self.assertNotEqual(startFmt,testFmt)
		self.testTableMakerA.createTable( self.testDataA, fmt=testFmt )
		self.assertEqual(startFmt, self.testTableMakerA.fmt)

	#These functions really test specific input options
	@mock.patch("gen_basis_helpers.shared.table_maker_base.tabulate.tabulate")
	def testFmtParameterPasedToTabulate(self, mockedTabulate):
		self.fmt = "html"
		self.createTestObjs()
		self.testTableMakerA.createTable( self.testDataA )
		mockedTabulate.assert_called_once_with( self.testDataA, tablefmt=self.fmt )

	@mock.patch("gen_basis_helpers.shared.table_maker_base.tabulate.tabulate")
	def testHeadersParameterPassedToTabulate(self, mockedTabulate):
		self.headers = ["hi","there"]
		self.createTestObjs()
		self.testTableMakerA.createTable( self.testDataA )
		mockedTabulate.assert_called_once_with( self.testDataA, headers=self.headers )


class TestStandardTableMakerDefaultMapping(unittest.TestCase):
	""" Functions to test that we correctly translate the default input into that required for tabulate"""

	def setUp(self):
		self.dSetA = [ [1,5], [2,4], [3,4], [4,4], [5,7] ]
		self.dSetB = [ [1,4], [2,5], [3,2], [4,5], [5,8] ]
		self.dSetC = [ [1,3], [2,3], [3,3], [4,6], [5,9] ]
		self.fmt = None
		self.headers = None
		self.mapFunct = None
		self.createTestObjs()

	def createTestObjs(self):
		self.testData = [ self.dSetA, self.dSetB, self.dSetC ]
		self.testTableMakerA = tCode.TableMakerStandard(fmt=self.fmt, headers=self.headers,
		                                                mapInputToRowsFunct=self.mapFunct)

	@mock.patch("gen_basis_helpers.shared.table_maker_base.tabulate.tabulate")
	def testExpectedValsPassedToTabulateA(self, mockedTabulate):
		self.testTableMakerA.createTable(self.testData)
		expVals = [ [1,5,4,3], [2,4,5,3], [3,4,2,3], [4,4,5,6], [5,7,8,9] ]
		mockedTabulate.assert_called_once_with( expVals )

	def testRaisesWhenDataSetsDiffLengthsSecondDataSet(self):
		self.dSetB.append( (4,5) ) #Force data to be diff lengths. 
		self.createTestObjs()
		with self.assertRaises(AssertionError):
			self.testTableMakerA.createTable(self.testData)

	def testRaisesWhenDataSetsDiffLengthsFirstDataSet(self):
		""" The first column is qualitatively different from the others, since its the source of xVals """
		self.dSetA.append( (5,5) )
		self.createTestObjs()
		with self.assertRaises(AssertionError):
			self.testTableMakerA.createTable(self.testData)

	def testRaisesWhenXValsNotShared(self):
		self.dSetB[0][0] += 2
		self.createTestObjs()
		with self.assertRaises(AssertionError):
			self.testTableMakerA.createTable(self.testData)


