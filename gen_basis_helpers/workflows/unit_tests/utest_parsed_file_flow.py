
import types
import unittest
import unittest.mock as mock

import gen_basis_helpers.workflows.parsed_file_workflow as tCode

class TestParsedFileWorkflow(unittest.TestCase):

	def setUp(self):
		self.calcObjA = mock.Mock()
		self.createTestObjs()

	def createTestObjs(self):
		self.testObjA = tCode.ParsedFileWorkflow(self.calcObjA)

	def testWriteInpFileCalledUponInitialisation(self):
		self.calcObjA.writeFile.assert_called_with()

	def testRunGivesExpectedResult(self):
		expParsedFile = mock.Mock()
		expOutput = types.SimpleNamespace( parsedFile=expParsedFile )
		self.calcObjA.parsedFile = expParsedFile

		self.testObjA.run()
		actOutputList = self.testObjA.output

		self.assertTrue( len(actOutputList)==1 )
		self.assertEqual(expOutput, actOutputList[0])

	def testExpectedPreRunShellComms(self):
		expCommA = "commA"
		expComms = [expCommA]
		self.calcObjA.runComm = expCommA
		self.assertEqual( expComms, self.testObjA.preRunShellComms )

