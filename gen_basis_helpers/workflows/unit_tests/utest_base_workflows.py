
import unittest


import gen_basis_helpers.workflows.base_flow as tCode


class DudWorkflow(tCode.BaseLabelledWorkflow):
	def __init__(self, runComms=None, workFolder=None,
	             output=None, label=None):
		self._runComms = runComms
		self._workFolder = workFolder
		self._label = label
		self.outputDict = output
		self._namespaceAttrs = [ [key for key in self.outputDict.keys()] ]

	def run(self):
		self.output = self._runFunct

	@property
	def namespaceAttrs(self):
		return self._namespaceAttrs

	@property
	def label(self):
		return self._label

	@property
	def preRunShellComms(self):
		return self._runComms


class TestStandardLabelledWorkflowComposite(unittest.TestCase):

	def setUp(self):
		#We use three leaf objevts here
		self.leafRunComms = [x for x in [["runComm_a"], list(), ["runComm_c"]] ]
		self.leafRunComms[0].append("run_comm_a2")
		self.leafLabels = ["a","b","c"]
		self.workFolders = ["fake_folder_{}".format(x) for x in ["a","b","c"]]
		self.outputDictA = {"attrA":"output_a"}
		self.outputDictB = {"attrB":"output_b"}
		self.outputDictC =  {"attrB":"output_c"}

		self.createTestObjs()

	def createTestObjs(self):
		#Create the leaf objects
		outputDicts = [self.outputDictA, self.outputDictB, self.outputDictC]
		leafNames = ["testLeafObj{}".format(x) for x in ["A","B","C"]]
		for idx,name in enumerate(leafNames):
			currObj = DudWorkflow( runComms=self.leafRunComms[idx], workFolder=self.workFolders[idx],
			                        output=outputDicts[idx], label=self.leafLabels[idx] )
			setattr(self, name, currObj)


		#Create a couple of composites
		self.testCompObjA = tCode.StandardLabelledWorkflowComposite( [self.testLeafObjA, self.testLeafObjB] )
		self.testCompObjB = tCode.StandardLabelledWorkflowComposite( [self.testCompObjA, self.testLeafObjC] )

	def testExpectedNamespaceAttrsCompA(self):
		expAttrs = [["attrA"],["attrB"]]
		actAttrs = self.testCompObjA.namespaceAttrs
		self.assertEqual(expAttrs, actAttrs)

	def testExpectedLabelsAttrsCompA(self):
		expLabels = self.leafLabels[:2]
		actLabels = self.testCompObjA.label
		self.assertEqual(expLabels,actLabels)

	def testExpRunCommsCompB(self):
		expRunComms = ["runComm_a","run_comm_a2", "runComm_c"]
		actRunComms = self.testCompObjB.preRunShellComms
		self.assertEqual(expRunComms, actRunComms)








