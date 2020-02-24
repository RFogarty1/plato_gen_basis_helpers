
import copy
import types
import unittest
import unittest.mock as mock

import gen_basis_helpers.cp2k.job_utils.converger_template as tCode

class TestStandardConvergerCreatorTemplate(unittest.TestCase):

	#This uses some highly coupled functions/values...
	def setUp(self):
		self.convObjs = [1,2]
		self.attrAOnAll = 4
		self.labelToSet = "just_a_str"
		self.initObject = types.SimpleNamespace( **{"attrAOnAll":self.attrAOnAll,"convParam":None,
		                                                               "label":None} )
		self.createBasicObjFunct = lambda x:copy.deepcopy( self.initObject )
		self.modCalcObjWithConvergenceParamsFunct = _stubModCalcObj
		self.attrDict = {"label": self.labelToSet} 
		self.createTestObj()

	def createTestObj(self):
		self.testObj = tCode.StandardConvergerStandardInputTemplate(self.convObjs, self.createBasicObjFunct,
		                                             self.modCalcObjWithConvergenceParamsFunct, self.attrDict)
	@mock.patch("gen_basis_helpers.cp2k.job_utils.converger_template.calcRunners.StandardInputObjComposite")
	def testExpectedOutputObjsCreated(self, mockedStandardInpComposite):
		mockedStandardInpComposite.side_effect = lambda x: x

		#Figure out what we're expecting
		expLabels = [self.labelToSet for x in range(len(self.convObjs))]
		expConvParams = self.convObjs
		expAttrAVals = [self.attrAOnAll for x in range(len(self.convObjs))]

		#Run the creator function
		outputObj = self.testObj.createStandardInput()
		actLabels = [x.label for x in outputObj]
		actConvParams = [x.convParam for x in outputObj]
		actAttrAVals = [x.attrAOnAll for x in outputObj]

		#Check for equality
		self.assertEqual(expLabels, actLabels)
		self.assertEqual(expAttrAVals, actAttrAVals)
		self.assertEqual(expConvParams, actConvParams)

def _stubModCalcObj(instance, calcObj, convObj):
	calcObj.label = instance.label
	calcObj.convParam = convObj
	print("convObj = {}".format(convObj))

	return calcObj

