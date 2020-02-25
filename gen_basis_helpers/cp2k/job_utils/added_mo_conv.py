
import os


from .conv_grid import mapWorkflowOutputToArray
from ...workflows import convergers as convFlow
from ...shared import calc_runners as calcRunners
from . import converger_template as convTemp




class AddedMoConverger(convTemp.StandardConvergerStandardInputTemplate):

	def __init__(self, baseFolder, createBasicObjFunct, numbAddedMos, label, mapFunction=None):
		""" Initilizer for creator for StandardInput to run calculations to converge energy w.r.t. number of added MOs used in CP2K
		
		Args:
			baseFolder: (str) Folder to carry out the calculations in
			createBasicObjFunct: f(self) function. See the class createBasicObject docstring for more info
			numbAddedMos: (int iter) Number of added MOs for each
			label: (StandardLabel object) This contains info on things like the structure, compound and method used.
			mapFunction: ( f(StandardInput obj) ) Function to convert workflow output to standardOutput data field. Leaving as None is recommended (in which case the final output SHOULD be an array of numbAddedMos vs deltaE per atom)
	
		"""
		self.baseFolder = baseFolder
		self._createBasicObjFunct = createBasicObjFunct
		self._convObjs = list(numbAddedMos)
		self._label = label
		if mapFunction is None:
			self._mapFunct = mapWorkflowOutputToArray
		else:
			self._mapFunct = mapFunction

	@property
	def convObjs(self):
		return self._convObjs

	def modCalcObjWithConvergenceParams(self, calcObj, convObj):
		calcObj.addedMOs = convObj
		calcObj.basePath = os.path.join( self.baseFolder, "added_mos_{}".format(convObj) )
		return calcObj


	def createStandardInputObjFromCalcObjs(self, calcObjs):
		workFlow = convFlow.GridConvergenceEnergyWorkflow( calcObjs, self.convObjs )
		outObj = calcRunners.StandardInputObj(workFlow, self._label, mapFunction=self._mapFunct)
		return outObj





