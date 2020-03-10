
import types


from . import base_flow as baseFlow

class SelfPointDefectWorkflow(baseFlow.BaseLabelledWorkflow):
	"""Workflow for calculating formation energies of self-defects. This means self-interstitials or vacancies for pure elements 
	"""

	def __init__(self, defectCalcObj, bulkCalcObj, label=None):
		""" Initializer
		
		Args:
			defectCalcObj: (CalcMethod object) Object used for calculating total energy for the structure with a defect
			bulkCalcObj: (CalcMethod object) Object used for calculating total energy for the bulk system
				
		"""

		self.eType = "electronicTotalE"
		self.defectCalcObj = defectCalcObj
		self.bulkCalcObj = bulkCalcObj
		startDict = {k:None for k in self.namespaceAttrs[0]}
		self.output = [ types.SimpleNamespace(**startDict) ]
		self._writeInpFiles()


	def _writeInpFiles(self):
		self.defectCalcObj.writeFile()
		self.bulkCalcObj.writeFile()

	def run(self):
		defectEnergy = self._defectNumbAtoms * (self._defectEPerAtom-self._bulkEPerAtom)
		self.output[0].bulkEPerAtom = self._bulkEPerAtom
		self.output[0].defectEPerAtom = self._defectEPerAtom
		self.output[0].defectE = defectEnergy

	@property
	def preRunShellComms(self):
		return [self.bulkCalcObj.runComm,self.defectCalcObj.runComm]

	@property
	def namespaceAttrs(self):
		return [ ["bulkEPerAtom","defectEPerAtom","defectE"] ]


	@property
	def _bulkEPerAtom(self):
		return self._getEPerAtomFromCalcObj(self.bulkCalcObj)

	@property
	def _defectEPerAtom(self):
		return self._getEPerAtomFromCalcObj(self.defectCalcObj)

	@property
	def _defectNumbAtoms(self):
		return self.defectCalcObj.parsedFile.numbAtoms

	def _getEPerAtomFromCalcObj(self, calcObj):
		parsedFile = calcObj.parsedFile
		totalEnergy = getattr(parsedFile.energies,self.eType)
		nAtoms = parsedFile.numbAtoms
		return totalEnergy/nAtoms



