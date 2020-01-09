
import os
import types
import plato_fit_integrals.core.workflow_coordinator as wFlowBase

from ..shared import misc_utils as misc


""" Goal of these objects is to make it simpler to calculate many-body e0 exchange-correlation contributions to the energy of a system """


class ManyBodyE0SingleCalcObj():
	"""Object that handles running a simple calculation for many body or 2-body E0. Must generate files, run command and energy (after job run)

	"""

	@property
	def e0(self):
		""" The e0 energy for this object (units implementation-dependent)"""
		raise NotImplementedError

	def writeFiles(self):
		raise NotImplementedError

	@property
	def workFolder(self):
		raise NotImplementedError()

	@workFolder.setter
	def workFolder(self, value):
		raise NotImplementedError()

	@property
	def runComm(self):
		raise NotImplementedError()


class PlatoManyBodyE0SingleCalcObj(ManyBodyE0SingleCalcObj):

	def __init__(self, platoCalcObj):
		""" Initialiser
		
		Args:
			platoCalcObj: (CalcObj, see shared/plato_calc_objs.py), a Plato calc-obj. This contains various useful functions such as that for writing the file. This class mainly is a wrapper for platoCalcObj
		
		"""
		self._calcObj = platoCalcObj

	def writeFiles(self):
		self._calcObj.writeFile()

	@property
	def workFolder(self):
		return os.path.split(self._calcObj.filePath)[0] 

	@property
	def runComm(self):
		return [self._calcObj.getRunComm()]

	@property
	def e0(self):
		parsedFile = self._calcObj.parseOutFile()
		e0Coh = parsedFile["energies"].e0Coh
		return e0Coh


class ManyBodyCorrWorkflow(wFlowBase.WorkFlowBase):

	def __init__(self, fullCalcObj, twoBodyCalcObj):
		""" Initializer
		
		Args:
			fullCalcObj: (ManyBodyE0SingleCalcObj) Object used to run calculation/get E0 energy for the full many-body calculation
			twoBodyCalcObj: (ManyBodyE0SingleCalcObj) Object used to run calculation/get E0 energy for the two-body calculation
 
		Returns
			What Function Returns
	 
		Raises:
			Errors
		"""
		self.fullCalcObj = fullCalcObj
		self.twoBodyCalcObj = twoBodyCalcObj
		self.output = types.SimpleNamespace()
		self._assertWorkFoldersAreTheSame()
		self._writeFilesUponInitialisation()

	def _writeFilesUponInitialisation(self):
		self.fullCalcObj.writeFiles()
		self.twoBodyCalcObj.writeFiles()

	@property
	def preRunShellComms(self):
		outList = list()
		outList.extend(self.fullCalcObj.runComm)
		outList.extend(self.twoBodyCalcObj.runComm)
		return outList

	@property
	def namespaceAttrs(self):
		return ["manyBodyE0Xc"]


	def run(self):
		energyDiff = self.fullCalcObj.e0 - self.twoBodyCalcObj.e0
		self.output.manyBodyE0Xc = energyDiff


	#Not using os.path.abspath since i want to avoid mocking it out; Shouldnt be hard to build abspath into the lower level objects
	def _assertWorkFoldersAreTheSame(self):
		manyBodyFolder, twoBodyFolder = self.fullCalcObj.workFolder, self.twoBodyCalcObj.workFolder
		assert manyBodyFolder == twoBodyFolder, "The same folder needs to be used to many-body and two-body calculations; 2-body={}\nmany body = {}".format(manyBodyFolder, twoBodyFolder)


class CompositePerfectCrystalManyBodyE0Calculator():
	volumes = misc.StandardComponentDescriptor("volumes")
	runComms = misc.StandardComponentDescriptor("runComms")
	label = misc.StandardComponentDescriptor("label")

	@misc.getObjectsWithComponentsInstanceWrapper(isComposite=True)
	def __init__(self, e0Calculators):
		""" Initializer 
		
		Args:
			e0Calculators: (iter of PerfectCrystalManyBodyE0Calculator objects)
		"""
		self.objs = list(e0Calculators)


	def createAnalyser(self):
		outList = list()
		for obj in self.objs:
			outList.append( obj.createAnalyser() )
		return CompositePerfectCrystalManyBodyE0Analyser( outList ) #Always returns a single object

class PerfectCrystalManyBodyE0Calculator():

	@misc.getObjectsWithComponentsInstanceWrapper(isComposite=False)
	def __init__(self, manyBodyWorkFlow, volume, label):
		""" Initializer
		
		Args:
			manyBodyWorkFlow: (ManyBodyCorrWorkflow object)
			volume: (float) Volume per atom of the current crystal 
			label: (StandardLabel object) Contains strings with element, structure and method used in this calculation
		"""

		self._workFlow = manyBodyWorkFlow
		self._volume = volume
		self._label = label

	@property
	def volumes(self):
		return [self._volume]

	@property
	def runComms(self):
		return self._workFlow.preRunShellComms

	@property
	def label(self):
		return [self._label]

	def createAnalyser(self):
		self._workFlow.run() #Figure out the many body energy
		mbEnergy = self._workFlow.output.manyBodyE0Xc
		return PerfectCrystalManyBodyE0Analyser( self._volume, mbEnergy, self._label )

class CompositePerfectCrystalManyBodyE0Analyser():

	data  = misc.StandardComponentDescriptor("data")
	label = misc.StandardComponentDescriptor("label")

	@misc.getObjectsWithComponentsInstanceWrapper(isComposite=True)
	def __init__(self, objs):
		self.objs = list(objs)


class PerfectCrystalManyBodyE0Analyser():

	@misc.getObjectsWithComponentsInstanceWrapper(isComposite=False)
	def __init__(self, volume, energy, label):
		self._volume = volume
		self._result = energy
		self._label = label 

	@property
	def data(self):
		return [(self._volume, self._result)]

	@property
	def label(self):
		return [self._label]






