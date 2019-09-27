
import os
from ..shared import calc_methods as calcMethods
from . import dos_helpers as dosHelp
import plato_pylib.utils.job_running_functs as jobRun


class DosLabel():
	"""Label holding strings used to differentiate different DosRunner/Analyser objects

	Attributes:
		eleKey (str): key representing the element (or possible elements) involved.
		structKey (str): key representing the structure involved
		methodKey (str): key representing the calculation method used

	"""
	def __init__(self, eleKey=None, structKey=None, methodKey=None):

		setattr(self, "eleKey", eleKey)
		setattr(self, "structKey", structKey)
		setattr(self, "methodKey", methodKey)

		self.reqArgs = ["eleKey", "structKey", "methodKey"]
		for x in self.reqArgs:
			if getattr(self,x) is None:
				raise ValueError("{} is a required argument for DosLabel".format(x))



	@property
	def components(self):
		""" List of components contributing to identity of this label. At time of writing this is eleKey, structKey and methodKey. Consistent ordering not garanteed.
		
		"""
		return [getattr(self,x) for x in self.reqArgs]


	def __eq__(self, other):
		for attr in self.reqArgs:
			try:
				if getattr(self,attr) != getattr(other,attr):
					return False
			except AttributeError: #If other is an inappropriate type (e.g an integer) this should trigger
				return False
		return True


	def __hash__(self):
		return hash( tuple([getattr(self,x) for x in self.reqArgs]) )


class DosRunnerBase():
	"""Class for running DoS calculations.

	Attributes:
		dosComms: Bash commands for creating density of states data
		singlePointEnergyComms: Bash commands for running the single-point energy calculations required

	"""

	@property
	def dosComms(self):
		""" iter, each entry is a string for the bash command to create density of states data
		
		"""
		raise NotImplementedError()


	@property
	def singlePointEnergyComms(self):
		""" iter, each entry is a string for the bash command to run the total energy calculation (which is needed to generate Dos Data
		
		"""
		raise NotImplementedError()

	def runSinglePointEnergyCalcs(self, nCores=1):
		""" Runs the single-point energy calculations required to get density of states data
		
		Args:
			nCores: (int) Number of cores to run calculations over (default=1)
				
		Returns
			Nothing
		
		"""
		raise NotImplementedError()


	def runDosGeneratingCalcs(self, nCores=1):
		""" Calculates (and saves somewhere accesible) the density of states data 
		
		Args:
			nCores: (int, optional) Number of cores to run calculations over (default=1)
				
		Returns
			Nothing
		
		"""
		raise NotImplementedError()

	def createAnalyser(self):
		""" Creates a DosAnalyser (possibly composite) object - see DosAnalyserBase or similar
		
		Returns
			DosAnalyser object
		
		"""
		raise NotImplementedError()



class DosRunnerComposite(DosRunnerBase):

	def __init__(self, inpObjs):
		pass


class DosRunnerPlato(DosRunnerBase):

	def __init__(self, platoCalcObj, smearWidth, stepSize, label):
		""" Initialiser for DosRunnerPlato
		
		Args:
			platoCalcObj (calc_methods.CalcObj): object which holds details for running the SPE calculation
			smearWidth (float): Smearing width for the DoS in eV
			stepSize (float): size of steps to take in eV (e.g. 0.1 means output a DoS value for every 0.1 eV)
			label (DosLabel): object which holds labels to help identify/classify/separate different jobs
				
		"""
		self._calcObj = platoCalcObj 
		self._smearWidth = smearWidth
		self._stepSize = stepSize
		self._label = label

	@property
	def dosComms(self):
		return [dosHelp.getDosPlotData(self._baseFilePath +'.occ', self._smearWidth, self._stepSize)]

	@property
	def singlePointEnergyComms(self):
		return [self._calcObj.getRunComm()]

	@property
	def _baseFilePath(self):
		filePath = self._calcObj.filePath
		basePath = os.path.splitext(filePath)[0]
		return basePath

	@property
	def label(self):
		return self.label

	def runSinglePointEnergyCalcs(self, nCores=1):
		self._calcObj.writeFile()
		jobRun.executeRunCommsParralel( self.singlePointEnergyComms, nCores )

	def runDosGeneratingCalcs(self, nCores=1):
		jobRun.executeRunCommsParralel( self.dosComms, nCores )

	def _getDosDosDataDictAfterRunningCalcs(self):
		return dosHelp.getDosPlotData(self._baseFilePath + ".occ", self._smearWidth, self._stepSize, runDos=False)


	def createAnalyser(self):
		dosDataDict = self._getDosDosDataDictAfterRunningCalcs()
		dosData = dosDataDict["dosdata"]
		eFermi = dosDataDict["efermi"]
		return DosAnalyserPlato(dosData, eFermi, self._label)



#TODO: Define the base class
class DosAnalyserPlato():

	def __init__(self, dosData, eFermi, label):
		""" Initialiser for DosAnalyser for plato
		
		Args:
			dosData (np array, nx2): first column is energies (assumed eV), second column is density of states
			eFermi (float): Fermi energy in same units as energies (which are first col of dosData)
			label (DosLabel object):
	
		"""
		self.dosData = dosData
		self.eFermi = eFermi
		self.label = label













class DosOptions():
    def __init__(self, methodStr, kpts, startFolder, uCellStruct, smearWidth, stepSize, modelDataFolder, integGrid = None):
        
        self.methodStr = methodStr
        self.kpts = kpts
        self.startFolder = startFolder
        self.smearWidth = smearWidth
        self.stepSize = stepSize
        self.integGrid = integGrid
        self.struct = uCellStruct
        
        #Create the method
        self.platoMethodObj = calcMethods.createPlatoMethodObj(self.methodStr)
        self.platoMethodObj.kpts = kpts
        self.platoMethodObj.dataSet = modelDataFolder
        if integGrid is not None:
            self.platoMethodObj.integGrid = integGrid

            
    def createDosObj(self):
        inpPath = os.path.abspath( os.path.join(self.startFolder,self.methodStr,"spe_for_dos.in") )
        strDict = self.platoMethodObj.getStrDictWithStruct(self.struct)
        runCommFunct = self.platoMethodObj.runCommFunction
        calcObj =  calcMethods.getPlatoCalcObjFromInpPathAndStrDictAndRunCommFunction(inpPath, strDict, runCommFunct)
        return DosCalcObj(calcObj, self.smearWidth, self.stepSize)
            
#class DosCalcObj():
#    def __init__(self, calcObj:"calcMethods.CalcObj object", smearWidth, stepSize):
#        self.fileCalcObj = calcObj
#        self.smearWidth = smearWidth
#        self.stepSize = stepSize
#    
#    def writeEnergyCalcFile(self):
#        self.fileCalcObj.writeFile()
#    
#    def getRunCommEnergyCalc(self):
#        return self.fileCalcObj.getRunComm()
#    
#    def getRunCommDos(self):
#        return dosHelp.getDosPlotData(self.baseFilePath, self.smearWidth, self.stepSize)
#    
#    def getDosPlotData(self):
#        return dosHelp.getDosPlotData(self.baseFilePath + ".occ", self.smearWidth, self.stepSize, runDos=False)["dosdata"]
#    
#    @property
#    def baseFilePath(self):
#        return self.fileCalcObj.filePath
#    
#    @property
#    def eFermi(self):
#        return dosHelp.getDosPlotData(self.baseFilePath + ".occ", self.smearWidth, self.stepSize, runDos=False)["efermi"]


