

import os

import plato_pylib.parseOther.parse_castep_files as parseCastep
import plato_pylib.plato.parse_plato_out_files as parsePlatoOut
import plato_pylib.utils.defects as defects

from plato_fit_integrals.initialise.base_objs import PointDefectRunnerBase

class StandardPlatoDefectRunner(PointDefectRunnerBase):

	def __init__(self, calcObj, energyAttr):
		self._calcObj = calcObj
		self.energyAttr = energyAttr

	def writeFiles(self):
		self._calcObj.writeFile()

	@property
	def workFolder(self):
		wFolder = os.path.split(self._calcObj.filePath)[0] 
		return wFolder

	@property
	def nAtoms(self):
		parsedFile = self._getParsedFile()
		return parsedFile["numbAtoms"]

	@property
	def ePerAtom(self):
		parsedFile = self._getParsedFile()
		nAtoms = parsedFile["numbAtoms"]
		totalEnergy = getattr( parsedFile["energies"], self.energyAttr )
		return totalEnergy / nAtoms

	@property
	def runComm(self):
		return [self._calcObj.getRunComm()]
 
	def _getParsedFile(self):
		return self._calcObj.parseOutFile()


class PlatoDefectRunnerForE0Only(StandardPlatoDefectRunner):
	"""Defect runner for the case where you want to fix E1/E2 and only calculate E0. Example use
is fitting a pair-potential while avoiding calculation of most expensive terms

	"""
	def __init__(self, calcObj, energyAttr="e0Tot", otherEnergies=0):
		""" Initialiser for PlatoDefectRunnerForE0Only.
		
		Args:
			calcObj: gen_basis_helpers PlatoCalcObj class (maybe called something slightly diff)
			energyAttr: The attribute to take from plato_pylib energies class
			otherEnergies: Energy to add to E0. Usually itd be Etot-E0 for an initial calculation
		Returns
			object implementing PointDefectRunnerBase interface (plato_fit_integrals library)
	 
		"""
		self._calcObj = calcObj
		self.energyAttr = energyAttr
		self.otherEnergies = otherEnergies

	@property
	def ePerAtom(self):
		parsedFile = self._getParsedFile()
		nAtoms = parsedFile["numbAtoms"]
		totalEnergy = getattr( parsedFile["energies"], self.energyAttr ) + self.otherEnergies
		return totalEnergy / nAtoms



#Older, less recommended functs
class VacancyFileObj:
	def __init__(self, inpPathNoVac, inpPathVac, outExt):
		self.inpPathNoVac = inpPathNoVac
		self.inpPathVac = inpPathVac
		self.outExt = outExt
	
	def getVacancyEnergy(self, eType=None):
		vacPath = self.inpPathVac + self.outExt
		noVacPath = self.inpPathNoVac + self.outExt
		eVac,nAtomsVac = self._getEnergyAndNumbAtoms(vacPath, eType=eType)
		eNoVac, nAtomsNoVac = self._getEnergyAndNumbAtoms(noVacPath, eType=eType)
		return defects.calcVacancyE(eNoVac, eVac, nAtomsNoVac)
	 
	def _getEnergyAndNumbAtoms(self,inpPath, eType=None):
		if eType is None:
			eType = "electronicTotalE"

		parsedFile = parseOutputFile(inpPath)
		energy, nAtoms = getattr(parsedFile["energies"], eType), parsedFile["numbAtoms"]
		return energy,nAtoms



def parseOutputFile(outPath):
	if outPath.endswith('.castep'):
		outParsed = parseCastep.parseCastepOutfile(outPath)
	elif outPath.endswith('.out'):
		outParsed = parsePlatoOut.parsePlatoOutFile(outPath)
	else:
		raise ValueError("{} is an invalid file extension")

	return outParsed




