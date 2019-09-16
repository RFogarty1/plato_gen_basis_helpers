


import plato_pylib.parseOther.parse_castep_files as parseCastep
import plato_pylib.plato.parse_plato_out_files as parsePlatoOut
import plato_pylib.utils.defects as defects

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

