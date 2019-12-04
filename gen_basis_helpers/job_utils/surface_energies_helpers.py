
import os

import plato_fit_integrals.initialise.create_surf_energies_workflows as surfFlows



#TODO: Make this inherit from a more general class, and simply add the surfaceArea property



class PlatoSurfRunner(surfFlows.SurfaceRunnerBase):

	def __init__(self, calcObj, surfaceArea, energyAttr):
		self._surfaceArea = surfaceArea
		self._calcObj = calcObj
		self.energyAttr = energyAttr

	@property
	def surfaceArea(self):
		return self._surfaceArea

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

