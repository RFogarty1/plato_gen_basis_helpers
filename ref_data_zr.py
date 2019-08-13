
""" Provides access to reference structures and data for pure Zr """

import os
import ref_elemental_objs as refEleObjs
import plato_pylib.plato.mod_plato_inp_files as modInp
import plato_pylib.plato.plato_paths as platoPaths

import ref_data_mg as refMg

tb1Model = os.path.join("Test","zr_reg_test")
dft2Model = str(tb1Model)
dftModel = str(tb1Model) #Note havnet even got this stub version yet


def createZrReferenceDataObj():
	basePath = "/media/ssd1/rf614/Work/Plato"
	tb1ModAbs = modInp.getAbsolutePathForPlatoTightBindingDataSet(tb1Model)
	dft2ModelAbs = modInp.getAbsolutePathForPlatoTightBindingDataSet(dft2Model)
	dftModelAbs = modInp.getAbsolutePathForPlatoTightBindingDataSet(dftModel,dtype="dft")
	modelHolder = platoPaths.PlatoModelFolders(tb1Path=tb1ModAbs, dft2Path=dft2ModelAbs, dftPath=dftModelAbs)
	return ZrReferenceDataObj(modelHolder)


#Overall object holding ref Data
class ZrReferenceDataObj(refEleObjs.RefElementalDataBase):

	def __init__(self, modelHolder):
		self._modelHolder = modelHolder

	@property
	def modelFiles(self):
		return self._modelHolder

	def getPlaneWaveGeom(self,key):
		return _get_FAKE_MADE_UP_PLANE_WAVE_GEOM(key)

	#TODO: THIS IS USING THE MG VALUES AT THE MOMENT
	def getStructsForEos(self,key):
		mgStructs = refMg.getUCellsForBulkModCalcs(key)
		for currStruct in mgStructs:
			fCoords = currStruct.fractCoords
			for currCoords in fCoords:
				currCoords[-1] = "Zr"
			currStruct.fractCoords = fCoords
		return mgStructs

	#TODO: THIS IS USING MG VALUES
	def getEosFitDict(self,key,eos="murnaghan"):
		return refMg.getPlaneWaveEosFitDict(key,eos=eos)

def _get_FAKE_MADE_UP_PLANE_WAVE_GEOM(key):
	structTypeToFunct = {"hcp":_getZrPlaneWaveHcpGeomAsUCell}
	return structTypeToFunct[key.lower()]()


def _getZrPlaneWaveHcpGeomAsUCell():
	lattParams = [3.23,3.23,5.15] #Angstroms
	basicHcp = refMg._getPerfectHcpMinimalUCell("Zr")
	basicHcp.setLattParams(lattParams)
	basicHcp.convAngToBohr()
	return basicHcp


