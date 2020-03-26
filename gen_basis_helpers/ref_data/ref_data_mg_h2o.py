
import os
from . import helpers_ref_data as helpers
from ..shared import config_vars as configVars

BASE_FOLDER = os.path.join( configVars.CASTEP_DB_PATH, "mg_h2o" )




def getPlaneWaveGasPhaseGeom(nH2o, charge=2):
	""" Get the optimised plane-wave geometry for an Mg-H2O system in the gas phase
	
	Args:
		nH2o: (int) Number of water molecules
		charge: (int) Total charge on the cluster
			 
	Returns
		 outGeom: (plato_pylib UnitCell object) Geometry of the cluster
 
	"""
	if (nH2o==6) and (charge==2):
		return _getMgH2o6PlaneWaveGeom_2pls()
	else:
		raise ValueError("Geometry unavailable for nH2o={} and charge={}".format(nH2o,charge))


def _getMgH2o6PlaneWaveGeom_2pls():
	filePath = os.path.join(BASE_FOLDER, "mg_h2o6","geom","mg_h2o6_geom_opt.castep")
	return helpers.getUCellInBohrFromCastepOutFile(filePath)


def getPlaneWaveGasPhaseTotalEnergy(nH2o, charge=2):
	""" Get the total energy for an Mg-H2O system in the gas phase from a plane-wave calculation
	
	Args:
		nH2O: (int) Number of H2O molecules
		charge: (int) Total charge on the cluster
			 
	Returns
		 outEnergy: (float) Total energy of the cluster
 
	"""
	if (nH2o==6) and (charge==2):
		return _getMgH2o6PlaneWaveTotalEnergy_2pls()
	else:
		raise ValueError("Geometry unavailable for nH2o={} and charge={}".format(nH2o,charge))
	

def _getMgH2o6PlaneWaveTotalEnergy_2pls():
	filePath = os.path.join(BASE_FOLDER, "mg_h2o6","geom","mg_h2o6_geom_opt.castep")
	return helpers.getEnergyFromCastepOutFile(filePath)
