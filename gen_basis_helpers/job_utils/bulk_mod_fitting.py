
import itertools
import functools

import numpy as np

from plato_pylib.shared.ucell_class import UnitCell
from plato_pylib.parseOther.parse_castep_files import parseCastepOutfile
from plato_pylib.plato.parse_plato_out_files import parsePlatoOutFile
from plato_pylib.parseOther.parse_cp2k_files import parseCpout


import ase.eos


import plato_pylib.utils.fit_eos as fitBMod



EV_TO_JOULE = 1.60218e-19
BOHR_TO_METRE = 5.2917724900001e-11
EV_TO_RYD =  1 / 13.6056980659
ANG_TO_BOHR = 1.88973

def getBulkModFromOutFilesAseWrapper(outFileList, **kwargs):
	return fitBMod.getBulkModFromOutFilesAseWrapper(outFileList, **kwargs)


def getVolAndEnergiesForASEFromOutFileList(outFileList, **kwargs):
	return fitBMod.getVolAndEnergiesForASEFromOutFileList(outFileList,**kwargs)





#Gets the value in GPa
def getBulkModUnitConv(atomicEUnits:str, atomicLengthUnits:str):
	return fitBMod.getBulkModUnitConv(atomicEUnits, atomicLengthUnits)

