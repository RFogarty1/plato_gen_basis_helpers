
import itertools
import functools

import numpy as np

from plato_pylib.shared.ucell_class import UnitCell
from plato_pylib.parseOther.parse_castep_files import parseCastepOutfile
from plato_pylib.plato.parse_plato_out_files import parsePlatoOutFile
from plato_pylib.parseOther.parse_cp2k_files import parseCpout


import ase.eos



EV_TO_JOULE = 1.60218e-19
BOHR_TO_METRE = 5.2917724900001e-11
EV_TO_RYD =  1 / 13.6056980659
ANG_TO_BOHR = 1.88973

def getBulkModFromOutFilesAseWrapper(outFileList, **kwargs):
	kwargs = {k.lower():v for k,v in kwargs.items()}
	eosModel = kwargs.get("eos", "murnaghan")
	maxFev = kwargs.get("maxfev", 10000)


	ase.eos.curve_fit = functools.partial(ase.eos.curve_fit, maxfev=maxFev)


	allVols, allEnergies = getVolAndEnergiesForASEFromOutFileList(outFileList,**kwargs)

	#And fit the bulk mod/deal with unit convs
	eos = ase.eos.EquationOfState(allVols, allEnergies, eos=eosModel)
	outDict = dict()
	outDict["v0"], outDict["e0"], outDict["b0"] = eos.fit()
	outDict["b0"] *= getBulkModUnitConv("ev","ang")
	outDict["bulkMod"] = outDict["b0"]
	outDict["eos"] = eos
	outDict["v0"] *= ANG_TO_BOHR**3

	allPlotData = eos.getplotdata()
	

#	outDict["data"] = np.array( [(x,y) for x,y in itertools.zip_longest(allPlotData[-2],allPlotData[-1])] )
	outDict["data"] = np.array( [(x,y) for x,y in itertools.zip_longest(allVols,allEnergies)] )
	outDict["fitdata"] = np.array( [(x,y) for x,y in itertools.zip_longest(allPlotData[4],allPlotData[5])] )


	outDict["data"] = outDict["data"][outDict["data"][:,0].argsort()]
	outDict["fitdata"] = outDict["fitdata"][outDict["fitdata"][:,0].argsort()]

	return outDict


def getVolAndEnergiesForASEFromOutFileList(outFileList, **kwargs):
	fileType = kwargs.get("fileType".lower(),None)
	eAttr = kwargs.get("energyType".lower(), "any")

	#Figure out filetype/parser to use
	if fileType is None:
		if outFileList[0].endswith(".castep"):
			fileType = "castep"
		elif outFileList[0].endswith(".cpout"):
			fileType = "cp2k"
		elif outFileList[0].endswith(".out"):
			fileType = "plato"
		else:
			raise ValueError("{} does not have a recognised file extension".format(outFileList[0]))

	fileTypeToParser = { "castep":parseCastepOutfile,
	                     "cp2k": parseCpout,
	                     "plato": parsePlatoOutFile}
	parser = fileTypeToParser[fileType]


	#Figure out the energy/Volume conversions needed
	fileTypeToVolConv = {"castep":1.0,
	                     "cp2k": 1.0,
	                     "plato": 1.0 * (1/(ANG_TO_BOHR**3))}

	fileTypeToEnergyConv = {"castep": 1.0,
	                        "cp2k": 1.0, #Could well be wrong here
	                        "plato": 1/EV_TO_RYD}

	energyConv = fileTypeToEnergyConv[fileType]
	volConv = fileTypeToVolConv[fileType]

	#Now parse the vol vs energies
	allEnergies, allVols = list(), list()
	for currFile in outFileList:
		parsedFile = parser(currFile)
		nAtoms = parsedFile["numbAtoms"]
		allVols.append( parsedFile["unitCell"].volume * volConv/nAtoms )
		if eAttr=="any": 
			try:
				currEnergy = getattr(parsedFile["energies"],"electronicTotalE")
			except ValueError:
				currEnergy = getattr(parsedFile["energies"], "electronicCohesiveE")
		else:
			currEnergy =  getattr(parsedFile["energies"], eAttr)
		allEnergies.append( (currEnergy/nAtoms)*energyConv )


	return allVols, allEnergies





#Gets the value in GPa
def getBulkModUnitConv(atomicEUnits:str, atomicLengthUnits:str):
	if atomicEUnits.lower() == "ev":
		eUnitConv = EV_TO_JOULE
	elif atomicEUnits.lower() == "ryd":
		eUnitConv = (1/EV_TO_RYD) * EV_TO_JOULE
	else:
		raise ValueError( "atomicEUnits value of {} is invalid".format(atomicEUnits))
	if atomicLengthUnits.lower() == "bohr":
		lengthUnitConv = BOHR_TO_METRE
	elif atomicLengthUnits.lower() == "ang":
		lengthUnitConv = ANG_TO_BOHR * BOHR_TO_METRE 
	else:
		raise ValueError( "atomicLengthUnits value of {} is invalid".format(atomicLengthUnits))

	bulkModAtomicToGPA = (1e-9) *  ( EV_TO_JOULE / (lengthUnitConv**3) )
	return bulkModAtomicToGPA

