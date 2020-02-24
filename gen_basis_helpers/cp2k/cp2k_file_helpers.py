#!/usr/bin/python3

import os
from pycp2k import CP2K
#import ase.calculators.cp2k as cp2kAseCalc


def createCp2kSinglePointObjFromUCellGeom(uCell, elementBasisInfo:"obj list each with .element,.basis,.potential", **kwargs):
	lattVects, fractCoords = uCell.lattVects, uCell.fractCoords
	outObj = createCp2kSinglePointFromMinimalGeom(lattVects, fractCoords, elementBasisInfo, **kwargs)	
	return outObj


def createCp2kSinglePointFromMinimalGeom(lattVects, fractCoords, elementBasisInfo, **kwargs):
	outObj = createDefaultCp2kCalcObj(**kwargs)
	addSubSysSectionCp2kObj(outObj, lattVects, fractCoords, elementBasisInfo)
	return outObj


#def __getCP2KCalculatorObjFromUCellGeom(uCell, elementBasisInfo:"obj list each with .element,.basis,.potential", **kwargs):
#	""" pass an optDict as **optDict to add whatever options you want to the default. AVOID for now. Untested + likely unworking """
#	initObj = createCp2kSinglePointObjFromUCellGeom(uCell, elementBasisInfo, **kwargs)
#	inpStr = initObj.get_input_string()
#	outCalc = cp2kAseCalc.CP2K(inp=inpStr, command="cp2k_shell.sopt")
##	outCalc = cp2kAseCalc.CP2K(command="cp2k_shell.sopt")
#	return outCalc


def createDefaultCp2kCalcObj(**kwargs):
	kwargs = {k.lower():v for k,v in kwargs.items()}

	outDir = kwargs.get("outdir",None)
	projName = kwargs.get("projname",None)

	outObj = CP2K()
	outObj.working_directory = os.path.abspath(os.getcwd())
	outObj.project_name = "cp2k_file"

	#Shortcutrs for below
	cp2kInput = outObj.CP2K_INPUT
	globSect = cp2kInput.GLOBAL
	forceEval = cp2kInput.FORCE_EVAL_add() #Generally the main section (at least for a SPE calc). Can be given more than once hence the add
	subSys = forceEval.SUBSYS
	dft = forceEval.DFT
	
	#Defining tons of params here
	globSect.Run_type = "ENERGY"
	globSect.Print_level = "MEDIUM"
	forceEval.Method = "Quickstep"

	forceEval.PRINT.FORCES.Section_parameters = "On"

	dft.Basis_set_file_name = "BASIS_SET"	
	dft.Potential_file_name = "GTH_POTENTIALS"

	dft.QS.Eps_default = "1.0E-10"
	dft.MGRID.Ngrids = "4"
	dft.MGRID.Cutoff = "[eV] 5000"
	dft.MGRID.Rel_cutoff = "[eV] 50000"

	dft.XC.XC_FUNCTIONAL.Section_parameters = "PBE"
	dft.KPOINTS.Scheme = "MONKHORST-PACK 1 1 1"

	dft.SCF.Scf_guess = "ATOMIC"
	dft.SCF.Eps_scf = "1.0E-7"
	dft.SCF.Max_scf = "300"
	dft.SCF.Added_mos = "4"	
	dft.SCF.DIAGONALIZATION.Section_parameters = "ON"
	dft.SCF.DIAGONALIZATION.Algorithm = "Standard"
	dft.SCF.MIXING.Section_parameters = "T"
	dft.SCF.MIXING.Method = "BROYDEN_MIXING"
	dft.SCF.MIXING.Alpha = "0.4"
	dft.SCF.MIXING.Nbroyden = "8"
	dft.SCF.SMEAR.Section_parameters = "ON"
	dft.SCF.SMEAR.Method = "FERMI_DIRAC"
	dft.SCF.SMEAR.Electronic_temperature = "[K] 157.9"


	modCp2kObjBasedOnDict(outObj,kwargs)

	return outObj

def modCp2kObjBasedOnDict(cp2kObj, optDict):
	useDict = {k.lower():v for k,v in optDict.items()}

	if useDict.get("kpts",None) is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.KPOINTS.Scheme = "MONKHORST-PACK " + " ".join([str(x) for x in useDict["kpts"]])

	if useDict.get("outdir",None) is not None:
		cp2kObj.working_directory = os.path.abspath( useDict["outdir"] )

	if useDict.get("filename") is not None:
		cp2kObj.project_name = useDict["filename"]

	if useDict.get("gridCutAbs".lower()) is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.MGRID.Cutoff = "[eV] " + str(useDict["gridCutAbs".lower()])

	if useDict.get("gridCutRel".lower()) is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.MGRID.Rel_cutoff = "[eV] " + str(useDict["gridCutRel".lower()])

	if useDict.get("basisfile") is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.Basis_set_file_name = useDict["basisfile"]

	if useDict.get("potfile") is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.Potential_file_name = useDict["potfile"]

	if useDict.get("maxscf") is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.SCF.Max_scf = str( useDict["maxscf"] )

	if useDict.get("addedMOs".lower()) is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.SCF.Added_mos = str( useDict["addedMOs".lower()] )

	if useDict.get("printMOs".lower()) is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.PRINT.MO.Section_parameters = "ON"
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.PRINT.MO.Occupation_numbers = True


def addGeomAndBasisInfoToSimpleCP2KObj(cp2kObj, uCell, elementBasisInfo, section="forceEval".lower()):
	""" Takes cp2kObj and adds in keywords for the basis set and the geometry (subsys section). NOTE: This should only be called ONCE on the object. Also it probably isnt general enough to always work
	
	Args:
		cp2kObj: (pycp2k CP2K object) Which holds all input settings for a CP2K calculation
		uCell: (plato_pylib UnitCell object) Which defines the geometry
		elementBasisInfo: (list of CP2KBasisObjBase objects) Each object defines the basis set/pseudopotential to use for one element type
	
	Returns
		Nothing; modifies the cp2kObj object in place
	
	"""

	addSubSysSectionCp2kObjFromUCell(cp2kObj, uCell, elementBasisInfo)
	basisFileStrs = [x.basisFile for x in elementBasisInfo]
	potFileStrs = list(set([x.potFile for x in elementBasisInfo])) #Should all be identical; CP2K cant support multiple PP files for one calculation

	assert len(potFileStrs)==1, "Only 1 psuedopotential datafile allowed; this is a restriction within CP2K (not just this python code"

	modDict = {"basisfile":basisFileStrs, "potfile":potFileStrs}
	modCp2kObjBasedOnDict(cp2kObj,modDict)

	return None

def addSubSysSectionCp2kObjFromUCell(cp2kObj, uCell, elementBasisInfo:"obj list each with .element,.basis,.potential" , section="forceEval".lower()):
	if section.lower() == "forceEval".lower():
		subSys = cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].SUBSYS	
	else:
		raise ValueError("{} is an invalid option for section".format(section))

	lattVects, fractCoords = uCell.lattVects, uCell.fractCoords
	addSubSysSectionCp2kObj(cp2kObj, lattVects, fractCoords, elementBasisInfo, section=section)

def addSubSysSectionCp2kObj(cp2kObj, lattVects, fractCoords, elementBasisInfo:"obj list each with .element,.basis,.potential" , section="forceEval".lower()):
	if section.lower() == "forceEval".lower():
		subSys = cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].SUBSYS	
	else:
		raise ValueError("{} is an invalid option for section".format(section))

	_addCellSectionToSubSysFromMinimalInterface(lattVects, fractCoords, subSys)
	_addBasisInfoSectionToSubSys(elementBasisInfo,subSys)



def _addCellSectionToSubSysFromMinimalInterface(lattVects, scaledCoords, subSysSection):
	vectForm = "[bohr] {:.8f} {:.8f} {:.8f}"
	subSysSection.CELL.A = vectForm.format( *lattVects[0] )
	subSysSection.CELL.B = vectForm.format( *lattVects[1] )
	subSysSection.CELL.C = vectForm.format( *lattVects[2] )

	formattedAtomList = list()
	for currAtom in scaledCoords:
		formattedAtomList.append( [currAtom[-1]] + currAtom[:3] )

	subSysSection.COORD.Scaled = True
	subSysSection.COORD.Default_keyword = formattedAtomList


def _addBasisInfoSectionToSubSys(elementBasisInfo, subSysSection):
	for x in elementBasisInfo:
		currSect = subSysSection.KIND_add()
		currSect.Section_parameters = x.element
		currSect.Element = x.element
		currSect.Basis_set = x.basis
		currSect.Potential = x.potential


