#!/usr/bin/python3

import itertools as it
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
	modderFunct = _getStandardPyCp2kModder()
	modderFunct(cp2kObj, optDict)


class Pycp2kModderStandard():
	""" Callable class for modifying a PyCP2K object with a dictionary of update options. See self.modObjBasedOnDict for call interface. Note this class is generally used as a backend.

	"""
	def __init__(self, extraKeys=None, extraFuncts=None, finalFuncts=None):
		""" Initializer
		
		Args:
			extraKeys: (iter of str) Each key represents an input key for the optDict
			extraFuncts: (iter of functs, same length as extraKeys) f(pyCP2KObj, val). Each must modify the cp2k object in place with {key:val} pair from an input dictionary. extraFunct[idx](obj) triggers if extraKeys[idx] found in useDict 
			finalFuncts: (iter of functs) Interface is f(cp2kObj,useDict). Probably more useful/flexible than the others, allows dealing with multiple keys simultaneously
				 
		"""
		self.extraKeys = list() if extraKeys is None else list(extraKeys)
		self.extraFuncts = list() if extraFuncts is None else list(extraFuncts)
		self.finalFuncts = list() if finalFuncts is None else list(finalFuncts)

		assert len(self.extraKeys)==len(self.extraFuncts)

	def modObjBasedOnDict(self, cp2kObj, optDict):
		useDict = {k.lower():v for k,v in optDict.items()}
		_standardModCp2kObjBasedOnDict(cp2kObj, useDict)
		for key,funct in it.zip_longest(self.extraKeys, self.extraFuncts):
			if useDict.get(key.lower(),None) is not None:
				funct(cp2kObj, useDict[key.lower()])

		for funct in self.finalFuncts:
			funct(cp2kObj, useDict)

	def __call__(self, cp2kObj, optDict):
		self.modObjBasedOnDict(cp2kObj,optDict)

def _getStandardPyCp2kModder():
	outModder = Pycp2kModderStandard()
	_attachXcFunctionalToModder(outModder)
	outModder.finalFuncts.append(_modCp2kObjBasedOnGrimmeDisp)
	return outModder

def _attachXcFunctionalToModder(modder):
	key = "xcFunctional".lower()
	def modXcFunctionalInObj(cp2kObj, val):
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.XC.XC_FUNCTIONAL.Section_parameters = val.upper()
	_attachFunctionToModderInstance(key, modXcFunctionalInObj, modder)


def _modCp2kObjBasedOnGrimmeDisp(cp2kObj, useDict):
	ppPart = cp2kObj.CP2K_INPUT.FORCE_EVAL_list[0].DFT.XC.VDW_POTENTIAL
	def initIfNeeded():
		if len(ppPart.PAIR_POTENTIAL_list)==0:
			ppPart.PAIR_POTENTIAL_add()
			cp2kObj.CP2K_INPUT.FORCE_EVAL_list[0].DFT.XC.VDW_POTENTIAL.Dispersion_functional = "Pair_Potential".upper()

	if useDict.get("grimme_disp_corrType".lower(), None) is not None:
		initIfNeeded()
		ppPart.PAIR_POTENTIAL_list[0].Type = useDict["grimme_disp_corrType".lower()]

	if useDict.get("grimme_disp_excludeKindsD3".lower(), None) is not None:
		initIfNeeded()
		ppPart.PAIR_POTENTIAL_list[0].D3_exclude_kind = " ".join([str(x) for x in useDict["grimme_disp_excludeKindsD3".lower()]])

	if useDict.get("grimme_disp_paramFile".lower(), None) is not None:
		initIfNeeded()
		ppPart.PAIR_POTENTIAL_list[0].Parameter_file_name = useDict["grimme_disp_paramfile"]

	if useDict.get("grimme_disp_refFunctional".lower(),None) is not None:
		initIfNeeded()
		ppPart.PAIR_POTENTIAL_list[0].Reference_functional = useDict["grimme_disp_reffunctional"].upper()

	if useDict.get("grimme_disp_printDFTD".lower(), None) is not None:
		initIfNeeded()
		if useDict["grimme_disp_printDFTD".lower()]:
			ppPart.PAIR_POTENTIAL_list[0].PRINT_DFTD.Section_parameters = "ON"


def _attachFunctionToModderInstance(key, function, instance):
	instance.extraKeys.append(key)
	instance.extraFuncts.append(function)

#Simply all options i messed with pre-refactor
def _standardModCp2kObjBasedOnDict(cp2kObj, useDict):

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

	if useDict.get("epsScf".lower()) is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.SCF.Eps_scf = str( useDict["epsScf".lower()] )

	if useDict.get("addedMOs".lower()) is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.SCF.Added_mos = str( useDict["addedMOs".lower()] )

	if useDict.get("printMOs".lower()) is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.PRINT.MO.Section_parameters = "ON"
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.PRINT.MO.Occupation_numbers = True

	if useDict.get("printOverlapConditionDiag".lower(), False) is True:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.PRINT.OVERLAP_CONDITION.Section_parameters = "ON"
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.PRINT.OVERLAP_CONDITION.Diagonalization = True

	if useDict.get("printAOMullikenPop".lower()) is not None:
		if useDict["printAOMullikenPop".lower()]:
			cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.PRINT.MULLIKEN.Section_parameters = "ON"
			cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.PRINT.MULLIKEN.Print_gop = True

	if useDict.get("epsDef".lower(),None) is not None:
		rawNumber = useDict["epsDef".lower()]
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.QS.Eps_default = "{:.2g}".format( rawNumber ).upper()

	if useDict.get("epsGvgRspace".lower(),None) is not None:
		rawNumber = useDict["epsGvgRspace".lower()]
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.QS.Eps_gvg_rspace = "{:.2g}".format( rawNumber ).upper()
	
	if useDict.get("epsPgfOrb".lower(),None) is not None:
		rawNumber = useDict["epsPgfOrb".lower()]
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.QS.Eps_pgf_orb = "{:.2g}".format( rawNumber ).upper()

	if useDict.get("epsPPNL".lower(),None) is not None:
		rawNumber = useDict["epsPPNL".lower()]
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.QS.Eps_ppnl = "{:.2g}".format( rawNumber ).upper()

	if useDict.get("epsRho".lower(), None) is not None:
		rawNumber = useDict["epsRho".lower()]
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.QS.Eps_rho = "{:.2g}".format( rawNumber ).upper()

	if useDict.get("epsCoreCharge".lower(),None) is not None:
		rawNumber = useDict["epsCoreCharge".lower()]
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.QS.Eps_core_charge = "{:.2g}".format( rawNumber ).upper()

	if useDict.get("charge") is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.Charge = useDict["charge"]

	if useDict.get("runType".lower()) is not None:
		cp2kObj.CP2K_INPUT.GLOBAL.Run_type = useDict["runtype"].upper()
		if useDict["runtype"].lower() == "cell_opt":
			cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].Stress_tensor = "analytical".upper()
		if useDict["runtype"].lower() == "bsse":
			for frag in useDict["fragmentsBSSE".lower()]:
				currFrag = cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].BSSE.FRAGMENT_add()
				currFrag.List = " ".join(["{}" for x in frag]).format(*frag) 

	if useDict.get("geo_constrain_cell_angles") is not None:
		angleConstraints = useDict["geo_constrain_cell_angles"]
		if all(angleConstraints):
			cp2kObj.CP2K_INPUT.MOTION.CELL_OPT.Keep_angles = True
		elif not any(angleConstraints): #All False
			pass
		else:
			raise ValueError("Constraining angles to {} is not currently supported".format(angleConstraints))


	if useDict.get("useSmearing".lower(),True) is False:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.SCF.SMEAR.Section_parameters = False


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
		currSect.Section_parameters = x.kind
		currSect.Basis_set = x.basis
		if x.ghost:
			currSect.Ghost = True
		else:
			currSect.Element = x.element
			currSect.Potential = x.potential




