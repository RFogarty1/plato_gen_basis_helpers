#!/usr/bin/python3

import itertools as it
import os

import plato_pylib.shared.unit_convs as uConvHelp

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
	outModder.finalFuncts.append(_modCp2kObjBasedOnDispNonLocalCorr)
	outModder.finalFuncts.append(_modCp2kObjBasedOnScfOptDict)
	outModder.finalFuncts.append(_modCp2kObjBasedOnSurfDipoleCorrOptDict)
	outModder.finalFuncts.append(_modCp2kObjBasedOnMolecularDynamicsOptDict)
	outModder.finalFuncts.append(_modCp2kObjBasedOnTrajPrintDict)
	outModder.finalFuncts.append(_modCp2kObjBasedOnRestartPrintDict)
	outModder.finalFuncts.append(_modCp2kObjBasedOnAtomicConstraints)
	outModder.finalFuncts.append(_modCp2kObjBasedOnCollectiveVariables)
	outModder.finalFuncts.append(_modCp2kObjBasedOnMetadynamicsOptions)
	outModder.finalFuncts.append(_modCp2kObjBasedOnNudgedBandReplicasSection)
	outModder.finalFuncts.append(_modCp2kObjBasedOnHirshfeldOptions)
	return outModder

def _attachXcFunctionalToModder(modder):
	key = "xcFunctional".lower()
	_attachFunctionToModderInstance(key, _modXcFunctionalInObjStd, modder)


def _modXcFunctionalInObjStd(cp2kObj, val):
	specialKeys = ["optB88_pw92".lower()]
	xcSection = cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.XC
	if val.lower() not in specialKeys:
		xcSection.XC_FUNCTIONAL.Section_parameters = val.upper()
	else:
		useVal = val.lower()
		if useVal=="optB88_pw92".lower():
			xcSection.XC_FUNCTIONAL.Section_parameters = None
			xcSection.XC_FUNCTIONAL.LIBXC.Functional = "XC_GGA_X_OPTB88_VDW"
			xcSection.XC_FUNCTIONAL.PW92_add()
			xcSection.XC_FUNCTIONAL.PW92_list[-1].Section_parameters=True
		else:
			raise ValueError("Dont know what to do with {}".format(val))



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


def _modCp2kObjBasedOnDispNonLocalCorr(cp2kObj, useDict):
	vdwPotPart = cp2kObj.CP2K_INPUT.FORCE_EVAL_list[0].DFT.XC.VDW_POTENTIAL
	def initIfNeeded():
		if len(vdwPotPart.NON_LOCAL_list)==0:
			vdwPotPart.NON_LOCAL_add()
			vdwPotPart.Dispersion_functional = "Non_local".upper()

	if useDict.get("disp_nl_corrtype",None) is not None:
		initIfNeeded()
		vdwPotPart.NON_LOCAL_list[-1].Type = useDict["disp_nl_corrtype"]

	if useDict.get("disp_nl_cutoff",None) is not None:
		initIfNeeded()
		vdwPotPart.NON_LOCAL_list[-1].Cutoff = useDict["disp_nl_cutoff"]

	if useDict.get("disp_nl_kernelFileName".lower(),None) is not None:
		initIfNeeded()
		vdwPotPart.NON_LOCAL_list[-1].Kernel_file_name = useDict["disp_nl_kernelFileName".lower()]

	if useDict.get("disp_nl_verboseOutput".lower(), None) is not None:
		initIfNeeded()
		vdwPotPart.NON_LOCAL_list[-1].Verbose_output = useDict["disp_nl_verboseOutput".lower()]

def _modCp2kObjBasedOnScfOptDict(cp2kObj, useDict):

	scfSection = cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.SCF
	mixSection = scfSection.MIXING
	if useDict.get("scfMixAlpha".lower(),None) is not None:
		mixSection.Alpha = useDict["scfMixAlpha".lower()]
	
	if useDict.get("scfMixMethod".lower(),None) is not None:
		mixSection.Method = useDict["scfMixMethod".lower()]

	if useDict.get("scfMixingOn".lower(),None) is not None:
		currVal = "T" if useDict["scfMixingOn".lower()] else "F"
		mixSection.Section_parameters = currVal

	if useDict.get("scfDiagAlgorithm".lower(),None) is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.SCF.DIAGONALIZATION.Algorithm = useDict["scfDiagAlgorithm".lower()].upper()

	if useDict.get("scfOTMinimizer".lower(),None) is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.SCF.OT.Minimizer = useDict["scfOTMinimizer".lower()]

	if useDict.get("scfOTEnergies".lower(),None) is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.SCF.OT.Energies = useDict["scfOTEnergies".lower()]

	if useDict.get("scfOTRotation".lower(),None) is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.SCF.OT.Rotation = useDict["scfOTRotation".lower()]

	if useDict.get("scfOTStepsize".lower(),None) is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.SCF.OT.Stepsize = useDict["scfOTStepsize".lower()]

	if useDict.get("scfOTPreconditioner".lower(),None) is not None:
		scfSection.OT.Preconditioner = useDict["scfOTPreconditioner".lower()]

	if useDict.get("scfOTEnergyGap".lower(),None) is not None:
		scfSection.OT.Energy_gap = useDict["scfOTEnergyGap".lower()]

	if useDict.get("scfOTSafeDIIS".lower(),None) is not None:
		scfSection.OT.Safe_diis = useDict["scfOTSafeDIIS".lower()]

	if useDict.get("scfGuess".lower(),None) is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.SCF.Scf_guess = useDict["scfGuess".lower()].upper()

	if useDict.get("scfPrintRestart_eachMD".lower(),None) is not None:
		scfSection.PRINT.RESTART.EACH.Md = int(useDict["scfPrintRestart_eachMD".lower()])

	if useDict.get("scfPrintRestart_eachSCF".lower(),None) is not None:
		scfSection.PRINT.RESTART.EACH.Qs_scf = int(useDict["scfPrintRestart_eachSCF".lower()])

	if useDict.get("scfPrintRestart_backupCopies".lower(),None) is not None:
		scfSection.PRINT.RESTART.Backup_copies = int(useDict["scfPrintRestart_backupCopies".lower()])

	if useDict.get("scfPrintRestartHistoryOn".lower(),None) is not None:
		val = "ON" if useDict["scfPrintRestartHistoryOn".lower()] is True else "OFF"
		scfSection.PRINT.RESTART_HISTORY.Section_parameters = val

	if useDict.get("scfPrintRestartHistory_eachMD".lower(),None) is not None:
		scfSection.PRINT.RESTART_HISTORY.EACH.Md = int(useDict["scfPrintRestartHistory_eachMD".lower()])

	if useDict.get("scfPrintRestartHistory_eachSCF".lower(),None) is not None:
		scfSection.PRINT.RESTART_HISTORY.EACH.Qs_scf = int(useDict["scfPrintRestartHistory_eachSCF".lower()])

	if useDict.get("scfPrintRestartHistory_backupCopies".lower(),None) is not None:
		scfSection.PRINT.RESTART_HISTORY.Backup_copies = int(useDict["scfPrintRestartHistory_backupCopies".lower()])

	if useDict.get("scfDiagOn".lower(),None) is not None:
		scfSection.DIAGONALIZATION.Section_parameters = useDict["scfDiagOn".lower()]

	if useDict.get("scfOuterEps".lower(),None) is not None:
		scfSection.OUTER_SCF.Eps_scf = useDict["scfOuterEps".lower()]

	if useDict.get("scfOuterMaxIters".lower(),None) is not None:
		scfSection.OUTER_SCF.Max_scf = useDict["scfOuterMaxIters".lower()]

	if useDict.get("scfMaxIterAfterHistoryFull".lower(),None) is not None:
		scfSection.Max_scf_history = useDict["scfMaxIterAfterHistoryFull".lower()]

	if useDict.get("scfHistoryEps".lower(),None) is not None:
		scfSection.Eps_scf_history = useDict["scfHistoryEps".lower()]

def _modCp2kObjBasedOnSurfDipoleCorrOptDict(cp2kObj, useDict):

	dftSection = cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT
	if useDict.get("surf_dipole_corr_useCorr".lower(),None) is not None:
		dftSection.Surface_dipole_correction = useDict["surf_dipole_corr_useCorr".lower()]

	if useDict.get("surf_dipole_corr_surfDipoleDir".lower(),None) is not None:
		dftSection.Surf_dip_dir = useDict["surf_dipole_corr_surfDipoleDir".lower()].upper()

def _modCp2kObjBasedOnMolecularDynamicsOptDict(cp2kObj, useDict):
	mdSection = cp2kObj.CP2K_INPUT.MOTION.MD

	if useDict.get("mdEnsemble".lower(),None) is not None:
		mdSection.Ensemble = useDict["mdEnsemble".lower()]

	if useDict.get("mdSteps".lower(),None) is not None:
		mdSection.Steps = useDict["mdSteps".lower()]

	if useDict.get("mdTimeStep".lower(),None) is not None:
		mdSection.Timestep = useDict["mdTimeStep".lower()]

	if useDict.get("mdTemperature".lower(),None) is not None:
		mdSection.Temperature = useDict["mdTemperature".lower()]

	if useDict.get("mdThermostatType".lower(),None) is not None:
		mdSection.THERMOSTAT.Type = useDict["mdThermostatType".lower()].upper()

	if useDict.get("mdThermoStatOpts".lower(),None) is not None:
		thermoOpts = useDict["mdThermoStatOpts".lower()]
		thermoOpts.addToPyCp2kObj(cp2kObj)

	if useDict.get("mdPrintKindTemps".lower(),None) is not None:
		val = "ON" if useDict["mdPrintKindTemps".lower()] is True else "OFF"
		mdSection.PRINT.TEMP_KIND.Section_parameters = val
		mdSection.Temp_kind = useDict["mdPrintKindTemps".lower()]

	if useDict.get("mdThermoRegions".lower(),None) is not None:
		thermoRegions = useDict["mdThermoRegions".lower()]
		[x.addToPyCp2kObj(cp2kObj) for x in thermoRegions]

	if useDict.get("mdStartTime".lower(),None) is not None:
		mdSection.Time_start_val = useDict["mdStartTime".lower()]

	if useDict.get("mdStartStep".lower(),None) is not None:
		mdSection.Step_start_val = useDict["mdStartStep".lower()]


def _modCp2kObjBasedOnTrajPrintDict(cp2kObj, useDict):
	motionPart = cp2kObj.CP2K_INPUT.MOTION
	def initIfNeeded():
		if len(motionPart.PRINT_list)==0:
			motionPart.PRINT_add()

	if useDict.get("trajPrintEachMd".lower(),None) is not None:
		initIfNeeded()
		motionPart.PRINT_list[-1].TRAJECTORY.EACH.Md = str( useDict["trajPrintEachMd".lower()] )

	if useDict.get("trajPrintEachScf".lower(),None) is not None:
		initIfNeeded()
		motionPart.PRINT_list[-1].TRAJECTORY.EACH.Qs_scf = str( useDict["trajPrintEachScf".lower()] )


def _modCp2kObjBasedOnRestartPrintDict(cp2kObj, useDict):
	motionPart = cp2kObj.CP2K_INPUT.MOTION

	def initIfNeeded():
		if len(motionPart.PRINT_list)==0:
			motionPart.PRINT_add()

	if useDict.get("restartPrintEachMd".lower(),None) is not None:
		initIfNeeded()
		motionPart.PRINT_list[-1].RESTART.EACH.Md = str( useDict["restartPrintEachMd".lower()] )


def _modCp2kObjBasedOnAtomicConstraints(cp2kObj, useDict):
	constraintPart = cp2kObj.CP2K_INPUT.MOTION.CONSTRAINT

	fixedIndices = useDict.get("atPosConstraint_fixIdxPositions".lower(), None)
	if fixedIndices is not None:
		fixComponents = useDict.get("atPosConstraint_fixComponents".lower())
		assert len(fixedIndices)==len(fixComponents)
		for indices, comps in zip(fixedIndices, fixComponents):
			constraintPart.FIXED_ATOMS_add()
			constraintPart.FIXED_ATOMS_list[-1].List = indices
			constraintPart.FIXED_ATOMS_list[-1].Components_to_fix = comps

def _modCp2kObjBasedOnCollectiveVariables(cp2kObj, useDict):
	colVars = useDict.get("colVars".lower(),None)
	if colVars is not None:
		for var in colVars:
			var.addColVarToSubsys(cp2kObj)

def _modCp2kObjBasedOnMetadynamicsOptions(cp2kObj, useDict):
	metaDynSection = cp2kObj.CP2K_INPUT.MOTION.FREE_ENERGY.METADYN

	def _initPrintListIfNeeded():
		if len(metaDynSection.PRINT_list)==0:
			metaDynSection.PRINT_add()

	if useDict.get("metaDyn_doHills".lower(),None) is not None:
		metaDynSection.Do_hills = useDict["metaDyn_doHills".lower()]

	if useDict.get("metaDyn_hillHeight".lower(),None) is not None:
		metaDynSection.Ww = useDict["metaDyn_hillHeight".lower()]

	if useDict.get("metaDyn_ntHills".lower(),None) is not None:
		metaDynSection.Nt_hills = int(useDict["metaDyn_ntHills".lower()])

	if useDict.get("metaDyn_printColvarCommonIterLevels".lower(),None) is not None:
		_initPrintListIfNeeded()
		metaDynSection.PRINT_list[-1].COLVAR.Common_iteration_levels = int(useDict["metaDyn_printColvarCommonIterLevels".lower()])

	if useDict.get("metaDyn_printHillsCommonIterLevels".lower(),None) is not None:
		_initPrintListIfNeeded()
		metaDynSection.PRINT_list[-1].HILLS.Common_iteration_levels = int(useDict["metaDyn_printHillsCommonIterLevels".lower()])

	if useDict.get("metaDyn_printHills".lower(),None) is True:
		_initPrintListIfNeeded()
		metaDynSection.PRINT_list[-1].HILLS.Section_parameters = "ON"

	if useDict.get("metaVars".lower(),None) is not None:
		for mVar in useDict.get("metaVars".lower()):
			mVar.addMetaVarToPyCp2kObj(cp2kObj)

	if useDict.get("metaDyn_spawnHillsHeight".lower(),None) is not None:
		metaDynSection.SPAWNED_HILLS_HEIGHT.Default_keyword = useDict.get("metaDyn_spawnHillsHeight".lower())

	if useDict.get("metaDyn_spawnHillsPos".lower(),None) is not None:
		metaDynSection.SPAWNED_HILLS_POS.Default_keyword = useDict.get("metaDyn_spawnHillsPos".lower())

	if useDict.get("metaDyn_spawnHillsScale".lower(),None) is not None:
		metaDynSection.SPAWNED_HILLS_SCALE.Default_keyword = useDict.get("metaDyn_spawnHillsScale".lower())

	#This might be doable automatically; but i think there are multiple ways to specify it
	if useDict.get("metaDyn_nHillsStartVal".lower(),None) is not None:
		metaDynSection.Nhills_start_val = useDict.get("metaDyn_nHillsStartVal".lower())

def _modCp2kObjBasedOnNudgedBandReplicasSection(cp2kObj, useDict):
	nebSection = cp2kObj.CP2K_INPUT.MOTION.BAND

	#TODO: May have to be cleverer with the replica list later; but this is fine for now
	# since co-ords are currently the only bit of info we give for the replica subsections
	if useDict.get("nudgedband_replica_coords".lower(),None) is not None:
		nebCoords = useDict["nudgedband_replica_coords".lower()]
		for coords in nebCoords:
			nebSection.REPLICA_add()
			convCoords = list()
			for currCoords in coords:
				convCoords.append(  [x*uConvHelp.BOHR_TO_ANG for x in currCoords] )
			nebSection.REPLICA_list[-1].COORD.Default_keyword = convCoords

	if useDict.get("nudgedband_numbReplicas".lower(),None) is not None:
		nebSection.Number_of_replica = useDict["nudgedband_numbReplicas".lower()]

	if useDict.get("nudgedband_procsPerReplicaEnv".lower(),None) is not None:
		nebSection.Nproc_rep = useDict["nudgedband_procsPerReplicaEnv".lower()]

	if useDict.get("nudgedBand_springConstant".lower(),None) is not None:
		nebSection.K_spring = useDict["nudgedBand_springConstant".lower()]

	if useDict.get("nudgedBand_type".lower(),None) is not None:
		nebSection.Band_type = useDict["nudgedBand_type".lower()]

def _modCp2kObjBasedOnHirshfeldOptions(cp2kObj, useDict):
	hirshSection = cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.PRINT.HIRSHFELD

	if useDict.get("hirshfeld_on".lower(), None) is not None:
		startVal = useDict["hirshfeld_on".lower()]
		if startVal is True:
			useVal = "ON"
		elif startVal is False:
			useVal = "OFF"
		else:
			useVal = startVal
		hirshSection.Section_parameters = useVal

	if useDict.get("hirshfeld_selfConsistent".lower(),None) is not None:
		hirshSection.Self_consistent = useDict["hirshfeld_selfConsistent".lower()]

	if useDict.get("hirshfeld_shapeFunction".lower(),None) is not None:
		hirshSection.Shape_function = useDict["hirshfeld_shapeFunction".lower()].upper()

def _attachFunctionToModderInstance(key, function, instance):
	instance.extraKeys.append(key)
	instance.extraFuncts.append(function)

#Simply all options i messed with pre-refactor and probably a few more random ones added later
def _standardModCp2kObjBasedOnDict(cp2kObj, useDict):

	if useDict.get("kpts",None) is not None:
		val = useDict["kpts"]
		try:
			lowerCaseVal = val.lower()
		except AttributeError:
			cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.KPOINTS.Scheme = "MONKHORST-PACK " + " ".join([str(x) for x in val])
		else:
			if lowerCaseVal=="none":
				cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.KPOINTS.Scheme = val
			else:
				raise ValueError("{} is an invalid value for kpts".format(val))

	if useDict.get("outdir",None) is not None:
		cp2kObj.working_directory = os.path.abspath( useDict["outdir"] )

	if useDict.get("filename") is not None:
		cp2kObj.project_name = useDict["filename"]

	if useDict.get("gridCutAbs".lower()) is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.MGRID.Cutoff = "[eV] " + str(useDict["gridCutAbs".lower()])

	if useDict.get("gridCutRel".lower()) is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.MGRID.Rel_cutoff = "[eV] " + str(useDict["gridCutRel".lower()])

	if useDict.get("nGrids".lower()) is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.MGRID.Ngrids = useDict.get("nGrids".lower())

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

	if useDict.get("scfPrintRestart".lower(),None) is not None:
		val = "ON" if useDict["scfPrintRestart".lower()] is True else "OFF"
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.SCF.PRINT.RESTART.Section_parameters = val

	if useDict.get("qsExtrapolationMethod".lower(),None) is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.QS.Extrapolation = useDict["qsExtrapolationMethod".lower()]

	if useDict.get("qsExtrapolationOrder".lower(),None) is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.QS.Extrapolation_order = useDict["qsExtrapolationOrder".lower()]

	if useDict.get("walltime",None) is not None:
		cp2kObj.CP2K_INPUT.GLOBAL.Walltime = useDict["walltime"]

	if useDict.get("useSmearing".lower(),True) is False:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.SCF.SMEAR.Section_parameters = False

	if useDict.get("prefDiagLib".lower(),None) is not None:
		cp2kObj.CP2K_INPUT.GLOBAL.Preferred_diag_library = useDict["prefDiagLib".lower()].upper()

	if useDict.get("rsGrid_distrib".lower(),None) is not None:
		if len(cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.MGRID.RS_GRID_list)==0:
			cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.MGRID.RS_GRID_add()
		currVal = useDict.get("rsGrid_distrib".lower())
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.MGRID.RS_GRID_list[-1].Distribution_layout = currVal

	if useDict.get("extRestartName".lower(),None) is not None:
		cp2kObj.CP2K_INPUT.EXT_RESTART.Restart_file_name = useDict["extRestartName".lower()]

	if useDict.get("dftInpWfnRestartFilename".lower(),None) is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.Restart_file_name = useDict["dftInpWfnRestartFilename".lower()]

	if useDict.get("printForces".lower(),None) is not None:
		startVal = useDict["printForces".lower()]
		if startVal is True:
			useVal = "ON"
		elif startVal is False:
			useVal = "OFF"
		else:
			useVal = startVal
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].PRINT.FORCES.Section_parameters = useVal

	if useDict.get("printPdos".lower(),None) is not None:
		startVal = useDict["printPdos".lower()]
		if startVal is True:
			useVal = "ON"
		elif startVal is False:
			useVal = "OFF"
		else:
			useVal = startVal
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.PRINT.PDOS.Section_parameters = useVal

	if useDict.get("ldosIndices".lower(),None) is not None:
		for indices in useDict["ldosIndices".lower()]:
			cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.PRINT.PDOS.LDOS_add()
			cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].DFT.PRINT.PDOS.LDOS_list[-1].List = " ".join([str(x) for x in indices])

	if useDict.get("motion_cellopt_constraint",None) is not None:
		cp2kObj.CP2K_INPUT.MOTION.CELL_OPT.Constraint = useDict["motion_cellopt_constraint"]


	if useDict.get("motionPrintForces".lower(),None) is not None:
		if len(cp2kObj.CP2K_INPUT.MOTION.PRINT_list)<1:
			cp2kObj.CP2K_INPUT.MOTION.PRINT_add()

		val = "ON" if useDict["motionPrintForces".lower()] is True else "OFF"
		cp2kObj.CP2K_INPUT.MOTION.PRINT_list[-1].FORCES.Section_parameters = val

	if useDict.get("motionPrintVelocities".lower(),None) is not None:
		if len(cp2kObj.CP2K_INPUT.MOTION.PRINT_list)<1:
			cp2kObj.CP2K_INPUT.MOTION.PRINT_add()

		val = "ON" if useDict["motionPrintVelocities".lower()] is True else "OFF"
		cp2kObj.CP2K_INPUT.MOTION.PRINT_list[-1].VELOCITIES.Section_parameters = val

	if useDict.get("subsysInpVelocities".lower(),None) is not None:
		cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].SUBSYS.VELOCITY.Default_keyword = useDict["subsysInpVelocities".lower()]


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

def addBasisInfoToSimpleCP2KObj(cp2kObj, elementBasisInfo):
	subSys = cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].SUBSYS
	_addBasisInfoSectionToSubSys(elementBasisInfo,subSys)
	basisFileStrs = [x.basisFile for x in elementBasisInfo]
	potFileStrs = list(set([x.potFile for x in elementBasisInfo])) #Should all be identical; CP2K cant support multiple PP files for one calculation

	assert len(potFileStrs)==1, "Only 1 psuedopotential datafile allowed; this is a restriction within CP2K (not just this python code"

	modDict = {"basisfile":basisFileStrs, "potfile":potFileStrs}
	modCp2kObjBasedOnDict(cp2kObj,modDict)

def addGeomInfoToSimpleCP2KObj(cp2kObj, uCell):
	subSys = cp2kObj.CP2K_INPUT.FORCE_EVAL_list[-1].SUBSYS
	_addCellSectionToSubSysFromMinimalInterface(uCell.lattVects, uCell.fractCoords, subSys)


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




