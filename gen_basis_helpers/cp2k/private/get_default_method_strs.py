
"""Private module, just used to initialise method strs for CP2K. Main initial use is to have some that can be used in unit-tests and never have to change """


#NOTE: The important variable is at the bottom; it needs to be there since the creation functions need to be read (by the interpreter) first

import os
from pycp2k import CP2K



def _createTestCP2KObject():

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

	return outObj


#NOTE: This is a separate function to test obj, since i want to be able to change this one (but NEVER test obj) in the future
def _createStandardSPEObj():
	return _createTestCP2KObject()


#NOTE: These keys HAVE to be in lower case
defaultMethodStrsToObjCreators = {"cp2k_test_object":_createTestCP2KObject,
                                  "spe_standard":_createStandardSPEObj}


