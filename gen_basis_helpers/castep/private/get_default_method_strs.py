

#NOTE: The important variable is at the bottom; it needs to be there since the creation functions need to be read (by the interpreter) first




#cut_off_energy not set by defualt; since the user should ALWAYS be choosing that
def _createTestCastepParamDict():
	outDict = dict()
	outDict["xc_functional"] = "pbe"
	outDict["calculate_hirshfeld"] = "true"
	outDict["calculate_densdiff"] = "false" #Density difference between final and atomic densities
	outDict["mix_spin_amp"] = "2.000000000000000" #Mixing of spin densities in density mixing
	outDict["task"] = "singlepoint"
	outDict["num_backup_iter"] = "2" #Number of geom iterations between backups of restart files
	outDict["grid_scale"] = "2.000000000"
	outDict["popn_calculate"] = "true"
	outDict["calculate_elf"] = "false" #Whether to calculate electron localization function
	outDict["fixed_npw"] = "false" #relevant only to geom opts
	outDict["geom_energy_tol"] = "1.000000000000000e-005" #Tightness of free-energy geometry convergence criteria
	outDict["mix_charge_amp"] = "0.500000000000000" #Mixing amplitude for charge density in density mixing
	outDict["spin_fix"] = "6" #Number of scf iterations where the spin is fixed
	outDict["pdos_calculate_weights"] = "false"
	outDict["metals_method"] = "dm" #dm=density mixing
	outDict["fix_occupancy"] = "false" #false means treat as metal; true means treat as insulator
	outDict["run_time"] = "259100" #Maximum number of seconds to run the job for; it should exit cleanly after this
	outDict["smearing_width"] = "0.0136 eV"
	outDict["mix_history_length"] = "20"
	outDict["calculate_stress"] = "false"
	outDict["charge"] = "0"
	outDict["num_dump_cycles"] = "0"
	outDict["mixing_scheme"] = "pulay"
	outDict["geom_method"] = "bfgs"
	outDict["fine_grid_scale"] = "2.300000"
	outDict["geom_max_iter"] = "200"
	outDict["perc_extra_bands"] = "40"
	outDict["elec_energy_tol"] = "1.000000000000000e-006"
	outDict["page_wvfns"] = "0" #Option relates to wriing wavefunctions to disk (not RAM) to save memory
	outDict["max_scf_cycles"] = "150"
	outDict["spin_polarized"] = "false"
	outDict["finite_basis_corr"] = "0"
	outDict["smearing_scheme"] = "fermidirac"
	outDict["geom_stress_tol"] = "0.05"
	outDict["geom_force_tol"] = "0.01"
	outDict["spin"] = "0"
	outDict["mix_charge_gmax"] = "1.5"
	outDict["opt_strategy"] = "speed"
	outDict["geom_disp_tol"] = "5.000000000000000e-004"
	outDict["mix_spin_gmax"] = "1.5"
	return outDict



def _createStandardSPEParamDict():
	return _createTestCastepParamDict()


defaultMethodStrsToObjCreators = {"castep_test_object":_createTestCastepParamDict,
                                  "spe_standard":_createStandardSPEParamDict}



