from . import two_cent_2019_paper_method_names as twoCentPaperMap



def createMethodStrToMethodNameDict_three_centre_2019_tight_binding_paper():

	outDict = twoCentPaperMap.createMethodStrToMethodNameDict_two_centre_2019_tight_binding_paper()

	outDict["dft2_mcwedahop_mcwedaxtal_pp_uniform"] = "H_Mc_XT_Mc_PP_uden0"
	outDict["dft2_mcwedahop_exactxtal_exact_e0"] = "H_Mc_XT_Ex_PP_Ex"
	outDict["dft2_2bodyhop_all_else_exact"] = "H_2b_XT_Ex_PP_Ex"
	outDict["dft2_exacthop_mcweda_xtal_exacte0"] = "H_Ex_XT_Mc_PP_Ex"
	outDict["dft2_exacthop_2bxtal_sncorr_exacte0"] = "H_Ex_XT_2b_MbSN_PP_Ex"
	outDict["dft2_exact_e1_pp_uniform_e0"] = "H_Ex_XT_Ex_PP_uden0"
	outDict["dft2_mcwedahop_mcwedaxtal_pp"] = "H_Mc_XT_Mc_PP_2b"
	outDict["dft2_mcwedahop_mcwedaxtal_pp_fit_uniform"] = "H_Mc_XT_Mc_PP_fit_uden0"
	outDict["dft2_exacthop_2bxtal_exacte0"] = "H_Ex_XT_2b_PP_Ex"
	outDict["dft2_exact_e1_pp"] = "H_Ex_XT_Ex_PP_2b"
	outDict["dft2_2bhopxc_xt_exact_pp_exact"] = "H_2bxc_XT_Ex_PP_Ex"
	outDict["dft2_mcwedahop_xtalsncorrd_exact_e0"] = "H_Mc_XT_2b_MbSN_PP_Ex"
	outDict["dft2_exact_e1_pp_uniform_second_order"] = "H_Ex_XT_Ex_PP_uden2"
	outDict["dft2_2bhopxc_xt_2b_pp_2b"] = "H_2bxc_XT_2b_PP_2b"
	outDict["dft2_2bodyhopvna_2bodyhopvnl_else_exact"] = "H_2bvna_2bvnl_XT_Ex_PP_2b"
	outDict["dft2_e1_2bxc_exacte0"] = "H_2bxc_XT_2bxc_PP_Ex"
	outDict["dft2_mcwedahop_mcwedaxtal_exact_e0"] = "H_Mc_XT_Mc_PP_Ex"
	outDict["plane-wave"] = "plane-wave"

	return outDict


def createMethodAliasToMethodNameDict_cp2k_basis_set_2020_paper_pure_mg():
	outDict = dict()

	#Pure Mg basis sets
	outDict["ucl-ssp"]                       = "mDZVP"
	outDict["ucl-sssp"]                      = "mTZVP"
	outDict["ucl-ssspp"]                     = "mTZV2P"
	outDict["rc_7pt0"]                       = "rc7pt0-ssppdd"
	outDict["rc_7pt5"]                       = "rc7pt5-ssppdd"
	outDict["rc_8pt0"]                       = "rc8pt0-ssppdd"
	outDict["rc_7pt5_ssppd_3expSets"]        = "rc7pt5-ssppd"
	outDict["rc_7pt5_spd_s"]                 = "rc7pt5-sspd"
	outDict["rc_7pt5_ssp"]                   = "rc7pt5-ssp"
	outDict["rc_7pt5_1z"]                    = "rc7pt5-spd"
	outDict["rc_7pt5_ssp_pls_s_diff_8pt0"]   = "rc7pt5-ssp+s"
	outDict["rc_7pt5_spd_s_pls_s_diff_8pt0"] = "rc7pt5-sspd+s"
	outDict["rc_8pt0_ssppd_3expSets"]        = "rc8pt0-ssppd"
	outDict["rc_8pt0_1z"]                    = "rc8pt0-spd"
	outDict["rc_8pt0_spd_s"]                 = "rc8pt0-sspd"

	#Mixed mg-oxygen basis sets
	outDict["mg_7pt5_ssspd-O-sp"]            = "Mg-rc7pt5-sspd+s-O-mSZV"
	outDict["mg_7pt5_ssspd-O-ssppd"]         = "Mg-rc7pt5-sspd+s-O-mDZVP"
	outDict["mg_7pt5_ssspd-O-ssspppd"]       = "Mg-rc7pt5-sspd+s-O-mTZVP"
	outDict["mg_7pt5_ssspd-O-ssspppdd"]      = "Mg-rc7pt5-sspd+s-O-mTZV2P"
	outDict["mg_ucl-ssspp-O-sp"]             = "Mg-TZV2P-O-mSZV"
	outDict["mg_ucl-ssspp-O-ssppd"]          = "Mg-TZV2P-O-mDZVP"
	outDict["mg_ucl-ssspp-O-ssspppd"]        = "Mg-TZV2P-O-mTZVP"
	outDict["mg_ucl-ssspp-O-ssspppdd"]       = "Mg-TZV2P-O-mTZV2P"

	#Pure oxygen basis sets
	outDict["O-sp"]                          = "O-mSZV"
	outDict["O-ssppd"]                       = "O-mDZVP"
	outDict["O-ssspppd"]                     = "O-mTZVP"
	outDict["O-ssspppdd"]                    = "O-mTZV2P"

	#mixed mg-oh basis sets
	outDict["mg_7pt5_ssspd-OH-SZV"]         = "Mg-rc7pt5-sspd+s-OH-mSZV"
	outDict["mg_7pt5_ssspd-OH-DZVP"]        = "Mg-rc7pt5-sspd+s-OH-mDZVP"
	outDict["mg_7pt5_ssspd-OH-TZVP"]        = "Mg-rc7pt5-sspd+s-OH-mTZVP"
	outDict["mg_7pt5_ssspd-OH-TZV2P"]       = "Mg-rc7pt5-sspd+s-OH-mTZV2P"
	outDict["mg_ucl-ssspp-OH-SZV"]          = "Mg-TZV2P-OH-mSZV"
	outDict["mg_ucl-ssspp-OH-DZVP"]         = "Mg-TZV2P-OH-mDZVP"
	outDict["mg_ucl-ssspp-OH-TZVP"]         = "Mg-TZV2P-OH-mTZVP"
	outDict["mg_ucl-ssspp-OH-TZV2P"]        = "Mg-TZV2P-OH-mTZV2P"

	return outDict



