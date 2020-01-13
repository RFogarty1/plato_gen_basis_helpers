from . import two_cent_2019_paper_method_names as twoCentPaperMap



def createMethodStrToMethodNameDict_three_centre_2019_tight_binding_paper():

	outDict = twoCentPaperMap.createMethodStrToMethodNameDict_two_centre_2019_tight_binding_paper()

	outDict["dft2_mcwedahop_mcwedaxtal_pp_uniform"] = "H_Mc_XT_Mc_PP_uden"
	outDict["dft2_mcwedahop_exactxtal_exact_e0"] = "H_Mc_XT_Ex_PP_Ex"
	outDict["dft2_2bodyhop_all_else_exact"] = "H_2b_XT_Ex_PP_Ex"
	outDict["dft2_exacthop_mcweda_xtal_exacte0"] = "H_Ex_XT_Mc_PP_Ex"
	outDict["dft2_exacthop_2bxtal_sncorr_exacte0"] = "H_Ex_XT_2b_MbSN_PP_Ex"
	outDict["dft2_exact_e1_pp_uniform_e0"] = "H_Ex_XT_Ex_PP_uden"
	outDict["dft2_mcwedahop_mcwedaxtal_pp"] = "H_Mc_XT_Mc_PP_2b"
	outDict["dft2_mcwedahop_mcwedaxtal_pp_fit_uniform"] = "H_Mc_XT_Mc_PP_fit_uden"
	outDict["dft2_exacthop_2bxtal_exacte0"] = "H_Ex_XT_2b_PP_Ex"
	outDict["dft2_exact_e1_pp"] = "H_Ex_XT_Ex_PP_2b"
	outDict["dft2_2bhopxc_xt_exact_pp_exact"] = "H_2bxc_XT_Ex_PP_Ex"
	outDict["dft2_mcwedahop_xtalsncorrd_exact_e0"] = "H_Mc_XT_2b_MbSN_PP_Ex"
	outDict["dft2_exact_e1_pp_uniform_second_order"] = "H_Ex_XT_Ex_PP_uden2"
	outDict["plane-wave"] = "plane-wave"

	return outDict

