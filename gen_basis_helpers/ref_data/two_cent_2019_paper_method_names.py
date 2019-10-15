

def createMethodStrToMethodNameDict_two_centre_2019_tight_binding_paper():
	outDict = dict()

	outDict["dft2_exact"] = "LCAO"
	outDict["tb1_2c_only"] = "H_2b_XT_2b_PP_2b"
	outDict["tb1_2bxtal_plus_sncorr"] = "H_2b_XT_2b_MbSN_PP_2b"
	outDict["tb1_snxtal_2bmb"] = "H_2b_XT_2bSN_MbSN_PP_2b"
	outDict["tb1_2c_fittedhop"] = "H_fit_XT_2b_PP_2b"
	outDict["tb1_2bxtal_plus_sncorr_fittedhop"] = "H_fit_XT_2b_MbSN_PP_2b"
	outDict["tb1_snxtal_2bmb_fittedHop"] = "H_fit_XT_2bSN_mbSN_PP_2b"
	outDict["tb1_2bxtal_plus_sncorr_fittedhop_uniform_dens_pp"] = "H_fit_XT_2b_MbSN_PP_2b_uden" 
	outDict["tb1_snxtal_2bmb_fittedHop_uniformdens_pp"] = "H_fit_XT_2bSN_MbSN_PP_2b_uden"
	outDict["tb1_2bxtal_plus_sncorr_fittedhop_uniform_dens_pp_fitted_pp"] = "H_fit_XT_2b_MbSN_PP_fit_uden"
	outDict["tb1_snxtal_2bmb_fittedHop_uniformdens_pp_fittedPP"] = "H_fit_XT_2bSN_MbSN_PP_fit_uden"
	outDict["dft_plato"] = "scf-LCAO"

	return outDict


