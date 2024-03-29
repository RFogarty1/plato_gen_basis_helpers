

""" File contains names to register for mg basis sets"""

from .. import basis_register as basRegister

from .. import cp2k_basis_obj as basObj


#Single-zeta basis sets (despite the file name)
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("GTH-PBE-q2-spd-1z-rc7pt0-r05pt5-1".lower())
def _getMgRc7pt0BasisObj_singleZeta():
    return basObj.CP2KBasisObjStandard( element="Mg", basis="spd-1z-rc7pt0-r05pt5-1",
                                 potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("GTH-PBE-q2-spd-1z-rc7pt5-r06pt0-1".lower())
def _getMgRc7pt5BasisObj_singleZeta():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="spd-1z-rc7pt5-r06pt0-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("GTH-PBE-q2-spd-1z-rc8pt0-r06pt5-1".lower())
def _getMgRc8pt0BasisObj_singleZeta():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="spd-1z-rc8pt0-r06pt5-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )


#Rc=7.0 basis set
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("GTH-PBE-q2-spd-2z-rc7pt0-r05pt5-1".lower())
def _getMgRc7pt0BasisObj():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="spd-2z-rc7pt0-r05pt5-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q2-sp-2z-spd-1zfit-rc7pt0-r05pt5-1".lower())
def _getMgRc7pt0BasisObj_sp2z_d1z_a():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="sp-2z-spd-1zfit-rc7pt0-r05pt5-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )


@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-spd-2z-rc7pt0-r05pt5-plus-molopt-sp-1".lower())
def _getMgRc7pt0PlusMoloptSPBasisObj():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="spd-2z-rc7pt0-r05pt5-plus-molopt-sp-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-spd-2z-rc7pt0-r05pt5-plus-rc9pt0-sp-1".lower())
def _getMgRc7pt0Plus9pt0SPBasisObj():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="spd-2z-rc7pt0-r05pt5-plus-rc9pt0-sp-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-spd-2z-rc7pt0-r05pt5-plus-rc12pt0-sp-1".lower())
def _getMgRc7pt0Plus12pt0SPBasisObj():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="spd-2z-rc7pt0-r05pt5-plus-rc12pt0-sp-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )


@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-spd-2z-rc7pt0-r05pt5-plus-diffuse-rc8pt0-sp-1".lower())
def _getMgRc7pt0PlusDiffuse8pt0SPBasisObj():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="spd-2z-rc7pt0-r05pt5-plus-diffuse-rc8pt0-sp-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-spd-2z-rc7pt0-r05pt5-plus-diffuse-rc8pt5-sp-1".lower())
def _getMgRc7pt0PlusDiffuse8pt5SPBasisObj():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="spd-2z-rc7pt0-r05pt5-plus-diffuse-rc8pt5-sp-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-spd-2z-rc7pt0-r05pt5-plus-diffuse-rc9pt0-sp-1".lower())
def _getMgRc7pt0PlusDiffuse9pt0SPBasisObj():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="spd-2z-rc7pt0-r05pt5-plus-diffuse-rc9pt0-sp-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-spd-2z-rc7pt0-r05pt5-plus-diffuse-rc10pt0-sp-1".lower())
def _getMgRc7pt0PlusDiffuse10pt0SPBasisObj():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="spd-2z-rc7pt0-r05pt5-plus-diffuse-rc10pt0-sp-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )







@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("GTH-PBE-q2-spd-2z-rc7pt5-r06pt0-1".lower())
def _getMgRc7pt5BasisObj():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="spd-2z-rc7pt5-r06pt0-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )


@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q2-sp-2z-spd-1zfit-rc7pt5-r06pt0-1".lower())
def _getMgRc7pt5BasisObj_sp2z_d1z_a():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="sp-2z-spd-1zfit-rc7pt5-r06pt0-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q2-spd-1z-s-2z-rc7pt5-r06pt0-1".lower())
def _getMgRc7pt5BasisObj_s2z_pd1z_a():
    return basObj.CP2KBasisObjStandard( element="Mg", basis="spd-1z-s-2z-rc7pt5-r06pt0-1",
                                 potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q2-sspd-plus-s-diffuse-single-exp-atom-opt".lower())
def _get_sspd_plsS_single_diffuse_opt():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="rc-7pt5-sspd+s-single-exp-atom-opt",
	                                    potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q2-sspd-plus-s-diffuse-six-exp-atom-opt".lower())
def _get_sspd_plsS_six_diffuse_opt():
	return basObj.CP2KBasisObjStandard( element="will be overwritten", basis="rc-7pt5-sspd+s-six-exp-atom-opt",
	                                    potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )


@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q2-sspd-7pt5-sspd-sharedExp-4gau")
def _unused():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="sspd-rc7pt5-sharedExp-4gauss",
	                                    potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q2-sspd-7pt5-sspd-sharedExp-5gau")
def _unused():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="sspd-rc7pt5-sharedExp-5gauss",
	                                    potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q2-sspd+s-7pt5-sharedExp-4gau-single-exp-diffuse")
def _unused():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="sspd+s-rc7pt5-sharedExp-4gauss-diffuse-single-exp-atom-opt",
	                                    potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q2-sspd+s-7pt5-sharedExp-5gau-single-exp-diffuse")
def _unused():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="sspd+s-rc7pt5-sharedExp-5gauss-diffuse-single-exp-atom-opt",
	                                    potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

#OTHER BASIS SETS I ACTUALLY USE REGULARLY
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-BLYP-q2-sspd-plus-s-diffuse-single-exp-atom-opt".lower())
def _unused():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="rc-7pt5-sspd+s-single-exp-atom-opt",
	                                    potential="GTH-BLYP-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q2-spd-1z-s-2z-rc7pt5-r06pt0-plus-s-diffuse-8pt0".lower())
def _getMgRc7pt5BasisObj_ssspd_diffuse():
    return basObj.CP2KBasisObjStandard( element="Mg", basis="spd-1z-s-2z-rc7pt5-r06pt0-plus-diffuse-s-8pt0",
                                 potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )


@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q2-spd-1z-s-2z-rc7pt5-r06pt0-plus-s-ucl".lower())
def _getMgRc7pt5BasisObj_ssspd_molopt():
    return basObj.CP2KBasisObjStandard( element="Mg", basis="spd-1z-s-2z-rc7pt5-r06pt0-plus-ucl-s-1",
                                 potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )



@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-spd-2z-rc7pt5-r06pt0-plus-molopt-sp-1".lower())
def _getMgRc7pt5PlusMoloptSPBasisObj():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="spd-2z-rc7pt5-r06pt0-plus-molopt-sp-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )


@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q2-ssp-rc7pt5-r06pt0".lower())
def _():
    return basObj.CP2KBasisObjStandard( element="Mg", basis="ssp-rc7pt5-r06pt0-1",
                                 potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )


@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q2-ssp-rc7pt5-r06pt0-plus-s-diffuse-8pt0".lower())
def _():
    return basObj.CP2KBasisObjStandard( element="Mg", basis="ssp-rc7pt5-r06pt0-plus-s-diffuse-8pt0",
	                                    potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q2-ssp-plus-s-diffuse-single-exp-atom-opt".lower())
def _():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="ssp-rc7pt5-r06pt0-plus-s-diffuse-single-exp-atom-opt",
	                                    potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )



#Rc = 8.0 basis set
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("GTH-PBE-q2-spd-2z-rc8pt0-r06pt5-1".lower())
def _getMgRc8pt0BasisObj():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="spd-2z-rc8pt0-r06pt5-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q2-sp-2z-spd-1zfit-rc8pt0-r06pt5-1".lower())
def _getMgRc8pt0BasisObj_sp2z_d1z_a():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="sp-2z-spd-1zfit-rc8pt0-r06pt5-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )


@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q2-sp-2z-spd-1z-2exp-sets-rc8pt0-r06pt5-1".lower())
def _getMgRc8pt0BasisObj_sp2z_d1z_b():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="sp-2z-spd-1z-2exp-sets-rc8pt0-r06pt5-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q2-spd-1z-s-2z-rc8pt0-r06pt5-1".lower())
def _getMgRc8pt0BasisObj_s2z_pd1z_a():
    return basObj.CP2KBasisObjStandard( element="Mg", basis="spd-1z-s-2z-rc8pt0-r06pt5-1",
                                 potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )



@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q2-spd-1z-s-2z-rc8pt0-r06pt5-plus-s-diffuse-8pt0".lower())
def _getMgRc8pt0BasisObj_ssspd_diffuse():
    return basObj.CP2KBasisObjStandard( element="Mg", basis="spd-1z-s-2z-rc8pt0-r06pt5-plus-diffuse-s-8pt0",
                                 potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )


@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q2-spd-1z-s-2z-rc8pt0-r06pt5-plus-s-ucl".lower())
def _getMgRc8pt0BasisObj_ssspd_molopt():
    return basObj.CP2KBasisObjStandard( element="Mg", basis="spd-1z-s-2z-rc8pt0-r06pt5-plus-ucl-s-1",
                                 potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )



@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-spd-2z-rc8pt0-r06pt5-plus-molopt-sp-1".lower())
def _getMgRc8pt0PlusMoloptSPBasisObj():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="spd-2z-rc8pt0-r06pt5-plus-molopt-sp-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

