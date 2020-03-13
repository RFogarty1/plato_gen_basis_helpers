

""" File contains names to register for mg basis sets"""

from .. import basis_register as basRegister

from .. import cp2k_basis_obj as basObj


#Rc=7.0 basis set
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("GTH-PBE-q2-spd-2z-rc7pt0-r05pt5-1".lower())
def _getMgRc7pt0BasisObj():
	return basObj.CP2KBasisObjStandard( element="will_be_overwritten", basis="spd-2z-rc7pt0-r05pt5-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-spd-2z-rc7pt0-r05pt5-plus-molopt-sp-1".lower())
def _getMgRc7pt0PlusMoloptSPBasisObj():
	return basObj.CP2KBasisObjStandard( element="will_be_overwritten", basis="spd-2z-rc7pt0-r05pt5-plus-molopt-sp-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-spd-2z-rc7pt0-r05pt5-plus-rc9pt0-sp-1".lower())
def _getMgRc7pt0Plus9pt0SPBasisObj():
	return basObj.CP2KBasisObjStandard( element="will_be_overwritten", basis="spd-2z-rc7pt0-r05pt5-plus-rc9pt0-sp-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-spd-2z-rc7pt0-r05pt5-plus-rc12pt0-sp-1".lower())
def _getMgRc7pt0Plus12pt0SPBasisObj():
	return basObj.CP2KBasisObjStandard( element="will_be_overwritten", basis="spd-2z-rc7pt0-r05pt5-plus-rc12pt0-sp-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )


@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-spd-2z-rc7pt0-r05pt5-plus-diffuse-rc8pt0-sp-1".lower())
def _getMgRc7pt0PlusDiffuse8pt0SPBasisObj():
	return basObj.CP2KBasisObjStandard( element="will_be_overwritten", basis="spd-2z-rc7pt0-r05pt5-plus-diffuse-rc8pt0-sp-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-spd-2z-rc7pt0-r05pt5-plus-diffuse-rc8pt5-sp-1".lower())
def _getMgRc7pt0PlusDiffuse8pt5SPBasisObj():
	return basObj.CP2KBasisObjStandard( element="will_be_overwritten", basis="spd-2z-rc7pt0-r05pt5-plus-diffuse-rc8pt5-sp-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-spd-2z-rc7pt0-r05pt5-plus-diffuse-rc9pt0-sp-1".lower())
def _getMgRc7pt0PlusDiffuse9pt0SPBasisObj():
	return basObj.CP2KBasisObjStandard( element="will_be_overwritten", basis="spd-2z-rc7pt0-r05pt5-plus-diffuse-rc9pt0-sp-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-spd-2z-rc7pt0-r05pt5-plus-diffuse-rc10pt0-sp-1".lower())
def _getMgRc7pt0PlusDiffuse10pt0SPBasisObj():
	return basObj.CP2KBasisObjStandard( element="will_be_overwritten", basis="spd-2z-rc7pt0-r05pt5-plus-diffuse-rc10pt0-sp-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )







@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("GTH-PBE-q2-spd-2z-rc7pt5-r06pt0-1".lower())
def _getMgRc7pt5BasisObj():
	return basObj.CP2KBasisObjStandard( element="will_be_overwritten", basis="spd-2z-rc7pt5-r06pt0-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-spd-2z-rc7pt5-r06pt0-plus-molopt-sp-1".lower())
def _getMgRc7pt5PlusMoloptSPBasisObj():
	return basObj.CP2KBasisObjStandard( element="will_be_overwritten", basis="spd-2z-rc7pt5-r06pt0-plus-molopt-sp-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )




#Rc = 7.5 basis set
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("GTH-PBE-q2-spd-2z-rc8pt0-r06pt5-1".lower())
def _getMgRc8pt0BasisObj():
	return basObj.CP2KBasisObjStandard( element="will_be_overwritten", basis="spd-2z-rc8pt0-r06pt5-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-spd-2z-rc8pt0-r06pt5-plus-molopt-sp-1".lower())
def _getMgRc8pt0PlusMoloptSPBasisObj():
	return basObj.CP2KBasisObjStandard( element="will_be_overwritten", basis="spd-2z-rc8pt0-r06pt5-plus-molopt-sp-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

