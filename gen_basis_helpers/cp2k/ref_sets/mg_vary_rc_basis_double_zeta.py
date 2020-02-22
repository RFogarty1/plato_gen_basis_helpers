

""" File contains names to register for mg basis sets"""

from .. import basis_register as basRegister

from .. import cp2k_basis_obj as basObj

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("GTH-PBE-q2-spd-2z-rc7pt0-r05pt5-1".lower())
def _getMgRc7pt0BasisObj():
	return basObj.CP2KBasisObjStandard( element="will_be_overwritten", basis="spd-2z-rc7pt0-r05pt5-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )



@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("GTH-PBE-q2-spd-2z-rc7pt5-r06pt0-1".lower())
def _getMgRc7pt5BasisObj():
	return basObj.CP2KBasisObjStandard( element="will_be_overwritten", basis="spd-2z-rc7pt5-r06pt0-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )


@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("GTH-PBE-q2-spd-2z-rc8pt0-r06pt5-1".lower())
def _getMgRc8pt0BasisObj():
	return basObj.CP2KBasisObjStandard( element="will_be_overwritten", basis="spd-2z-rc8pt0-r06pt5-1",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="PLATO_BASIS" )

