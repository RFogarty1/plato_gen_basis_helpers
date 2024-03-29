
from .. import basis_register as basRegister
from .. import cp2k_basis_obj as basObj


#1 basis function
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("H-GTH-PBE-q6-MolOpt-s".lower())
def _getHydrogenMolOptSZVBasis():
	return basObj.CP2KBasisObjStandard( element="H", basis="SZV-MOLOPT-GTH-q1",
	                             potential="GTH-PBE-q1", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT" )

#5 basis functions
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("H-GTH-PBE-q1-MolOpt-DZVP".lower())
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("H-GTH-PBE-q6-MolOpt-ssp".lower())
def _getHydrogenMolOptDZVPBasis():
	return basObj.CP2KBasisObjStandard( element="H", basis="DZVP-MOLOPT-GTH",
	                             potential="GTH-PBE-q1", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("H-GTH-BLYP-q1-MolOpt-DZVP".lower())
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("H-GTH-BLYP-q1-MolOpt-ssp".lower())
def _():
	return basObj.CP2KBasisObjStandard( element="H", basis="DZVP-MOLOPT-GTH",
	                                    potential="GTH-BLYP-q1", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT" )

#6 basis functions
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("H-GTH-PBE-q6-MolOpt-sssp".lower())
def _getHydrogenMolOptTZVPBasis():
	return basObj.CP2KBasisObjStandard( element="H", basis="TZVP-MOLOPT-GTH",
	                             potential="GTH-PBE-q1", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT" )


#9 basis functions
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("H-GTH-PBE-q1-MolOpt-TZV2P".lower())
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("H-GTH-PBE-q6-MolOpt-ssspp".lower())
def _getHydrogenMolOptTZV2PBasis():
	return basObj.CP2KBasisObjStandard( element="H", basis="TZV2P-MOLOPT-GTH",
	                             potential="GTH-PBE-q1", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("H-GTH-BLYP-q1-MolOpt-TZV2P".lower())
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("H-GTH-BLYP-q1-MolOpt-ssspp".lower())
def _():
    return basObj.CP2KBasisObjStandard( element="H", basis="TZV2P-MOLOPT-GTH",
                                        potential="GTH-BLYP-q1", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT" )

#14 basis functions
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("H-GTH-PBE-q6-MolOpt-sssppd".lower())
def _getHydrogenMolOptTZV2PXBasis():
	return basObj.CP2KBasisObjStandard( element="H", basis="TZV2PX-MOLOPT-GTH",
	                             potential="GTH-PBE-q1", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT" )



#GTH (non-molopt) basis sets
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("H-GTH-PADE-q1-DZVP".lower())
def _():
	return basObj.CP2KBasisObjStandard( element="H", basis="DZVP-GTH-q1",
	                                    potential="GTH-PADE-q1", potFile="GTH_POTENTIALS", basisFile="GTH_BASIS_SETS" )


@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("H-GTH-PBE-q1-DZVP".lower())
def _():
	return basObj.CP2KBasisObjStandard( element="H", basis="DZVP-GTH-q1",
	                                    potential="GTH-PBE-q1", potFile="GTH_POTENTIALS", basisFile="GTH_BASIS_SETS" )


@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("H-GTH-PBE-q1-TZV2P".lower())
def _():
	return basObj.CP2KBasisObjStandard( element="H", basis="TZV2P-GTH-q1",
	                                    potential="GTH-PBE-q1", potFile="GTH_POTENTIALS", basisFile="GTH_BASIS_SETS" )

