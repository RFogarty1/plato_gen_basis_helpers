


from .. import basis_register as basRegister
from .. import cp2k_basis_obj as basObj

#4 basis functions
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("O-GTH-PBE-q6-MolOpt-sp".lower())
def _getOxygenMolOptUclSZVBasis():
	return basObj.CP2KBasisObjStandard( element="O", basis="SZV-MOLOPT-GTH-q6",
	                             potential="GTH-PBE-q6", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT" )

#13 basis functions
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("O-GTH-PBE-q6-MolOpt-DZVP".lower())
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("O-GTH-PBE-q6-MolOpt-ssppd".lower())
def _getOxygenMolOptUclDZPVBasis():
	return basObj.CP2KBasisObjStandard( element="O", basis="DZVP-MOLOPT-GTH-q6",
	                             potential="GTH-PBE-q6", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("O-GTH-BLYP-q6-MolOpt-DZVP".lower())
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("O-GTH-BLYP-q6-MolOpt-ssppd".lower())
def _():
	return basObj.CP2KBasisObjStandard( element="O", basis="DZVP-MOLOPT-GTH-q6",
                                 potential="GTH-BLYP-q6", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT" )

#17 basis functions
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("O-GTH-PBE-q6-MolOpt-ssspppd".lower())
def _getOxygenMolOptUclTZVPBasis():
	return basObj.CP2KBasisObjStandard( element="O", basis="TZVP-MOLOPT-GTH-q6",
	                             potential="GTH-PBE-q6", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT" )

#22 basis functions
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("O-GTH-PBE-q6-MolOpt-TZV2P".lower())
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("O-GTH-PBE-q6-MolOpt-ssspppdd".lower())
def _getOxygenMolOptUclTZV2PBasis():
	return basObj.CP2KBasisObjStandard( element="O", basis="TZV2P-MOLOPT-GTH-q6",
	                             potential="GTH-PBE-q6", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("O-GTH-BLYP-q6-MolOpt-TZV2P".lower())
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("O-GTH-BLYP-q6-MolOpt-ssspppdd".lower())
def _():
    return basObj.CP2KBasisObjStandard( element="O", basis="TZV2P-MOLOPT-GTH-q6",
                                 potential="GTH-BLYP-q6", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT" )


#29 basis functions
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("O-GTH-PBE-q6-MolOpt-ssspppddf".lower())
def _getOxygenMolOptUclTZV2P_plusFOribtal_Basis():
	return basObj.CP2KBasisObjStandard( element="O", basis="TZV2PX-MOLOPT-GTH-q6",
	                             potential="GTH-PBE-q6", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT" )

