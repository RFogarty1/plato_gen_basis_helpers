


from .. import basis_register as basRegister
from .. import cp2k_basis_obj as basObj

#4 basis functions
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("O-GTH-PBE-q6-MolOpt-UCL-sp".lower())
def _getOxygenMolOptUclSZVBasis():
	return basObj.CP2KBasisObjStandard( element="will_be_overwritten", basis="SZV-MOLOPT-GTH-q6",
	                             potential="GTH-PBE-q6", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT_UCL" )

#13 basis functions
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("O-GTH-PBE-q6-MolOpt-UCL-ssppd".lower())
def _getOxygenMolOptUclDZPVBasis():
	return basObj.CP2KBasisObjStandard( element="will_be_overwritten", basis="DZVP-MOLOPT-GTH-q6",
	                             potential="GTH-PBE-q6", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT_UCL" )

#17 basis functions
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("O-GTH-PBE-q6-MolOpt-UCL-ssspppd".lower())
def _getOxygenMolOptUclTZVPBasis():
	return basObj.CP2KBasisObjStandard( element="will_be_overwritten", basis="TZVP-MOLOPT-GTH-q6",
	                             potential="GTH-PBE-q6", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT_UCL" )

#22 basis functions
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("O-GTH-PBE-q6-MolOpt-UCL-ssspppdd".lower())
def _getOxygenMolOptUclTZV2PBasis():
	return basObj.CP2KBasisObjStandard( element="will_be_overwritten", basis="TZV2P-MOLOPT-GTH-q6",
	                             potential="GTH-PBE-q6", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT_UCL" )

#29 basis functions
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("O-GTH-PBE-q6-MolOpt-UCL-ssspppddf".lower())
def _getOxygenMolOptUclTZV2P_plusFOribtal_Basis():
	return basObj.CP2KBasisObjStandard( element="will_be_overwritten", basis="TZV2PX-MOLOPT-GTH-q6",
	                             potential="GTH-PBE-q6", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT_UCL" )

