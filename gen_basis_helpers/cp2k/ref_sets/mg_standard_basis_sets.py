
""" Contains names of basis sets for registering """

from .. import basis_register as basRegister
from .. import cp2k_basis_obj as basObj


@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("GTH-PBE-q2-MolOpt-UCL-SZV".lower())
def _getMgMolOptUclSZVBasis():
	return basObj.CP2KBasisObjStandard( element="will_be_overwritten", basis="SZV-MOLOPT-SR-GTH-q2",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT_UCL" )


#ssp i think
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("GTH-PBE-q2-MolOpt-UCL-ssp".lower())
def _getMgMolOptUclDZVP_basis():
	return basObj.CP2KBasisObjStandard( element="will_be_overwritten", basis="DZVP-MOLOPT-SR-GTH-q2",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT_UCL" )

#sssp i think
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("GTH-PBE-q2-MolOpt-UCL-sssp".lower())
def _getMgMolOptUclTZVP_basis():
	return basObj.CP2KBasisObjStandard( element="will_be_overwritten", basis="TZVP-MOLOPT-SR-GTH-q2",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT_UCL" )


#ssspp i think
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("GTH-PBE-q2-MolOpt-UCL-ssspp".lower())
def _getMgMolOptUclTZVP_basis():
	return basObj.CP2KBasisObjStandard( element="will_be_overwritten", basis="TZV2P-MOLOPT-SR-GTH-q2",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT_UCL" )




