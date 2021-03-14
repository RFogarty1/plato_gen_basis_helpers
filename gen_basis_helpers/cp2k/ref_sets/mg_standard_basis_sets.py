
""" Contains names of basis sets for registering """

from .. import basis_register as basRegister
from .. import cp2k_basis_obj as basObj

#MOLOPT BASIS SETS
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("GTH-PBE-q2-MolOpt-UCL-SZV".lower())
def _getMgMolOptUclSZVBasis():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="SZV-MOLOPT-SR-GTH-q2",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT_UCL" )


#ssp i think
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q2-MolOpt-UCL-DZVP".lower())
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("GTH-PBE-q2-MolOpt-UCL-ssp".lower())
def _getMgMolOptUclDZVP_basis():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="DZVP-MOLOPT-SR-GTH-q2",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT_UCL" )

#sssp i think
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q2-MolOpt-UCL-TZVP".lower())
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("GTH-PBE-q2-MolOpt-UCL-sssp".lower())
def _getMgMolOptUclTZVP_basis():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="TZVP-MOLOPT-SR-GTH-q2",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT_UCL" )


#ssspp i think
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q2-MolOpt-UCL-TZV2P".lower())
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("GTH-PBE-q2-MolOpt-UCL-ssspp".lower())
def _getMgMolOptUclTZVP_basis():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="TZV2P-MOLOPT-SR-GTH-q2",
	                             potential="GTH-PBE-q2", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT_UCL" )




#q10 Pseudopotential molopt basis sets
@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q10-MolOpt-DZVP".lower())
def _unused():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="DZVP-MOLOPT-SR-GTH-q10",
	                                    potential="GTH-PBE-q10", potFile="GTH_POTENTIALS", basisFile="BASIS_MOLOPT" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q10-DZVP".lower())
def _unused():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="DZVP-GTH-q10",
	                                    potential="GTH-PBE-q10", potFile="GTH_POTENTIALS", basisFile="GTH_BASIS_SETS" )


@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q10-TZVP".lower())
def _unused():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="TZVP-GTH-q10",
	                                    potential="GTH-PBE-q10", potFile="GTH_POTENTIALS", basisFile="GTH_BASIS_SETS" )

@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Mg-GTH-PBE-q10-TZV2P".lower())
def _unused():
	return basObj.CP2KBasisObjStandard( element="Mg", basis="TZV2P-GTH-q10",
	                                    potential="GTH-PBE-q10", potFile="GTH_POTENTIALS", basisFile="GTH_BASIS_SETS" )








