
from .. import basis_register as basRegister
from .. import cp2k_basis_obj as basObj


@basRegister.decoRegisterCP2KBasisCreatorToMethodStr("Si-GTH-PBE-q4-DZVP")
def _unused():
	return basObj.CP2KBasisObjStandard( element="Si", basis="DZVP-GTH-q4",
	                                    potential="GTH-PBE-q4", potFile="GTH_POTENTIALS", basisFile="GTH_BASIS_SETS" )


