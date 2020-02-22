


import gen_basis_helpers.cp2k.cp2k_basis_obj as stdObjs




def _createStandardBasisObjForUnitTests():
	return stdObjs.CP2KBasisObjStandard( element="should_be_overwritten", basis="test_basis",
	                                     potential="test_potential", basisFile="test_basFile", potFile="test_potFile" )



defaultBasisStrsToCreators = {"utests-basis":_createStandardBasisObjForUnitTests}

