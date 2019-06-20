
import itertools as it
import os
import sys
sys.path.append('/media/ssd1/rf614/Work/usr_scripts/coding/Plato_Analysis_Lib_Functions')

import parse_ham_files as parseHam
import parse_plato_out_files as parseOut

def getVolVsOnSiteDiagTerms(inpFileList):
	''' Note the diags have the struct [atomIdx][row or col(obviously same number)] '''
	allVolumes = list()
	allDiagTerms = list()
	for fPath in inpFileList:
		baseFilePath = os.path.splitext(fPath)[0]
		outFilePath, hamFilePath = baseFilePath + ".out", baseFilePath + ".ham"
		allVolumes.append( parseOut.parsePlatoOutFile( outFilePath )["unitCell"].getVolume() )
		allDiagTerms.append( parseHam.getOnSiteDiagHamilTermsFromHamFile(hamFilePath) )

	outStruct = list()
	for vol,diags in it.zip_longest(allVolumes, allDiagTerms):
		outStruct.append( [vol,diags] )

	return outStruct
