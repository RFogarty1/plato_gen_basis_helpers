
import itertools as it
import os
import sys
sys.path.append('/media/ssd1/rf614/Work/usr_scripts/coding/Plato_Analysis_Lib_Functions')

import plato_pylib.plato.parse_matrices_out_files as parseMatrix
import plato_pylib.plato.parse_plato_out_files as parseOut

def getVolVsOnSiteDiagTerms(inpFileList, volPerAtom=False):
	''' Note the diags have the struct [atomIdx][row or col(obviously same number)] '''
	allVolumes = list()
	allDiagTerms = list()
	for fPath in inpFileList:
		baseFilePath = os.path.splitext(fPath)[0]
		outFilePath, hamFilePath = baseFilePath + ".out", baseFilePath + ".ham"
		currVol = _getVolFromParsedFile( parseOut.parsePlatoOutFile( outFilePath ), volPerAtom )
		allVolumes.append( currVol )
		allDiagTerms.append( parseMatrix.getOnSiteDiagHamilTermsFromHamFile(hamFilePath) )

	outStruct = list()
	for vol,diags in it.zip_longest(allVolumes, allDiagTerms):
		outStruct.append( [vol,diags] )

	return outStruct



def _getVolFromParsedFile(parsedFile, volPerAtom):
	totalVol = parsedFile["unitCell"].getVolume()
	if volPerAtom:
		return totalVol / parsedFile["numbAtoms"]
	else:
		return totalVol

