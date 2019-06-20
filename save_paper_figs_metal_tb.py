
import os
import pathlib

''' Code to save figure in a sensible place for working on the tight binding paper '''

def saveFigForMetalTbPaper(figHandle, saveFolderName, saveFileName, fmt="eps"):
	basePath = "/media/ssd1/rf614/Work/Documents/papers/janas_metal_tb/figures"
	saveFolder = os.path.join(basePath,saveFolderName)
	savePath = os.path.join(basePath,saveFolderName, saveFileName + "." + fmt)
	pathlib.Path(saveFolder).mkdir(exist_ok=True, parents=True)
	figHandle.savefig( savePath,format=fmt,bbox_inches='tight')


