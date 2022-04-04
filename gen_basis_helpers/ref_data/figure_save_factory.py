
import os

from . import figure_savers_base as baseObjs 
from ..shared import config_vars as configVars


def createFigureSaver_twoCentreTightBindingPaper():
	baseFolder = configVars.TWO_CENTRE_PAPER_2019_FIG_SAVE_FOLDER
	return baseObjs.StandardFigureSaver.fromDefaultPlusKwargs(baseFolder=baseFolder)

def createFigureSaver_threeCentreTightBindingPaper():
	baseFolder = configVars.THREE_CENTRE_PAPER_2019_FIG_SAVE_FOLDER
	return baseObjs.StandardFigureSaver.fromDefaultPlusKwargs(baseFolder=baseFolder)

def createFigureSaver_threeCentreTightBindingBostomTalk():
	outSaver = createFigureSaver_threeCentreTightBindingPaper()
	outSaver.baseFolder = os.path.join( outSaver.baseFolder, "mrs_talk" )
	return outSaver

def createFigureSaver_cp2kBasisPaper2020():
	baseFolder = configVars.CP2K_BASIS_PAPER_2020_FIG_SAVE_FOLDER
	return baseObjs.StandardFigureSaver.fromDefaultPlusKwargs(baseFolder=baseFolder)

def createFigureSaver_mgWaterMdPaper2021():
	baseFolder = configVars.MG_WATER_MD_PAPER_2021_FIG_SAVE_FOLDER
	return baseObjs.StandardFigureSaver.fromDefaultPlusKwargs(baseFolder=baseFolder)

def createFigureSaver_hydroxylMgWaterMdPaper2022():
	baseFolder = configVars.HYDROXYL_MG_WATER_MD_PAPER_2022_FIG_SAVE_FOLDER
	return baseObjs.StandardFigureSaver.fromDefaultPlusKwargs(baseFolder=baseFolder)
