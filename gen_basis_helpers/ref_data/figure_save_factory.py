
from . import figure_savers_base as baseObjs 



def createFigureSaver_twoCentreTightBindingPaper():
	baseFolder = "/media/ssd1/rf614/Work/Documents/jobs/Corrosion_Work/testing_models/paper_figures/two_centre_tight_binding"
	return baseObjs.StandardFigureSaver.fromDefaultPlusKwargs(baseFolder=baseFolder)


