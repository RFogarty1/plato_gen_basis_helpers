
import os

from ..shared import config_vars as cfgVars
from ..shared import data_plot_base as dPlotHelp

def getPlottersFromRecordStandard_mgWaterRepo(record, baseDbDir=cfgVars.MG_WATER_REPO_PATH):
	""" Gets iter of plotters from a record in the plotters collection for mg_water_work. 
	
	Args:
		record: (dict) The record containing the data. Most importantly it should have a "plotter_path_ext" key, which holds info on where the data file is
		baseDbDir: (str) Base directory for the repo.
			 
	Returns
		outPlotters: (iter of DataPlotterStandard)
 
	"""
	outPath = os.path.join(baseDbDir, record["plotter_path_ext"])
	return dPlotHelp.readStandardDataPlottersFromJson(outPath)

