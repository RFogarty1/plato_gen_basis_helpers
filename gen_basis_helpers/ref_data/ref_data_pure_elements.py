
""" Combines data for various pure-elements into a small number of objects. Effectively the interface module """

from types import SimpleNamespace
import os
import plato_pylib.plato.plato_paths as platoPaths
from . import ref_elemental_objs as refEleObjs
from . import ref_data_mg as refMg
from . import ref_data_zr as refZr

elementalData = SimpleNamespace(Mg=refMg.createMgReferenceDataObj(),
                                Zr=refZr.createZrReferenceDataObj())


shellAngMomMaps = refEleObjs.ShellAngMomMapper({"Mg": refMg.createMgAngMomShellIndices(),
                                                "Zr": refZr.createZrAngMomShellIndices()})


#Swap all datasets to use the McWeda models
elementalDataMcWeda = SimpleNamespace(Mg=refMg.createMgReferenceDataObj(),
                                Zr=refZr.createZrReferenceDataObj())


elementalDataMcWeda.Mg.modelFiles.tb1Model = platoPaths.getAbsolutePathForPlatoTightBindingDataSet( os.path.join("Mg_models","three_body_2019") )
elementalDataMcWeda.Mg.modelFiles.dft2Model = platoPaths.getAbsolutePathForPlatoTightBindingDataSet( os.path.join("Mg_models","three_body_2019") )
elementalDataMcWeda.Mg.modelFiles.dftModel = platoPaths.getAbsolutePathForPlatoTightBindingDataSet( os.path.join("Mg_models","three_body_2019"), dtype="dft")

elementalDataMcWeda.Zr.modelFiles.tb1Model = platoPaths.getAbsolutePathForPlatoTightBindingDataSet( os.path.join("Zr_models","three_body_2019") )
elementalDataMcWeda.Zr.modelFiles.dft2Model = platoPaths.getAbsolutePathForPlatoTightBindingDataSet( os.path.join("Zr_models","three_body_2019") )
elementalDataMcWeda.Zr.modelFiles.dftModel = platoPaths.getAbsolutePathForPlatoTightBindingDataSet( os.path.join("Zr_models","three_body_2019"), dtype="dft")


#Swap relevant data-sets to use LDA
elementalDataLDA = SimpleNamespace(Mg=refMg.createMgReferenceDataObj(),
                                   Zr=refZr.createZrReferenceDataObj())

elementalDataLDA.Mg.modelFiles.tb1Model  = platoPaths.getAbsolutePathForPlatoTightBindingDataSet( os.path.join("Mg_models", "two_body_2019", "lda_version") )
elementalDataLDA.Mg.modelFiles.dft2Model = platoPaths.getAbsolutePathForPlatoTightBindingDataSet( os.path.join("Mg_models", "two_body_2019", "lda_version") )

elementalDataLDA.Zr.modelFiles.tb1Model  = platoPaths.getAbsolutePathForPlatoTightBindingDataSet( os.path.join("Zr_models", "two_body_2019", "lda_version") )
elementalDataLDA.Zr.modelFiles.dft2Model = platoPaths.getAbsolutePathForPlatoTightBindingDataSet( os.path.join("Zr_models", "two_body_2019", "lda_version") )

