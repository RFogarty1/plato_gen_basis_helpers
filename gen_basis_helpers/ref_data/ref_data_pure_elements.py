
""" Combines data for various pure-elements into a small number of objects. Effectively the interface module """

from types import SimpleNamespace
from . import ref_elemental_objs as refEleObjs
from . import ref_data_mg as refMg
from . import ref_data_zr as refZr

elementalData = SimpleNamespace(Mg=refMg.createMgReferenceDataObj(),
                                Zr=refZr.createZrReferenceDataObj())


shellAngMomMaps = refEleObjs.ShellAngMomMapper({"Mg": refMg.createMgAngMomShellIndices(),
                                                "Zr": refZr.createZrAngMomShellIndices()})




