
""" Combines data for various pure-elements into a single object. Effectively the interface object """

from types import SimpleNamespace
import ref_data_mg as refMg
import ref_data_zr as refZr

elementalData = SimpleNamespace(Mg=refMg.createMgReferenceDataObj(),
                                Zr=refZr.createZrReferenceDataObj())



