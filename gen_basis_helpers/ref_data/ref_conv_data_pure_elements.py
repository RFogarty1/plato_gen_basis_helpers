
""" Stub right now, will be the interface eventually though """

from types import SimpleNamespace
from . import ref_conv_data_mg as refMg
from . import ref_conv_data_zr as refZr

elementalConvData = SimpleNamespace(Mg=refMg.createConvObj(),
                                    Zr=refZr.createConvObj())

