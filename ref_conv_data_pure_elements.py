
""" Stub right now, will be the interface eventually though """

from types import SimpleNamespace
import ref_conv_data_mg as refMg
import ref_conv_data_zr as refZr

elementalConvData = SimpleNamespace(Mg=refMg.createConvObj(),
                                    Zr=refZr.createConvObj())

