#!/usr/bin/env bash
cp config_vars.py gen_basis_helpers/shared
python3 setup.py install --user --record files.txt
