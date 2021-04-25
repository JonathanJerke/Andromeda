#!/bin/bash

##compile first
cd Object
export LDFLAGS=-L`pwd`
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`
cd ..

cd andromedaPy
export PYTHONPATH=$PYTHONPATH:`pwd`
python3 setup.py build_ext -i
cd ..
