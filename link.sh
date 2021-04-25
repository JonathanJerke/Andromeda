#!/bin/bash

##MEANT FOR new bash shell without recompiling binary
##RELINKS python so that installations issues do not arise


#run in Andromeda directory


export PATH=$PATH:`pwd`/bash:`pwd`/Object:`pwd`/python3

cd Object
export LDFLAGS=-L`pwd`
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`
cd ..

cd andromedaPy
export PYTHONPATH=$PYTHONPATH:`pwd`
python3 setup.py build_ext -i
cd ..
