#!/bin/bash

#run in Andromeda directory for best effect

##compile first
cd Object
make olympics
cp olympics andromeda

make libandromeda.so
export LDFLAGS=-L`pwd`
export LD_LIBRARY_PATH=`pwd`
cd ..

cd andromedaPy
python3 setup.py build_ext -i
cd ..

##run python3 
##import andromedaPy as apy
##code still under development
