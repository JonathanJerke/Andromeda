#!/bin/bash
#source this file

# apt-get update
# apt-get install python3.6
# pip3 install cython



#run in Andromeda directory for best effect

##compile first
cd Object
make olympics
cp olympics andromeda

make libandromeda.so
export LDFLAGS=-L`pwd`
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`
cd ..

cd andromedaPy
export PYTHONPATH=$PYTHONPATH:`pwd`
python3 setup.py build_ext -i
cd ..

##run python3 
