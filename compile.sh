#!/bin/bash
##source compile.sh
##if you run this a second time just run 
##bash compile.sh

#run in Andromeda directory 

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

cd bash
chmod +x *
cd ..

cd python3
chmod +x *
cd ..


##run python3 
