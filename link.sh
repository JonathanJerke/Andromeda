#!/bin/bash
export LDFLAGS=-L`pwd`
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`
export PYTHONPATH=$PYTHONPATH:`pwd`
export PATH=$PATH:`pwd`/bash:`pwd`/Object:`pwd`/python3
