#!/bin/bash

echo "LAUNCH = `pwd` in PATH"
export LAUNCH=`pwd`
export PATH=$LAUNCH/csh:$PATH
cd csh
chmod +x *
