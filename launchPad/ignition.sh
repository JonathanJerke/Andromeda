#!/bin/bash

echo "LAUNCH = `pwd` in PATH"
export LAUNCH=`pwd`
export PATH=$LAUNCH/../bash:$LAUNCH/../Object:$PATH
