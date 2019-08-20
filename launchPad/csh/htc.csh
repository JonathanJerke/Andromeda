#!/bin/bash

module load launcher
export LAUNCHER_JOB_FILE=launcher
export LAUNCHER_NHOSTS=1
export LAUNCHER_WORKDIR=`pwd`
export LAUNCHER_PPN=12

$LAUNCHER_DIR/paramrun
