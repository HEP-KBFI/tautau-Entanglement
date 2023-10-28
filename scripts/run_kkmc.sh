#!/bin/bash

# cmd line arguments
OUTFN=$1
SEED=$2
NOF_EVENTS=$3
MODE=$4

# set up b2
source /belle2/tools/b2setup $(ls /belle2/releases)

# generate the events
$CMSSW_BASE/bin/$SCRAM_ARCH/generate_kkmc.py $OUTFN $SEED $NOF_EVENTS $MODE
