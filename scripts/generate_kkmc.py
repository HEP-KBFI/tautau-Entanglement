#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##########################################################################
# basf2 (Belle II Analysis Software Framework)                           #
# Author: The Belle II Collaboration                                     #
#                                                                        #
# See git log for contributors and copyright holders.                    #
# This file is licensed under LGPL-3.0, see LICENSE.md.                  #
##########################################################################

# Original files:
# https://github.com/belle2/basf2/blob/eee37dce021b29d78d3d266db0f3cc2c7863d8d4/generators/examples/KKGenGenerationOnly.py ("default")
# https://github.com/belle2/basf2/blob/eee37dce021b29d78d3d266db0f3cc2c7863d8d4/generators/examples/KKGenGenerationOnly_tauola_orig.py ("orig")
# https://github.com/belle2/basf2/blob/eee37dce021b29d78d3d266db0f3cc2c7863d8d4/generators/examples/KKGenGenerationOnly_tauola_bbb.py ("bbb")

# Usage:
#

import sys
import basf2 as b2
from beamparameters import add_beamparameters

# cmd line arguments
out_filename = sys.argv[1]
seed         = int(sys.argv[2])
nof_events   = int(sys.argv[3])
mode         = sys.argv[4]

assert(mode in [ 'default', 'orig', 'bbb' ])

b2.set_random_seed(seed)
b2.set_log_level(b2.LogLevel.INFO)

# main path
main = b2.create_path()

# event info setter
main.add_module("EventInfoSetter", expList=0, runList=1, evtNumList=nof_events)

# beam parameters
beamparameters = add_beamparameters(main, "Y4S")

if mode == 'default':
  tauInputFile = 'data/generators/kkmc/tau.input.dat'
  taudecaytableFile = 'data/generators/kkmc/tau_decaytable.dat'
elif mode == 'orig':
  tauInputFile = 'data/generators/kkmc/tauola_orig.input.dat'
  taudecaytableFile = ''
elif mode == 'bbb':
  tauInputFile = 'data/generators/kkmc/tauola_bbb.input.dat'
  taudecaytableFile = ''
else:
  assert(False)

# to run the framework the used modules need to be registered
kkgeninput = b2.register_module('KKGenInput')
kkgeninput.param('tauinputFile', b2.find_file(tauInputFile))
kkgeninput.param('KKdefaultFile', b2.find_file('data/generators/kkmc/KK2f_defaults.dat'))
kkgeninput.param('taudecaytableFile', b2.find_file(taudecaytableFile) if taudecaytableFile else '')
kkgeninput.param('kkmcoutputfilename', f'{out_filename}.txt')

# run
main.add_module("Progress")
main.add_module(kkgeninput)
# main.add_module("RootOutput", outputFileName=f"{out_filename}.root") #NB! not a flat tree
main.add_module("HepMCOutput", OutputFilename=f'{out_filename}.hepmc', StoreVirtualParticles = True)
# main.add_module("PrintTauTauMCParticles", logLevel=LogLevel.INFO, onlyPrimaries=False)
# main.add_module("PrintMCParticles", logLevel=LogLevel.INFO, onlyPrimaries=False)

# generate events
b2.process(main)

# show call statistics
print(b2.statistics)
