#!/usr/bin/env cmsRun

# Adapted from: IOMC/Input/test/hepmc2gen.py
# Usage (maxEvents not required):
# ./hepmc2gen.py inputFiles=/some/input/file.hepmc outputFile=/some/output/file.root maxEvents=100

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.parseArguments()

import FWCore.ParameterSet.Config as cms

import subprocess
import threading
import time
import re
import os
import stat

FILE_PFX = 'file:'

def file_pfx(path):
  return FILE_PFX + path
re_file_pfx = re.compile('^{}'.format(FILE_PFX))

inputFiles = []
for inputFile_ in options.inputFiles:
  inputFile = os.path.abspath(re_file_pfx.sub('', inputFile_))
  if inputFile.endswith('.hepmc'):
    inputFiles.append(file_pfx(inputFile))
  elif inputFile.endswith(('.tgz', '.tar.gz')):

    def call_exe(cmd):
      print("Executing command: '{}'".format(cmd))
      status = subprocess.call(cmd, shell = True)
      print("End command: status = {}".format(status))
      return

    def create_fifo(fn):
      cmd = "mkfifo {}".format(fn)
      subprocess.call(cmd, shell = True)

    # create a HepMC file in current working dir
    tmpdir = os.getcwd()
    tmpfn = os.path.join(tmpdir, re.sub('(tar.gz|tgz)$', 'hepmc', os.path.basename(inputFile)))
    tmpfn_exists = os.path.exists(tmpfn)

    if not tmpfn_exists:
      create_fifo(tmpfn)
    else:
      if not stat.S_ISFIFO(os.stat(tmpfn).st_mode):
        print("Removing {}".format(tmpfn))
        os.remove(tmpfn)
      create_fifo(tmpfn)

    # unpack the tarball
    exe = "cat {} | gunzip -c > {} &".format(inputFile, tmpfn)
    t = threading.Thread(target = call_exe, args = ([ exe ]))
    t.start()

    # sleep a little
    sleep_duration = 1
    print("Sleeping {}s to allow start of pipes".format(sleep_duration))
    time.sleep(sleep_duration)

    print("Reading from {} instead of {}".format(inputFile, tmpfn))
    inputFiles.append(file_pfx(tmpfn))

outputFile = options.outputFile
if not outputFile.startswith(FILE_PFX):
  outputFile = file_pfx(outputFile)

process = cms.Process("GEN")

process.source = cms.Source("MCFileSource",
  fileNames = cms.untracked.vstring(inputFiles),
  firstLuminosityBlockForEachRun = cms.untracked.VLuminosityBlockID(), # why is this required?!
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load('Configuration.EventContent.EventContent_cff')
process.AODSIMoutput = cms.OutputModule("PoolOutputModule",
  compressionAlgorithm = cms.untracked.string('LZMA'),
  compressionLevel = cms.untracked.int32(4),
  dataset = cms.untracked.PSet(
    dataTier = cms.untracked.string('AODSIM'),
    filterName = cms.untracked.string('')
  ),
  eventAutoFlushCompressedSize = cms.untracked.int32(31457280),
  fileName = cms.untracked.string(outputFile),
  outputCommands = process.AODSIMEventContent.outputCommands
)

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.genParticles.src = cms.InputTag("source", "generator")


process.p = cms.Path(process.genParticles)
process.outpath = cms.EndPath(process.AODSIMoutput)

process.AODSIMoutput.outputCommands.extend([
  'keep GenRunInfoProduct_*_*_*',
  'keep GenLumiInfoProduct_*_*_*',
  'keep GenEventInfoProduct_*_*_*',
])
