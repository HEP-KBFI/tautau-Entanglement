
import os
import re
import subprocess

def getInputFileNames(inputFilePath, inputFile_regex = r"[a-zA-Z0-9-_]+.root"):
    inputFile_matcher = re.compile(inputFile_regex)
    inputFileNames = []
    files = os.listdir(inputFilePath)
    for file in files:
        if os.path.isdir(os.path.join(inputFilePath, file)):
            inputFileNames.extend(getInputFileNames(os.path.join(inputFilePath, file)))
        else:
            # check if name of inputFile matches regular expression
            if inputFile_matcher.match(file):
                inputFileNames.append("file:%s" % os.path.join(inputFilePath, file))
    return inputFileNames

def build_Makefile(makeFileName, jobOptions):
  print("Building Makefile = '%s'" % makeFileName)

  lines_Makefile = []
  lines_Makefile.append(".DEFAULT_GOAL := all")
  lines_Makefile.append("")
  lines_Makefile.append("SHELL := /bin/bash")
  lines_Makefile.append("")
  lines_Makefile.append("all: %s" % " ".join([ job['target'] for job in jobOptions ]))
  lines_Makefile.append("")
  for job in jobOptions:
    lines_Makefile.append("%s: %s" % (job['target'], " ".join(job['dependencies'])))
    for command in job['commands']:
      lines_Makefile.append("\t%s" % command)
    lines_Makefile.append("")
  lines_Makefile.append("clean:")
  outputFileNames = []
  for job in jobOptions:
    outputFileNames.extend(job['outputFileNames'])
  for outputFileName in outputFileNames:
    lines_Makefile.append("\trm -f %s" % outputFileName)
  lines_Makefile.append("")

  makeFile = open(makeFileName, "w") 
  for line in lines_Makefile:
    makeFile.write("%s\n" % line)
  makeFile.close()
 
