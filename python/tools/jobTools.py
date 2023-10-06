
import os
import re
import sys

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
 
def build_sbatchSubmission(sbatchSubmissionFileName, jobOptions, job_type):
  OUTPUT_PATTERN = '^/scratch/persistent/'

  with open(sbatchSubmissionFileName, 'w') as sbatchSubmissionFile:
    sbatchSubmissionFile.write('#!/bin/bash\n\n')
    for job_key, job in jobOptions.items():

      # cannot copy to /scratch/persistent from comp nodes -> copy to local instead and create a symlink
      output_dir = jobOptions[job_key]['outputFilePath']
      assert(re.match(OUTPUT_PATTERN, output_dir))
      output_dir_old = output_dir
      output_dir = re.sub(OUTPUT_PATTERN, "/local/", output_dir)
      if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
      if os.path.isdir(output_dir_old) and not os.path.islink(output_dir_old):
        os.rmdir(output_dir_old)
        os.symlink(output_dir, output_dir_old)

      sbatchSubmissionFile.write(
        "sbatch --partition=main --output={log} run_cluster_job.sh {tmp} {out} {cfg} {ins}\n".format(
          log = jobOptions[job_key]['logFileName'],
          tmp = os.path.join('/scratch', 'local', os.getenv('USER'), job_type, job_key),
          out = os.path.join(output_dir, jobOptions[job_key]['outputFileName']),
          cfg = jobOptions[job_key]['cfgFileName'],
          ins = ' '.join(jobOptions[job_key]['inputFileNames']),
        )
      )
    sbatchSubmissionFile.write('\n')

def query_yes_no(question, default = "yes"):
  """Prompts user yes/no

  Args:
    question: question to ask from the user
    default: default option to use; acceptable values: "yes", "no" or None
  """
  default = default.lower()
  valid = { "yes": True, "y": True, "ye": True, "no": False, "n": False }
  if default is None:    prompt = " [y/n] "
  elif default == "yes": prompt = " [Y/n] "
  elif default == "no":  prompt = " [y/N] "
  else:
    raise ValueError("Invalid default answer: '%s'" % default)

  while True:
    sys.stdout.write(question + prompt)
    choice = input().lower()
    if default is not None and choice == "": return valid[default]
    elif choice in valid:                    return valid[choice]
    else:
      sys.stdout.write("Please respond with 'yes' or 'no' (or 'y' or 'n').\n")
