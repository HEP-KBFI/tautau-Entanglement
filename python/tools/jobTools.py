
import os
import re
import sys
import subprocess
import jinja2
import argparse

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
  lines_Makefile.append("all: \\")
  for job_idx, job in enumerate(jobOptions):
    target = f"  {job['target']}"
    if job_idx < len(jobOptions) - 1:
      target += " \\"
    lines_Makefile.append(target)
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
 
def build_sbatchSubmission(sbatchSubmissionFileName, jobOptions, job_type, version):
  OUTPUT_PATTERN = '^/scratch/persistent/'

  with open(sbatchSubmissionFileName, 'w') as sbatchSubmissionFile:
    sbatchSubmissionFile.write('#!/bin/bash\n\n')
    for job_key, job in jobOptions.items():

      # cannot copy to /scratch/persistent from comp nodes -> copy to local instead and create a symlink
      output_dir = jobOptions[job_key]['outputFilePath']
      if re.match(OUTPUT_PATTERN, output_dir):
        output_dir_old = output_dir
        output_dir = re.sub(OUTPUT_PATTERN, "/local/", output_dir)
        if not os.path.isdir(output_dir):
          os.makedirs(output_dir)
        if os.path.isdir(output_dir_old) and not os.path.islink(output_dir_old):
          os.rmdir(output_dir_old)
          os.symlink(output_dir, output_dir_old)

      sbatchSubmissionFile.write(
        "sbatch --partition=main --output={log} {opt} run_cluster_job.sh {tmp} {out} {cmd} {cfg}\n".format(
          cmd = jobOptions[job_key]['cmd'],
          log = jobOptions[job_key]['logFileName'],
          tmp = os.path.join('/scratch', 'local', os.getenv('USER'), version, job_type, job_key),
          out = os.path.join(output_dir, jobOptions[job_key]['outputFileName']),
          cfg = jobOptions[job_key]['cfgFileName'],
          opt = jobOptions[job_key]['options'] if 'options' in jobOptions[job_key] else '',
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

def run_command_os(cmd, verbose = False):
  if verbose:
    print(f"Executing: {cmd}")
  os.system(cmd)

def run_command_subprocess(cmd, verbose = False):
  if verbose:
    print(f"Executing: {cmd}")
  try:
    process = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
    stdout, stderr = process.communicate()
    exit_code = process.returncode
    if exit_code != 0 or stderr:
      print(f"Got an exit code of {exit_code} with the following error message: '{stderr.decode()}'")
  except Exception as e:
    print(e)
    raise

def read_contents(filename):
  with open(filename, 'r') as fileptr:
    lines = fileptr.readlines()
  contents = '\n'.join([ line.rstrip('\n') for line in lines if not (line.startswith('#') or not line.strip()) ])
  return contents

def save_cmd(filename):
  cmd = ' '.join(sys.argv)
  if os.path.isfile(filename):
    file_contents = ' '.join(read_contents(filename).split())
    msg = f"It seems that you've run this workflow before with command '{file_contents}'.\n"
    if file_contents != cmd:
      msg += "Also, the command looks different from what you're currently trying to run.\n"
    proceed = query_yes_no(f"{msg}Do you still want to proceed?", default = 'no')
    if not proceed:
      return False
  with open(filename, 'w') as fileptr:
    fileptr.write('{}\n'.format(cmd))
  return True

def build_cfg(template_contents, output_filename, args):
  print(f"Building configFile = '{output_filename}'")
  with open(output_filename, 'w') as output_file:
    output_file.write(jinja2.Template(template_contents).render(**args))

def mkdir(dirname, verbose = False):
  if not os.path.isdir(dirname):
    if verbose:
      print(f"Creating directory: {dirname}")
    os.makedirs(dirname)

def positive_int_type(value):
  try:
    value_int = int(value)
  except ValueError:
    raise argparse.ArgumentTypeError('Not an integer: %s' % value)
  if value_int <= 0:
    raise argparse.ArgumentTypeError('Must be a positive integer: %d' % value_int)
  return value_int
