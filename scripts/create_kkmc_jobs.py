#!/usr/bin/env python3

# example usage:
# create_jobs.py default 10 10000 /local/$USER/belle_eeToTauTau/kkmc/default_10jobs_10Kevents run

from TauAnalysis.Entanglement.tools.jobTools import mkdir

import sys
import os
import re

mode               = sys.argv[1]
nof_jobs           = int(sys.argv[2])
nof_events_per_job = int(sys.argv[3])
dst_dir            = sys.argv[4]
job_type           = sys.argv[5]

submission_filename = 'submit.sh'
if len(sys.argv) > 6:
  submission_filename = sys.argv[6]
assert(not os.path.exists(submission_filename))

assert(job_type in [ 'run', 'test' ])

mkdir(dst_dir, True)

assert(dst_dir.startswith('/local/'))
log_dir = re.sub(r'^/local/', '/home/', dst_dir)
if job_type == 'run':
  mkdir(log_dir)

with open(submission_filename, 'w') as submission_file:
  submission_file.write('#!/bin/bash\n\n')
  for job_idx in range(nof_jobs):
    submission_cmd = f'sbatch --partition=io --output={log_dir}/{mode}_{job_idx}.log' if job_type == 'run' else ''
    submission_file.write(f"{submission_cmd} run_kkmc_job.sh {mode} {job_idx} {nof_events_per_job} {dst_dir}\n".lstrip())

os.chmod(submission_filename, 0o755)
print("Run: {}".format(os.path.abspath(submission_filename)))
