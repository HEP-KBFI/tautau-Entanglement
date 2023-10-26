#!/bin/bash

# Example usage:
# run_cluster_job.sh /some/tmp/dir output.root cfg.py out.log in-1.root in-2.root ... in-N.root

if [ "$#" -lt 4 ]; then
  echo "Need to provide at least four arguments: ${0} [tmp-dir] [output-file] [cmd] [config-file]";
  exit 1;
fi

TMP_DIRNAME=$1
OUT_FILENAME=$2
CMD=$3
CFG_FILENAME=$4

if ! command -v "$CMD" &> /dev/null; then
  echo "Not a valid command: $CMD";
  exit 2;
fi

if [ ! -f $CFG_FILENAME ]; then
  echo "No such file: $CFG_FILENAME";
  exit 3;
fi

OUT_DIRNAME=$(dirname $OUT_FILENAME)
if [ ! -d $OUT_DIRNAME ]; then
  echo "No such directory: $OUT_DIRNAME";
  exit 4;
fi

echo "Temporary directory is: $TMP_DIRNAME"
echo "Hostname: $HOSTNAME"
echo "Time: `date`"

mkdir -pv $TMP_DIRNAME
if [ ! -d $TMP_DIRNAME ]; then
  echo "Cannot create directory: $TMP_DIRNAME"
  exit 5;
fi
PREV_DIRNAME=$PWD;
cd $TMP_DIRNAME

/usr/bin/time --verbose $CMD $CFG_FILENAME
EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
  OUT_BASENAME=$(basename $OUT_FILENAME)
  if [ -f $OUT_BASENAME ]; then
    cp -v $OUT_BASENAME $OUT_DIRNAME
  else
    echo "Missing file: $OUT_BASENAME?"
  fi

  for ext in pdf png; do
    if ls *.$ext 1> /dev/null 2>&1; then
      cp -v *.$ext $OUT_DIRNAME
    fi
  done

else
  echo "Job finished with non-zero exit code ($EXIT_CODE) -> not copying anything"
fi

cd -
if [ "$PREV_DIRNAME" != "$TMP_DIRNAME" ]; then
  # remove only if we actually changed dirs
  rm -rfv $TMP_DIRNAME
fi

exit $EXIT_CODE
