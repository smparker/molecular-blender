#!/bin/bash

set -e

usage() { echo "Usage: $0 [-t | -p] [-b /path/to/blender]" 1>&2; exit 1; }

run_blender() {
  py=$1
  shift
  ${BLENDER} --background --python-exit-code 1 --python $py -- $args
}

test_default="yes"
prof_default="no"

do_test="default"
do_prof="default"
blender="auto"

while getopts "tpb:" o; do
    case "${o}" in
        t)
          do_test="yes"
          ;;
        p)
          do_prof="yes"
          ;;
        b)
          blender=${OPTARG}
          ;;
        *)
          usage
          ;;
    esac
done
shift $((OPTIND-1))

# Check raw location
if [ "$blender" == "auto" ]; then
  if [ -x "$(which blender)" ]; then
    BLENDER=blender
  elif [ -x "/Applications/Blender.app/Contents/MacOS/blender" ]; then
    BLENDER=/Applications/Blender.app/Contents/MacOS/blender
  fi
else
  echo "Checking for blender at $blender"
  if [ -x "$blender" ]; then
    BLENDER="$blender"
  else
    echo "Failed to find an executable at $blender"
  fi
fi

# check that $BLENDER was set
if [ -z $BLENDER ]; then
  echo "Failed to find a suitable Blender! Add it to your PATH or provide it with the -b option!"
  exit 1
fi

# if profiling turned on, change default behavior of test
if [ "$do_prof" == "yes" ]; then
  test_default="no"
fi

# apply defaults
[ "$do_test" == "default" ] && do_test="$test_default"
[ "$do_prof" == "default" ] && do_prof="$prof_default"

# jump to location of this script
cd $(dirname $0)

# pass on arguments to unittest
args="$@"

if [ "$do_test" == "yes" ]; then
  run_blender test.py $args
fi

if [ "$do_prof" == "yes" ]; then
  run_blender marching_profile.py $args
fi
