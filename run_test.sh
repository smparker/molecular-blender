#!/bin/bash

# Check raw location
if [ -x "$1" ]; then
  BLENDER="$1"
  shift
elif [ -x "$(which blender)" ]; then
  BLENDER=blender
elif [ -x "/Applications/blender.app/Contents/MacOS/blender" ]; then
  BLENDER=/Applications/blender.app/Contents/MacOS/blender
else
  echo "Failed to find a suitable Blender! Add it to your PATH and try again!"
  exit 1
fi

# jump to location of this script
cd $(dirname $0)

# pass on arguments to unittest
args="$@"
${BLENDER} --background --python-exit-code 1 --python test.py -- $args
