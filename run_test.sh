#!/bin/bash

# MacOSX location
blender=/Applications/blender.app/Contents/MacOS/blender

if [ ! -f $blenderdir ]; then
  echo "Blender not found!"
  exit 1
fi

# jump to location of this script
cd $(dirname $0)

# pass on arguments to unittest
args="$@"
${blender} --background --python test.py -- $args
