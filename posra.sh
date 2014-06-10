#!/bin/sh

# @DESC Wrapper script (bash only!) to clean up an image before passing it to osra.

if [ $# -eq 0 ]; then echo "Usage: $0 [options]." >&2; exit 1; fi

image=${!#} # image assumed the last arg, e.g. MMA.gif  
prefix=${image%.*}  # prefix
suffix=${image##*.} # suffix

imagea0="${prefix}a0.${suffix}" # image w/o alpha channel

/bin/convert $image -background white -alpha remove -alpha off $imagea0 # ImageMagick convert

set -- "${@:1:2}" $imagea0 # replace the last arg

smiles=$(/Projects/Applications/GitHub/UCSC-IBM-POSRA-CORE/src/osra.exe $*)
RC=$?
echo "$smiles"

exit $RC
