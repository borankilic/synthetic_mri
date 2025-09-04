#!/bin/bash

export LD_LIBRARY_PATH=""
"$1" /data/u_kilicb_software/romeo.jl -p "$2" -m "$3" -t $4 -B -o "$5"
echo ROMEO execution through bash script has finished.
