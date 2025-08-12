#!/bin/bash
SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
echo "$SCRIPTPATH"
cd "$SCRIPTPATH"
./prim2matlab.out < input.inp >/dev/null
mkdir -p mat_vec_out
mv mat_vec0* mat_vec_out/
echo "Extraction complete"
rm input.inp 
rm mat_maxmin 
#mv mat_vec.cdf mat_vec_1.cdf
