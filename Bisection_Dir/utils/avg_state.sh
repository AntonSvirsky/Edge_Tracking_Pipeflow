#!/bin/bash

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
echo 'Extracting IC mat_vec files'
echo "$SCRIPTPATH"
cd "$SCRIPTPATH"

../utils/avg_state.out < input0.inp > out.txt /dev/null
mv state0000.cdf.dat s_avg.cdf.dat
echo 'avg state created'
rm input0.inp

