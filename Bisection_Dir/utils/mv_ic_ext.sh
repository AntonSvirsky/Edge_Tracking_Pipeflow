#!/bin/bash
SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
echo 'Extracting IC mat_vec files'
echo "$SCRIPTPATH"
cd "$SCRIPTPATH"

../utils/p2m.out < input1.inp > /dev/null
mv mat_vec.cdf mv_1p.cdf
echo 'mv_1p.cdf Done'

../utils/p2m.out < input2.inp >  /dev/null
mv mat_vec.cdf mv_2p.cdf
echo 'mv_2p.cdf Done'
rm mat_maxmin
rm input1.inp
rm input2.inp

