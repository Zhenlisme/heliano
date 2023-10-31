#!/bin/bash

mkdir -p $PREFIX/bin/
cp -rf $SRC_DIR/ $PREFIX/bin/

## set for heliano environment
cd $PREFIX/bin/work

SCRIT_DIR_PATH=`pwd`

BCHECK=${SCRIT_DIR_PATH}/heliano_bcheck.R
FISHER=${SCRIT_DIR_PATH}/heliano_fisher.R
HMMmodel=${SCRIT_DIR_PATH}/RepHel.hmm
Headermodel=${SCRIT_DIR_PATH}/tclcv.txt
myPYTHON_PATH=${PYTHON}
SPLIT=${SCRIT_DIR_PATH}/SplitJoint.R
SORT=${SCRIT_DIR_PATH}/Sort.sh

cp heliano.py heliano

sed -i "s|_INTERPRETERPYTHON_PATH_|${myPYTHON_PATH}|" heliano

sed -i "s|_HMM_|${HMMmodel}|" heliano

sed -i "s|_HEADER_|${Headermodel}|" heliano

sed -i "s|_FISHER_|${FISHER}|" heliano

sed -i "s|_BOUNDARY_|${BCHECK}|" heliano

sed -i "s|_SPLIT_JOINT_|${SPLIT}|" heliano

sed -i "s|_SORTPRO_|${SORT}|" heliano

mv heliano $PREFIX/bin/

chmod 755 $PREFIX/bin/heliano

## set for heliano_cons environment
cp heliano_cons.py heliano_cons

sed -i "s|_INTERPRETERPYTHON_PATH_|${myPYTHON_PATH}|" heliano_cons

mv heliano_cons $PREFIX/bin/

chmod 755 $PREFIX/bin/heliano_cons

