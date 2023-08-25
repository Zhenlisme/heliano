#!/bin/bash

mkdir -p $PREFIX/bin/
cp -rf $SRC_DIR/ $PREFIX/bin/

## set for HELA environment
cd $PREFIX/bin/work

SCRIT_DIR_PATH=`pwd`

BCHECK=${SCRIT_DIR_PATH}/HELA_bcheck.R
FISHER=${SCRIT_DIR_PATH}/HELA_fisher.R
HMMmodel=${SCRIT_DIR_PATH}/RepHel.hmm
Headermodel=${SCRIT_DIR_PATH}/tclcv.txt
myPYTHON_PATH=${PYTHON}
SPLIT=${SCRIT_DIR_PATH}/SplitJoint.R

cp HELA.py HELA

sed -i "s|_INTERPRETERPYTHON_PATH_|${myPYTHON_PATH}|" HELA

sed -i "s|_HMM_|${HMMmodel}|" HELA

sed -i "s|_HEADER_|${Headermodel}|" HELA

sed -i "s|_FISHER_|${FISHER}|" HELA

sed -i "s|_BOUNDARY_|${BCHECK}|" HELA

sed -i "s|_SPLIT_JOINT_|${SPLIT}|" HELA

mv HELA $PREFIX/bin/

chmod 755 $PREFIX/bin/HELA

## set for HELA_cons environment
cp HELA_cons.py HELA_cons

sed -i "s|_INTERPRETERPYTHON_PATH_|${myPYTHON_PATH}|" HELA_cons

mv HELA_cons $PREFIX/bin/

chmod 755 $PREFIX/bin/HELA_cons

