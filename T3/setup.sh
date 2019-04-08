#!/bin/bash

export PATH=${PATH}:${CMSSW_BASE}/src/PandaCore/bin/

export PANDA="${CMSSW_BASE}/src/PandaAnalysis"
export PANDA_CFG="http://t3serv001.mit.edu/~bmaier/stuff/DoubleMuon_2018.cfg" 
export PANDA_FLATDIR="/scratch5/bmaier/hbb/2018/v_013_full"
mkdir -p $PANDA_FLATDIR

export SUBMIT_USER=$USER
export SUBMIT_TMPL="skim_wlnhbb2018_tmpl.py"
export SUBMIT_NAME="v_013_full"
export SUBMIT_WORKDIR="/data/t3home000/bmaier/condor/"${SUBMIT_NAME}"/work/"
export  SUBMIT_LOGDIR="/data/t3home000/bmaier/condor/"${SUBMIT_NAME}"/logs/"
export SUBMIT_LOCKDIR="/data/t3home000/bmaier/condor/"${SUBMIT_NAME}"/locks/"
export  SUBMIT_OUTDIR="/mnt/hadoop/scratch/bmaier/panda/"${SUBMIT_NAME}"/batch/"
mkdir -p $SUBMIT_WORKDIR $SUBMIT_LOCKDIR $SUBMIT_LOGDIR $SUBMIT_OUTDIR
export SUBMIT_CONFIG=T2


#export SUBMIT_NPY="/mnt/hadoop/scratch/snarayan/deep/"${SUBMIT_NAME}"/"
#export SUBMIT_NPY="/data/t3serv014/snarayan/deep/"${SUBMIT_NAME}"/"
#mkdir -p $SUBMIT_NPY/train $SUBMIT_NPY/test $SUBMIT_NPY/validate

export SUBMIT_URGENT=0

export SUBMIT_USER=$USER
export SUBMIT_REPORT="t3serv004.mit.edu:5000"
export SUBMIT_TEXTLOCK=0
