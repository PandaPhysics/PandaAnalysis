#!/bin/bash

export PATH=${PATH}:${CMSSW_BASE}/src/PandaCore/bin/

export PANDA="${CMSSW_BASE}/src/PandaAnalysis"
export PANDA_CFG="http://t3serv001.mit.edu/~bmaier/stuff/DoubleMuon_2017.cfg" 
export PANDA_FLATDIR="/scratch5/bmaier/hbb/2017/vbfg_v_012_v1"
mkdir -p $PANDA_FLATDIR

export SUBMIT_USER=$USER
export SUBMIT_TMPL="skim_dimuon2017_tmpl.py"
export SUBMIT_NAME="dimuon_v_012_v1_test"
export SUBMIT_WORKDIR="/data/t3home000/bmaier/condor/"${SUBMIT_NAME}"/work/"
export  SUBMIT_LOGDIR="/data/t3home000/bmaier/condor/"${SUBMIT_NAME}"/logs/"
export SUBMIT_LOCKDIR="/data/t3home000/bmaier/condor/"${SUBMIT_NAME}"/locks/"
export  SUBMIT_OUTDIR="/mnt/hadoop/scratch/bmaier/panda/"${SUBMIT_NAME}"/batch/"
mkdir -p $SUBMIT_WORKDIR $SUBMIT_LOCKDIR $SUBMIT_LOGDIR $SUBMIT_OUTDIR
export SUBMIT_CONFIG=T2

export SUBMIT_URGENT=0

export SUBMIT_USER=$USER
export SUBMIT_REPORT="t3serv004.mit.edu:5000"
export SUBMIT_TEXTLOCK=0
