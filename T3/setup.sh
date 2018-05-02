#!/bin/bash

export PATH=${PATH}:${CMSSW_BASE}/src/PandaCore/bin/

export PANDA="${CMSSW_BASE}/src/PandaAnalysis"
export PANDA_CFG="http://t3serv001.mit.edu/~snarayan/histcatalog/010_data.cfg" 
export PANDA_FLATDIR="/home/snarayan/home000/store/panda/v_010_0/"
mkdir -p $PANDA_FLATDIR
export PANDA_XSECS="/home/snarayan/cms/cmssw/analysis/MonoTop_Xsec/"

export SUBMIT_TMPL="skim_monojet_tmpl.py"
export SUBMIT_NAME="v_010_0"
export SUBMIT_WORKDIR="/data/t3serv014/snarayan/jobs/"${SUBMIT_NAME}"/work/"
export SUBMIT_LOGDIR="/data/t3serv014/snarayan/jobs/"${SUBMIT_NAME}"/logs/"
export SUBMIT_LOCKDIR="/data/t3serv014/snarayan/jobs/"${SUBMIT_NAME}"/locks/"
export SUBMIT_OUTDIR="/mnt/hadoop/scratch/snarayan/panda/"${SUBMIT_NAME}"/batch/"
mkdir -p $SUBMIT_WORKDIR $SUBMIT_OUTDIR $SUBMIT_LOGDIR $SUBMIT_LOCKDIR

export SUBMIT_CONFIG=T2

export SUBMIT_URGENT=0

