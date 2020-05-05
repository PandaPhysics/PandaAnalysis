#!/usr/bin/env python

from glob import glob
from os import stat,getenv,system,path
from PandaCore.Tools.script import *
from re import sub, match
from sys import argv
import argparse
import os

args = parse(
             '--Name',
             '--dataset',
             '--output',
             ('--type',{'default':'MC'}),
             ('--xces',{'default':'1'}))
             
das_host='https://cmsweb.cern.ch'

if not args.dataset:
  args.dataset="/ZJetsToNuNu_HT-400To600_13TeV-madgraph/RunIIAutumn18NanoAODv5-Nano1June2019_102X_upgrade2018_realistic_v19-v1/NANOAODSIM"

if not args.output:
  args.output ="{0}{1}".format(args.Name,'.cfg')

query = "file dataset="+args.dataset

os.system('dasgoclient -query \'%s\' >& tem'%query)

s = args.Name
site = 'root://cmsxrootd.fnal.gov/'
dtype=args.type
xsec='1'
out=args.output
with open(out, 'a') as out_file:
    with open('tem', 'r') as in_file:
        cnt=0;
        for line in in_file:
            cnt+=1;
            s="{0}{1}".format(args.Name,'_%d'%cnt)
            out_file.write('{0:<28} {2:<10} {3:<15} {1}'.format(s,site+line,dtype,xsec))
out_file.close()
os.system('rm tem')
