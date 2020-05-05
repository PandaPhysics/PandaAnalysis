#!/usr/bin/env python

from glob import glob
from os import stat,getenv,system,path
from PandaCore.Tools.script import *
from re import sub, match
from sys import argv
import argparse
import os

s = '123'
b = '678'
site = 'root://cmsxrootd.fnal.gov/'
dtype='MC'
xsec='100'

with open('out', 'w') as out_file:
    with open('tem', 'r') as in_file:
        for line in in_file:
            out_file.write('{0:<25} {2:<10} {3:<15} {1}\n'.format(s,site+line,dtype,xsec))
out_file.close()
     
