#!/usr/bin/env python

from re import sub
from os import system,getenv,path
from time import clock,time
import json

from PandaCore.Utils.root import root
from PandaCore.Tools.Misc import *
from PandaCore.Utils.load import *
import PandaAnalysis.T3.job_utilities as utils
from PandaAnalysis.Flat.analysis import * 

Load('PandaAnalyzer')
data_dir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/'

def fn(input_name, isData, full_path):
    
    # now we instantiate and configure the analyzer
    a = monotop()
    a.inpath = input_name
    a.outpath = utils.input_to_output(input_name)
    a.datapath = data_dir
    a.isData = isData
    utils.set_year(a, 2016)
    a.processType = utils.classify_sample(full_path, isData)	

    skimmer = root.pa.PandaAnalyzer(a)
    skimmer.AddPresel(root.pa.MonotopSel())

    return utils.run_PandaAnalyzer(skimmer, isData, a.outpath)


if __name__ == "__main__":
    utils.wrapper(fn, post_fn=utils.BDTAdder())

