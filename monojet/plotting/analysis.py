#!/usr/bin/env python

from PandaCore.Tools.script import *
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
import PandaAnalysis.VBF.PandaSelection as sel
#import PandaAnalysis.VBF.TriggerSelection as sel
from PandaCore.Drawers.plot_utility import *

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_FLATDIR')+'/' 
dataDir = baseDir

args = parse('--outdir', ('--cut', {'default':'1==1'}),
             '--region', ('--cat', {'default':'loose'}))
lumi = 36000.
blind=True
region = args.region

### DEFINE REGIONS ###

cut = tAND(sel.cuts[args.region],args.cut)
if args.cat!='loose':
    cut = tAND(cut,getattr(sel,args.cat))

### LOAD PLOTTING UTILITY ###
BLIND=False

plot = PlotUtility()
plot.Stack(True)
if not(BLIND and 'signal' in region):
    plot.Ratio(True)
    plot.FixRatio(1)
plot.SetTDRStyle("vbf")
plot.InitLegend()
plot.DrawMCErrors(True)
plot.AddCMSLabel()
plot.cut = cut
plot.SetEvtNum("eventNumber")
plot.SetLumi(lumi/1000)
plot.AddLumiLabel(True)
plot.do_overflow = True
plot.do_underflow = True

if args.cat == 'cnc':
    weight = sel.weights_cnc[region]%lumi
else:
    weight = sel.weights[region]%lumi
plot.mc_weight = weight

### DEFINE PROCESSES ###
if 'signal' in args.region:
    zjets         = Process('Z+jets [QCD]',root.kZjets)
else:
    zjets         = Process('Z+jets [QCD]',root.kExtra1)
wjets         = Process('W+jets [QCD]',root.kWjets)
zjets_ewk     = Process('Z+jets [EWK]',root.kExtra3)
wjets_ewk     = Process('W+jets [EWK]',root.kExtra2)
diboson       = Process('Diboson',root.kDiboson)
top           = Process('Top',root.kTTbar)
qcd           = Process("QCD",root.kQCD)
data          = Process("Data",root.kData)
vbf           = Process('VBF H(inv)',root.kSignal1)

### ASSIGN FILES TO PROCESSES ###
if 'signal' in region:
    zjets.add_file(baseDir+'ZtoNuNu.root')
    zjets_ewk.add_file(baseDir+'ZtoNuNu_EWK.root')
    data.add_file(baseDir+'MET.root')
    data.additional_cut = sel.triggers['met']
else:
    zjets.add_file(baseDir+'ZJets.root')
    zjets_ewk.add_file(baseDir+'ZJets_EWK.root')
    if 'muon' in region:
        data.add_file(baseDir+'MET.root')
        data.additional_cut = sel.triggers['met']
    elif 'electron' in region:
        data.add_file(baseDir+'SingleElectron.root')
        data.additional_cut = sel.triggers['ele']
wjets.add_file(baseDir+'WJets.root')
wjets_ewk.add_file(baseDir+'WJets_EWK.root')
diboson.add_file(baseDir+'Diboson.root')
top.add_file(baseDir+'TTbar.root');
top.add_file(baseDir+'SingleTop.root');
qcd.add_file(baseDir+'QCD.root')
vbf.add_file(baseDir+'vbfHinv_m125.root')

if 'single' in region:
    processes = [qcd,diboson,zjets,zjets_ewk,top,wjets_ewk,wjets,data]
elif 'di' in region:
    processes = [qcd,diboson,wjets,wjets_ewk,top,zjets_ewk,zjets,data]
else:
    processes = [qcd,diboson,top,wjets_ewk,zjets_ewk,wjets,zjets,vbf]
    if not BLIND:
        processes.append(data)

for p in processes:
    plot.add_process(p)

recoilBins = [200., 230., 260.0, 290.0, 320.0, 350.0, 390.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0]
nRecoilBins = len(recoilBins)-1

if 'electron' in region:
    lep = 'e'
    lepton = 'electron'
else:
    lep = '#mu'
    lepton = 'muon'

### CHOOSE DISTRIBUTIONS, LABELS ###
if 'signal' in region:
    recoil=VDistribution("pfmet",recoilBins,"PF MET [GeV]","Events/GeV")
    plot.add_distribution(FDistribution('dphipfmet',0.5,3.2,20,'min #Delta#phi(U,jets)','Events'))
elif 'single' in region:
    recoil=VDistribution("pfUWmag",recoilBins,"PF U(%s) [GeV]"%(lep),"Events/GeV")
    plot.add_distribution(FDistribution('%sPt[0]'%lepton,0,1000,20,'Leading %s p_{T} [GeV]'%lep,'Events/40 GeV'))
    plot.add_distribution(FDistribution('%sEta[0]'%lepton,-2.5,2.5,20,'Leading %s #eta'%lep,'Events/bin'))
#    plot.add_distribution(FDistribution('dphipfUW',0.5,3.2,20,'min #Delta#phi(U,jets)','Events'))
elif any([x in region for x in ['dielectron','dimuon']]):
    recoil=VDistribution("pfUZmag",recoilBins,"PF U(%s%s) [GeV]"%(lep,lep),"Events/GeV")
    plot.add_distribution(FDistribution('diLepMass',60,120,20,'m_{ll} [GeV]','Events/3 GeV'))
    plot.add_distribution(FDistribution('%sPt[0]'%lepton,0,1000,20,'Leading %s p_{T} [GeV]'%lep,'Events/40 GeV'))
    plot.add_distribution(FDistribution('%sEta[0]'%lepton,-2.5,2.5,20,'Leading %s #eta'%lep,'Events/bin'))
    plot.add_distribution(FDistribution('%sPt[1]'%lepton,0,1000,20,'Subleading %s p_{T} [GeV]'%lep,'Events/40 GeV'))
    plot.add_distribution(FDistribution('%sEta[1]'%lepton,-2.5,2.5,20,'Subleading %s #eta'%lep,'Events/bin'))
#    plot.add_distribution(FDistribution('dphipfUZ',0.5,3.2,20,'min #Delta#phi(U,jets)','Events'))

plot.add_distribution(recoil)

#plot.add_distribution(FDistribution('barrelHT',0,1000,20,'Barrel H_{T} [GeV]','Events/50 GeV'))
#plot.add_distribution(FDistribution('barrelJet1Pt',0,1000,20,'Barrel jet 1 p_{T} [GeV]','Events/50 GeV'))
plot.add_distribution(FDistribution('jot12Mass',0,4000,20,'m_{jj} [GeV]','Events/200 GeV'))
plot.add_distribution(FDistribution('jot12DEta',0,10,20,'#Delta#eta(j_{1},j_{2})','Events'))
plot.add_distribution(FDistribution("fabs(jot12DPhi)",0,3.142,20,"#Delta #phi leading jets","Events",filename='jot12DPhi'))
plot.add_distribution(FDistribution("jotEta[0]",-5,5,20,"Jet 1 #eta","Events"))
plot.add_distribution(FDistribution("jotEta[1]",-5,5,20,"Jet 2 #eta","Events"))
plot.add_distribution(FDistribution("jotPt[0]",80,500,20,"Jet 1 p_{T} [GeV]","Events"))
plot.add_distribution(FDistribution("jotPt[1]",40,500,20,"Jet 2 p_{T} [GeV]","Events"))
plot.add_distribution(FDistribution("1",0,2,1,"dummy","dummy"))

### DRAW AND CATALOGUE ###
region = '%s/%s'%(args.cat,region)
plot.draw_all(args.outdir+'/'+region+'_')
