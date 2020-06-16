#include "../interface/CommonOps.h"

using namespace pa;
using namespace std;
using namespace panda;

/*
void shiftMET(const RecoMet& met, TLorentzVector& v, shiftjes shift) 
{
  float pt;
  float phi;
  switch (shift) {
    case shiftjes::kNominal:
      pt = met.pt;
      phi = met.phi;
      break;
    case shiftjes::kJESTotalUp:
      pt = met.ptCorrUp;
      phi = met.phiCorrUp;
      break;
    case shiftjes::kJESTotalDown:
      pt = met.ptCorrDown;
      phi = met.phiCorrDown;
      break;
    default:
      logger.error("shiftMET", "Unknown JES type!");
      exit(1);
  }

  v.SetPtEtaPhiM(pt, 0, phi, 0);
}
*/

void TriggerOp::do_init(Registry& registry) 
{
  vector<TString> paths; 
  if (analysis.isData || analysis.mcTriggers) {
    // MET
    if (analysis.year == 2016) {
      paths = {
            "HLT_PFMET170_NoiseCleaned",
            "HLT_PFMET170_HBHECleaned",
            "HLT_PFMET170_JetIdCleaned",
            "HLT_PFMET170_NotCleaned",
            "HLT_PFMET170_HBHE_BeamHaloCleaned",
            "HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight",
            "HLT_PFMETNoMu110_NoiseCleaned_PFMHTNoMu110_IDTight",
            "HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight",
            "HLT_PFMETNoMu90_PFMHTNoMu90_IDTight",
            "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight",
            "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight",
            "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"
      };
    } else if (analysis.year == 2017) { 
        paths = {
          "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",
          "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
          "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight",
          "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight",
        };
    } else if (analysis.year == 2018) { 
        paths = {
          "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60",
          "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",
          "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight",
          "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
          "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight",
          "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight",
        };
    }
    triggerHandlers[kMETTrig].addTriggers(paths);

    // SingleEle    
    if (analysis.complicatedLeptons) {
      if (analysis.year == 2016) {
        paths = {
          "HLT_Ele25_eta2p1_WPTight_Gsf",
          "HLT_Ele27_eta2p1_WPLoose_Gsf",
          "HLT_Ele27_WPTight_Gsf",
          "HLT_Ele30_WPTight_Gsf",
          "HLT_Ele35_WPLoose_Gsf",
          "HLT_Ele27_WP85_Gsf",
          "HLT_Ele27_WPLoose_Gsf",
          "HLT_Ele105_CaloIdVT_GsfTrkIdT",
          "HLT_Ele115_CaloIdVT_GsfTrkIdT",
          "HLT_Ele27_eta2p1_WPTight_Gsf",
          "HLT_Ele32_eta2p1_WPTight_Gsf",
          "HLT_ECALHT800"
        };
      } else if (analysis.year == 2017) {
        paths = {
          "HLT_Ele115_CaloIdVT_GsfTrkIdT",
          "HLT_Ele27_WPTight_Gsf",
          "HLT_Ele32_WPTight_Gsf",
          "HLT_Ele35_WPTight_Gsf",
          "HLT_Ele32_WPTight_Gsf_L1DoubleEG",
          "HLT_Photon200"
        };
      } else if (analysis.year == 2018) {
        paths = {
          "HLT_Ele115_CaloIdVT_GsfTrkIdT",
          "HLT_Ele27_WPTight_Gsf",
          "HLT_Ele28_WPTight_Gsf",
          "HLT_Ele32_WPTight_Gsf",
          "HLT_Ele35_WPTight_Gsf",
          "HLT_Ele38_WPTight_Gsf",
          "HLT_Ele40_WPTight_Gsf"
          "HLT_Ele32_WPTight_Gsf_L1DoubleEG",
          "HLT_Photon200"
        };
      }
    } else {
      if (analysis.year == 2016) {
        paths = {
          "HLT_Ele27_WP85_Gsf",
          "HLT_Ele27_WPLoose_Gsf",
          "HLT_Ele105_CaloIdVT_GsfTrkIdT",
          "HLT_Ele27_WPTight_Gsf",
          "HLT_Ele30_WPTight_Gsf",
          "HLT_Ele27_eta2p1_WPTight_Gsf",
          "HLT_Ele32_eta2p1_WPTight_Gsf",
          "HLT_Ele35_WPLoose_Gsf",
          "HLT_ECALHT800"
        };
      } else if (analysis.year == 2017) {
        paths = {
          "HLT_Ele35_WPTight_Gsf",
          "HLT_Ele38_WPTight_Gsf",
          "HLT_Ele40_WPTight_Gsf"
        };
      } else if (analysis.year == 2018) {
        paths = {
          "HLT_Ele35_WPTight_Gsf",
          "HLT_Ele38_WPTight_Gsf",
          "HLT_Ele40_WPTight_Gsf"
        };
      }
    }
    triggerHandlers[kSingleEleTrig].addTriggers(paths);
    
    // single pho
    paths = {
          "HLT_Photon175",
          "HLT_Photon200",
          "HLT_Photon165_HE10",
          "HLT_Photon36_R9Id90_HE10_IsoM",
          "HLT_Photon50_R9Id90_HE10_IsoM",
          "HLT_Photon75_R9Id90_HE10_IsoM",
          "HLT_Photon90_R9Id90_HE10_IsoM",
          "HLT_Photon120_R9Id90_HE10_IsoM",
          "HLT_Photon165_R9Id90_HE10_IsoM",
          "HLT_Photon300_NoHE",
          "HLT_ECALHT800",
          "HLT_CaloJet500_NoJetID"
    };
    triggerHandlers[kSinglePhoTrig].addTriggers(paths);

    // Single muon
    if (analysis.complicatedLeptons || analysis.recalcECF) { // either comp lepton or tnp
      if (analysis.year == 2016) {
        paths = {
          "HLT_IsoMu24",
          "HLT_IsoTkMu24",
          "HLT_IsoMu22",
          "HLT_IsoTkMu22",
          "HLT_Mu45_eta2p1",
          "HLT_Mu50"
        };
      } else if (analysis.year == 2017) {
        paths = {
          "HLT_IsoMu24",      
          "HLT_IsoMu27",
          "HLT_IsoMu30"
          "HLT_Mu50"
        };
      } else if (analysis.year == 2018) {
        paths = {
          "HLT_IsoMu24",      
          "HLT_IsoMu27",
          "HLT_IsoMu30"
          "HLT_Mu50"
        };
      }
    } else {
      if (analysis.year == 2016) {
        paths = {
          "HLT_IsoMu20",
          "HLT_IsoMu22",
          "HLT_IsoMu24",
        };
      } else if (analysis.year == 2017) {
        paths = {
          "HLT_IsoMu24",
          "HLT_IsoMu27"
        };
      } else if (analysis.year == 2018) {
        paths = {
          "HLT_IsoMu24",
          "HLT_IsoMu27"
        };
      }
    }
    triggerHandlers[kSingleMuTrig].addTriggers(paths);

    // double muon
    if (analysis.year==2016) { 
      paths = {
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
        "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
	"HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL"
      };
    } else if (analysis.year==2017) {
      paths = {
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
        "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8",
        "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8"
      };
    } else if (analysis.year==2018) {
      paths = {
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
        "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8",
        "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8"
      };
    }
    triggerHandlers[kDoubleMuTrig].addTriggers(paths);

    // double ele
    if (analysis.year==2016) 
      paths = {
            "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            "HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf"
      };
    else if (analysis.year==2017)
      paths = {
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
        "HLT_DiEle27_WPTightCaloOnly_L1DoubleEG",
        "HLT_DoubleEle25_CaloIdL_MW",
        "HLT_DoubleEle33_CaloIdL_MW",
        "HLT_DoublePhoton70"
      };
    else if (analysis.year==2018)
      paths = {
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
        "HLT_DiEle27_WPTightCaloOnly_L1DoubleEG",
        "HLT_DoubleEle25_CaloIdL_MW",
        "HLT_DoubleEle27_CaloIdL_MW",
        "HLT_DoubleEle33_CaloIdL_MW",
        "HLT_DoublePhoton70"
      };
    triggerHandlers[kDoubleEleTrig].addTriggers(paths);
    
    // emu
    if (analysis.year==2016) {
      paths = {
        "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
        "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL"
      };
    } else if (analysis.year==2017) {
      paths = {
        "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
        "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL"
      };
    } else if (analysis.year==2018) {
      paths = {
        "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL"
      };
    }
    triggerHandlers[kEMuTrig].addTriggers(paths);

    // JetHT
    paths = {
          "HLT_PFHT650",
          "HLT_PFHT900",
          "HLT_PFJet500",
          "HLT_PFJet450",
          "HLT_PFJet320",
    };
    triggerHandlers[kJetHTTrig].addTriggers(paths);

    paths = {
          "HLT_Mu8_TrkIsoVVL",
          "HLT_Mu17_TrkIsoVVL"
    };
    triggerHandlers[kMuFakeTrig].addTriggers(paths);

    paths = {
          "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30",
          "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30",
          "HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30",
          "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30"
    };
    triggerHandlers[kEleFakeTrig].addTriggers(paths);


    // VBF+Pho
    if (analysis.year == 2016){
      paths = {
	"HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF"
      };
    }
    else if (analysis.year == 2017){
      paths = {
	"HLT_DiJet110_35_Mjj650_PFMET110",
	"HLT_DiJet110_35_Mjj650_PFMET120",
	"HLT_DiJet110_35_Mjj650_PFMET130",
	"HLT_PFMET120_PFMHT120_IDTight",
	"HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
	"HLT_Photon200"
      };
    }
    else if (analysis.year == 2018){
	paths = {
	"HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3",
	"HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3",
	"HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3",
	"HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3"
      };
    }
    triggerHandlers[kVBFPhoTrig].addTriggers(paths);


    // PhoFakeRate
    paths = {
      "HLT_Photon75"
    };
    triggerHandlers[kPhoFakeRate].addTriggers(paths);



/*
    for (auto &th : triggerHandlers) {
      unsigned N = th.paths.size();

      for (unsigned i = 0; i != N; i++) {
        unsigned panda_idx = event.registerTrigger(th.paths.at(i));
        th.indices[i] = panda_idx;
      }
    }
*/
  }
}
/*
void TriggerOp::checkEle32()
{
//  auto& filter1Objects = event.triggerObjects.filterObjects("hltEle32L1DoubleEGWPTightGsfTrackIsoFilter");
//  auto& filter2Objects = event.triggerObjects.filterObjects("hltEGL1SingleEGOrFilter");
  if (event.nTrigObj<=0) 
    return;

  bool ele_pass=false;
  TLorentzVector filter1ObjectP4;
  for(auto &trgobj : event.TrigObj){
      if(trgobj.id==11 && (trgobj.filterBits>>1)%2==1) {
         filter1ObjectP4.SetPtEtaPhiM(trgobj.pt,trgobj.eta,trgobj.phi,511e-6);
         ele_pass=true;
      }
  }
  if (!ele_pass) return;

  bool matchedToTriggerObject = false;
  for (auto& ele : event.Electron) {
    TLorentzVector temE;
    temE.SetPtEtaPhiM(ele.pt,ele.eta,ele.phi,ele.mass);
    if (!ele.mvaFall17V1Iso_WP90) 
      continue;
    if (filter1ObjectP4.DeltaR(temE)<0.1) {
      matchedToTriggerObject = true;
      break;
    }
  }
  if (matchedToTriggerObject) { 
    gt.trigger |= (1 << kSingleEleTrig);
  }
}
*/
void TriggerOp::do_execute()
{

  if(analysis.year==2018){ if(event.HLT.Ele32_WPTight_Gsf || event.HLT.Ele115_CaloIdVT_GsfTrkIdT) gt.trigger+=1;}
  else { if(event.HLT.Ele35_WPTight_Gsf || event.HLT.Ele115_CaloIdVT_GsfTrkIdT) gt.trigger+=1;}
  if(event.HLT.IsoMu27) gt.trigger+=2;
  if(event.HLT.PFMETNoMu120_PFMHTNoMu120_IDTight || event.HLT.PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60) gt.trigger+=4;
  if(event.HLT.Photon200) gt.trigger+=8;  
  if(event.HLT.PFHT1050) gt.trigger+=16;
  if(event.HLT.PFHT180) gt.trigger+=32;
  if(event.HLT.PFHT430) gt.trigger+=64;
  if(event.HLT.PFHT590) gt.trigger+=128;
  if(event.HLT.PFHT780) gt.trigger+=256;
  if(event.HLT.PFHT250) gt.trigger+=512;
  if(event.HLT.PFHT350) gt.trigger+=1024;  
}

void GlobalOp::do_execute()
{
  if (cfg.DEBUG > 5) {
    logger.debug("PandaAnalyzer::Run::Dump","");
    event.print(std::cout, 2);
    std::cout << std::endl;
    logger.debug("PandaAnalyzer::Run::Dump","");
    event.Photon.print(std::cout, 2);
    std::cout << std::endl;
    logger.debug("PandaAnalyzer::Run::Dump","");
    event.Muon.print(std::cout, 2);
    std::cout << std::endl;
    logger.debug("PandaAnalyzer::Run::Dump","");
    event.Electron.print(std::cout, 2);
    std::cout << std::endl;
    logger.debug("PandaAnalyzer::Run::Dump","");
    event.Jet.print(std::cout, 2);
    std::cout << std::endl;
    logger.debug("PandaAnalyzer::Run::Dump","");
    event.MET.print(std::cout, 2);
    std::cout << std::endl;
  }
  // event info


//  cout << "npdf=" << event.nLHEPdfWeight << " nqcd=" << event.nLHEScaleWeight <<" nobj=" << event.nTrigObj << endl;

  
  if(!analysis.isData){
  gt.npdf = event.nLHEPdfWeight;
  gt.nscale = event.nLHEScaleWeight;
  
  float temx2=0;
  float ave=0;

//  for(int k=0; k<gt.nscale; k++){
//    cout << "qcd" << k << "=" << event.LHEScaleWeight[k].VAL << endl;
//  }


  gt.pdfup=1;
  gt.pdfdow=1;
  for(int k=0; k<9; k++)
     gt.qcd[k]=1;

  for(int k=0; k<gt.npdf; k++){
    ave=ave+event.LHEPdfWeight[k].VAL;
  }
  ave=ave/gt.npdf;

  for(int k=0; k<gt.npdf; k++){
    temx2=temx2+pow((event.LHEPdfWeight[k].VAL-ave),2);
  }

  if(gt.npdf>1) temx2=sqrt(temx2/float(gt.npdf-1));
  else temx2=0;

    gt.pdfup=1+temx2;
    gt.pdfdow=1-temx2;


  if(gt.nscale<=9){
   for(int k=0; k<gt.nscale; k++){
      gt.qcd[k]=event.LHEScaleWeight[k].VAL;
    }
  }

  if(gt.nscale==44){
//   int mmap[9]={0,3,6,1,4,7,2,5,8}; //this sample is messed up....(WZTo3LNu_3Jets_MLL4-50)
   int ll=0;
   for(int k=4; k<gt.nscale; k=k+5){
      ll++; 
      gt.qcd[ll]=event.LHEScaleWeight[k].VAL;
    }
  } 

/*
   if(gt.qcd[0]!=1){
    for(int k=0; k<9; k++){
      cout << "qcd" << k <<"=" << event.LHEScaleWeight[k].VAL << endl;
    }
   }
*/ 
  }


  gt.mcWeight = event.genWeight;
  gt.runNumber = event.run;
  gt.lumiNumber = event.luminosityBlock;
  gt.eventNumber = event.event;
  gt.isData = analysis.isData ?  1 : 0; 
  gt.npv = event.PV.npvs;
  gt.rho = event.fixedGridRhoFastjetAll;
  gt.metFilter = (event.Flag.METFilters);
  gt.comFilter = 0;
  if(analysis.year==2016){
    if(analysis.isData){
      if(event.Flag.goodVertices && event.Flag.globalSuperTightHalo2016Filter && event.Flag.HBHENoiseFilter && event.Flag.HBHENoiseIsoFilter && event.Flag.EcalDeadCellTriggerPrimitiveFilter && event.Flag.BadPFMuonFilter && event.Flag.eeBadScFilter) gt.comFilter = 1;
    }
   else{
      if(event.Flag.goodVertices && event.Flag.globalSuperTightHalo2016Filter && event.Flag.HBHENoiseFilter && event.Flag.HBHENoiseIsoFilter && event.Flag.EcalDeadCellTriggerPrimitiveFilter && event.Flag.BadPFMuonFilter) gt.comFilter = 1;
   }
  }
  if(analysis.year==2017 || analysis.year==2018){
    if(analysis.isData){
      if(event.Flag.goodVertices && event.Flag.globalSuperTightHalo2016Filter && event.Flag.HBHENoiseFilter && event.Flag.HBHENoiseIsoFilter && event.Flag.EcalDeadCellTriggerPrimitiveFilter && event.Flag.BadPFMuonFilter && event.Flag.eeBadScFilter && event.Flag.ecalBadCalibFilterV2) gt.comFilter = 1;
    }
    else{
      if(event.Flag.goodVertices && event.Flag.globalSuperTightHalo2016Filter && event.Flag.HBHENoiseFilter && event.Flag.HBHENoiseIsoFilter && event.Flag.EcalDeadCellTriggerPrimitiveFilter && event.Flag.BadPFMuonFilter && event.Flag.ecalBadCalibFilterV2) gt.comFilter =1;
    }
  }

//  gt.metFilter = (gt.metFilter==1 && !event.metFilters.badPFMuons) ? 1 : 0;
//  gt.metFilter = (gt.metFilter==1 && !event.metFilters.badChargedHadrons) ? 1 : 0;
  if (!analysis.isData) {
    gt.pu = event.Pileup.nTrueInt;
    gt.sf_npv = utils.getCorr(cNPV, gt.npv);
    gt.sf_pu = utils.getCorr(cPU, gt.pu);
    gt.lhevpt=event.LHE.Vpt;
    gt.lhenjet = event.LHE.Njets;
    gt.lheht = event.LHE.HT;
  }

  
  
}
/*
template <typename TREE>
void BaseGenPOp<TREE>::do_execute()
{
  if (this->event.genParticles.size() > 0) {
    merge_particles(this->event.genParticles);
  } else {
    merge_particles(this->event.genParticlesU);
  }
}
*/
//template class BaseGenPOp<GeneralTree>;

void RecoilOp::do_execute()
{
  TLorentzVector vObj1, vObj2;
  gt.whichRecoil = 0; // -1=photon, 0=MET, 1,2=nLep
//   cout << "in recoil" << endl;
  if (gt.nLooseLep>0) {
    Lepton *lep1 = looseLeps->at(0);
    vObj1.SetPtEtaPhiM(lep1->pt,lep1->eta,lep1->phi,lep1->mass);

    // one lep => W
    METLOOP {
      auto& jets = (*jesShifts)[shift]; 

      jets.vpuppiUW = jets.vpuppiMET + vObj1;
      gt.puppiUWmag[shift] = jets.vpuppiUW.Pt();
      gt.puppiUWphi[shift] = jets.vpuppiUW.Phi();

      jets.vpfUW = jets.vpfMET + vObj1;
      gt.pfUWmag[shift] = jets.vpfUW.Pt();
      gt.pfUWphi[shift] = jets.vpfUW.Phi();
    }

    if (gt.nLooseLep>1 && (*lepPdgId)[0]+(*lepPdgId)[1]==0) {
      // two OS lep => Z
      Lepton *lep2 = looseLeps->at(1);
      vObj2.SetPtEtaPhiM(lep2->pt,lep2->eta,lep2->phi,lep2->mass);

      METLOOP {
        auto& jets = (*jesShifts)[shift]; 
        gt.puppiUZmag[shift] = jets.vpuppiUZ.Pt();
        gt.puppiUZphi[shift] = jets.vpuppiUZ.Phi();

        jets.vpuppiUZ = jets.vpuppiUW + vObj2;

        jets.vpfUZ = jets.vpfUW + vObj2;
        gt.pfUZmag[shift] = jets.vpfUZ.Pt(); 
        gt.pfUZphi[shift] = jets.vpfUZ.Phi(); 
      }

      gt.whichRecoil = 2;
    } else {
      gt.whichRecoil = 1;
    }
  }
  if (gt.nLoosePhoton>0) {
    Photon *pho = (*loosePhos)[0];
    vObj1.SetPtEtaPhiM(pho->pt,pho->eta,pho->phi,0.);

    METLOOP {
      auto& jets = (*jesShifts)[shift]; 

      jets.vpuppiUA = jets.vpuppiMET + vObj1;

      jets.vpfUA = jets.vpfMET + vObj1;
      gt.pfUAmag[shift] = jets.vpfUA.Pt(); 
      gt.pfUAphi[shift] = jets.vpfUA.Phi(); 
    }

    if (gt.nLooseLep==0) {
      gt.whichRecoil = -1;
    }
  }
  if (gt.nLooseLep==0 && gt.nLoosePhoton==0) {
    gt.whichRecoil = 0;
  }
}

void TriggerEffOp::do_execute() 
{
  // trigger efficiencies
  gt.sf_metTrig = utils.getCorr(cTrigMET,gt.pfmetnomu[jes2i(shiftjes::kNominal)]);
  gt.sf_metTrigZmm = utils.getCorr(cTrigMETZmm,gt.pfmetnomu[jes2i(shiftjes::kNominal)]);

  auto* lep0 = looseLeps->size()>0 ? (*looseLeps)[0] : nullptr;
  auto* lep1 = looseLeps->size()>1 ? (*looseLeps)[1] : nullptr;

  if(analysis.year == 2017 || analysis.year == 2018){
  if (gt.nLooseElectron>1 && analysis.complicatedLeptons) {
     double data_eff1=utils.getCorr(cTrigEleDataEff, gt.electronEta[0], gt.electronPt[0]);
     double data_eff2=utils.getCorr(cTrigEleDataEff, gt.electronEta[1], gt.electronPt[1]);
     double mc_eff1 = utils.getCorr(cTrigEleMCEff, gt.electronEta[0], gt.electronPt[0]);
     double mc_eff2 = utils.getCorr(cTrigEleMCEff, gt.electronEta[1], gt.electronPt[1]);
    gt.sf_eleTrig = (1-(1-data_eff1)*(1-data_eff2))/(1-(1-mc_eff1)*(1-mc_eff2));
    if(gt.eleTrigMatch[0]==1)
        gt.sf_eleTrig2 = utils.getCorr(cTrigEle, gt.electronEta[0], gt.electronPt[0]);
    else if(gt.eleTrigMatch[1]==1)
        gt.sf_eleTrig2 = utils.getCorr(cTrigEle, gt.electronEta[1], gt.electronPt[1]);
  } else if (gt.nLooseElectron==1) {
      gt.sf_eleTrig = utils.getCorr(cTrigEle, gt.electronEta[0], gt.electronPt[0]);
      if(gt.eleTrigMatch[0]==1)
          gt.sf_eleTrig2 = utils.getCorr(cTrigEle, gt.electronEta[0], gt.electronPt[0]);
  } // done with ele trig SF
  if (gt.nLooseMuon>=1 && analysis.complicatedLeptons) {
    Muon *mu1=nullptr, *mu2=nullptr;
    float eff1=0;
      eff1 = utils.getCorr(
        cTrigMu,
        fabs(gt.muonEta[0]),
        TMath::Max((float)26.,TMath::Min((float)499.99,gt.muonPt[0]))
      );
      gt.sf_muTrig = eff1;
  } // done with mu trig SF
  }
  if (gt.nLoosePhoton>0 && gt.loosePho1IsTight)
    gt.sf_phoTrig = utils.getCorr(cTrigPho,gt.loosePho1Pt);

//  if (analysis.vbf) {
//    gt.sf_metTrigVBF = utils.getCorr(cVBF_TrigMET,gt.barrelHTMiss);
//    gt.sf_metTrigZmmVBF = utils.getCorr(cVBF_TrigMETZmm,gt.barrelHTMiss);
//  }
}
