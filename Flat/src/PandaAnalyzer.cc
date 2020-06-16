#include "../interface/PandaAnalyzer.h"
#include "TMath.h"
#include <algorithm>
#include <vector>


using namespace panda;
using namespace std;
using namespace pa;

#define ADDOP(X) ops_all.back()->addSubOp<X>()


PandaAnalyzer::PandaAnalyzer(Analysis* a, int debug_/*=1*/) :
  Analyzer("PandaAnalyzer", a, debug_),
  cfgop(analysis, gt, DEBUG),
  wIDs(v_make_shared<TString>())
{
  if (DEBUG) logger.debug("PandaAnalyzer::PandaAnalyzer","Calling constructor");

  Config& cfg = cfgop.cfg;
  Utils& utils = cfgop.utils;

  gblop = new GlobalOp(event, cfg, utils, gt);
  ops_all.emplace_back(gblop);

  if (DEBUG) logger.debug("PandaAnalyzer::PandaAnalyzer","Adding AnalysisOps");

  // Define analyses - user should not touch below code. 
  preselop = new ContainerOp("pre-sel", event, cfg, utils, gt);
  ops_all.emplace_back(preselop);

  ADDOP(JetCorrOp);
  ADDOP(TriggerOp);
  ADDOP(SimpleLeptonOp);
  ADDOP(ComplicatedLeptonOp);
  ADDOP(SimplePhotonOp);
  ADDOP(ComplicatedPhotonOp);

  ADDOP(FatJetOp);
  ADDOP(JetOp);

  ADDOP(RecoilOp);
  ADDOP(TauOp);

  ADDOP(BTagSFOp);
  ADDOP(BTagWeightOp);
  ADDOP(TriggerEffOp);


  postselop = new ContainerOp("post-sel", event, cfg, utils, gt);
  ops_all.emplace_back(postselop);


  for (auto& op : ops_all)
    op->print();

  if (DEBUG) logger.debug("PandaAnalyzer::PandaAnalyzer","Reading inputs");
  // Read inputs
  getInput();

  cout << "Setting Branch status" << endl;
  event.setStatus(*tIn, {"!*"});
  event.setAddress(*tIn, cfgop.get_inputBranches());

/*
  TH1D* hDTotalMCWeight = static_cast<TH1D*>(static_cast<TH1D*>(fIn->Get("hSumW"))->Clone("hDTotalMCWeight"));
  hDTotalMCWeight->SetDirectory(0);
  TH1D* hDNPUWeight = nullptr; {
    TH1D* hbase = static_cast<TH1D*>(fIn->Get("hNPVTrue"));
    if (hbase == nullptr)
      hbase = static_cast<TH1D*>(fIn->Get("hNPVReco"));
    if (hbase != nullptr) {
      hDNPUWeight = static_cast<TH1D*>(hbase->Clone("hDNPUWeight"));
      hDNPUWeight->SetDirectory(0);
    }
  }

  TTree* tW = static_cast<TTree*>(fIn->Get("weights"));
  if (tW && analysis.processType == kSignal) {
    if (tW->GetEntries()!=377 && tW->GetEntries()!=22) {
      logger.error("PandaAnalyzer::PandaAnalyzer",
          TString::Format("Reweighting failed because only found %u weights!",
                          unsigned(tW->GetEntries())));
      throw runtime_error("");
    }
    TString *id = new TString();
    tW->SetBranchAddress("id",&id);
    unsigned nW = tW->GetEntriesFast();
    for (unsigned iW=0; iW!=nW; ++iW) {
      tW->GetEntry(iW);
      wIDs->push_back(*id);
    }
  } else if (analysis.processType==kSignal) {
    logger.error("PandaAnalyzer::PandaAnalyzer","This is a signal file, but the weights are missing!");
    throw runtime_error("");
  }
  registry.publishConst("wIDs", wIDs);
*/
  // Define outputs
  if (DEBUG) logger.debug("PandaAnalyzer::PandaAnalyzer","Writing outputs");

  makeOutput(); 
/*
  fOut->WriteTObject(hDTotalMCWeight); delete hDTotalMCWeight; hDTotalMCWeight = nullptr;
  if (hDNPUWeight != nullptr) {
    fOut->WriteTObject(hDNPUWeight); delete hDNPUWeight; hDNPUWeight = nullptr;
  }
*/
//  event.rng.setSize(20);

  // read input data
  cfgop.readData(analysis.datapath);
  for (auto& op : ops_all)
    op->readData(analysis.datapath);

  if (DEBUG) logger.debug("PandaAnalyzer::PandaAnalyzer","Called constructor");
}


PandaAnalyzer::~PandaAnalyzer()
{
}

void PandaAnalyzer::AddGoodLumiRange(int run, int l0, int l1)
{
  auto run_ = goodLumis.find(run);
  if (run_==goodLumis.end()) { // don't know about this run yet
    vector<LumiRange> newLumiList;
    newLumiList.emplace_back(l0,l1);
    goodLumis[run] = newLumiList;
  } else {
    run_->second.emplace_back(l0,l1);
  }
}


bool PandaAnalyzer::PassGoodLumis(int run, int lumi)
{
  auto run_ = goodLumis.find(run);
  if (run_==goodLumis.end()) {
    // matched no run
    if (DEBUG)
      logger.debug("PandaAnalyzer::PassGoodLumis",TString::Format("Failing run=%i",run));
    return false;
  }

  // found the run, now look for a lumi range
  for (auto &range : run_->second) {
    if (range.Contains(lumi)) {
      if (DEBUG)
        logger.debug("PandaAnalyzer::PassGoodLumis",TString::Format("Accepting run=%i, lumi=%i",run,lumi));
      return true;
    }
  }

  // matched no lumi range
  if (DEBUG)
    logger.debug("PandaAnalyzer::PassGoodLumis",TString::Format("Failing run=%i, lumi=%i",run,lumi));
  return false;
}


bool PandaAnalyzer::PassPresel(Selection::Stage stage)
{
  if (selections.size() == 0)
    return true;

  bool pass = false;
  for (auto& s : selections) {
    if (s->anded())
      continue;
    if (DEBUG>1)
      logger.debug("PandaAnalyzer::PassPresel",s->get_name());
    if (s->accept(stage)) {
      pass = true;
      break;
    }
  }

  for (auto& s : selections) {
    if (s->anded()) {
      if (DEBUG>1)
        logger.debug("PandaAnalyzer::PassPresel",s->get_name());
      pass = pass && s->accept(stage);
    }
  }

  return pass;
}



void PandaAnalyzer::Reset()
{
  for (auto& op : ops_all)
    op->reset();

  Analyzer::Reset();
}



void PandaAnalyzer::Terminate()
{
  for (auto& op : ops_all)
    op->terminate();

  Analyzer::Terminate();
}


// run
void PandaAnalyzer::Run()
{

  // INITIALIZE --------------------------------------------------------------------------
  unsigned nZero, nEvents, iE=0;
  setupRun(nZero, nEvents); 

  for (auto& op : ops_all){
    cout << "Printing" << endl;
    op->print();
    op->initialize(registry);
  }
  cout << "Done printing" << endl;

  ProgressReporter pr("PandaAnalyzer::Run",&iE,&nEvents,100);
  TimeReporter& tr = cfgop.cfg.tr;
  tr.TriggerEvent("configuration"); 

  std::cout <<"main run, nEvent=" << nEvents << "\n";
  // EVENTLOOP --------------------------------------------------------------------------
  TH1D* hDTotalMCWeight = new TH1D("MC_weight","MC_weight",7,0,7);

  float summcweight=0;
  for (iE=nZero; iE!=nEvents; ++iE) {
    pr.Report();
    if(iE%10000==0) std::cout << "entry=" << iE << "\n";
    Reset();
    event.getEntry(*tIn,iE);
    tr.TriggerEvent(TString::Format("GetEntry %u",iE));
    gblop->execute();
      if(gt.mcWeight>=0) summcweight++;
      else {summcweight--;}

//    cout << "run:lumi:event=" << gt.runNumber << ":"<< gt.lumiNumber <<":"<<gt.eventNumber << endl;
    if (analysis.isData && !PassGoodLumis(gt.runNumber,gt.lumiNumber))
        continue;
    preselop->execute();
    if (!PassPresel(Selection::sReco))
      continue;
    postselop->execute();
//    if (!PassPresel(Selection::sGen)) // only check gen presel here
//      continue;
    gt.Fill();
    tr.TriggerEvent("fill");
  }

  cout << "summcweight=" << summcweight <<endl; 
  hDTotalMCWeight->SetBinContent(1,summcweight);
  hDTotalMCWeight->Write("hDTotalMCWeight",TObject::kOverwrite);

  pr.Done();

  tr.Summary();
  for (auto& s : selections)
    s->report();

  if (DEBUG) { logger.debug("PandaAnalyzer::Run","Done with entry loop"); }

} // Run()
