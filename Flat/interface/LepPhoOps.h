#ifndef LEPPHOOPS
#define LEPPHOOPS

#include "Operator.h"
#include "AnalyzerUtilities.h"
#include "array"
#include "PandaAnalysis/Utilities/interface/RoccoR.h"
#include "PandaAnalysis/Utilities/interface/RelIso.h"

namespace pa {

  class SimpleLeptonOp : public AnalysisOp {
  public:
    SimpleLeptonOp(panda::Event& event_,
                    Config& cfg_,
                    Utils& utils_,
                    GeneralTree& gt_,
                    int level_=0) :
      AnalysisOp("simplep", event_, cfg_, utils_, gt_, level_),
      looseLeps(std::v_make_shared<panda::Lepton*>()),
      matchLeps(std::v_make_shared<panda::Particle*>()),
      lepPdgId(std::make_shared<std::array<int,4>>()),
      dilep(std::make_shared<TLorentzVector>()) { }
    virtual ~SimpleLeptonOp () { }

    virtual bool on() { return !analysis.genOnly && !analysis.complicatedLeptons; }

  protected:
    virtual void do_init(Registry& registry) {
      registry.publishConst("looseLeps", looseLeps); // sink sheps
      registry.publishConst("matchLeps", matchLeps);
      registry.publishConst("lepPdgId", lepPdgId);
      registry.publishConst("dilep", dilep);
      jesShifts = registry.access<std::vector<JESHandler>>("jesShifts");
//        loosePhos = registry.accessConst<std::vector<panda::Photon*>>("loosePhos");
//        tightPhos = registry.accessConst<std::vector<panda::Photon*>>("tightPhos");

    }
    virtual void do_execute();
    virtual void do_reset() {
      looseLeps->clear();
      matchLeps->clear();
      for (auto& id : *lepPdgId)
        id = 0;
      dilep->SetPtEtaPhiM(0,0,0,0);
    }
    void scaleFactors();

    std::shared_ptr<std::vector<panda::Lepton*>> looseLeps;
    std::shared_ptr<std::vector<panda::Particle*>> matchLeps;
    std::shared_ptr<std::array<int,4>> lepPdgId;
    std::shared_ptr<std::vector<JESHandler>> jesShifts{nullptr};
    std::shared_ptr<TLorentzVector> dilep;
  private:
//    std::shared_ptr<const std::vector<panda::Photon*>> loosePhos{nullptr};
//    std::shared_ptr<const std::vector<panda::Photon*>> tightPhos{nullptr};
  };

  class ComplicatedLeptonOp : public SimpleLeptonOp {
  public:
    ComplicatedLeptonOp(panda::Event& event_,
                         Config& cfg_,
                         Utils& utils_,
                         GeneralTree& gt_,
                         int level_=0) :
      SimpleLeptonOp(event_, cfg_, utils_, gt_, level_) { name = "complep"; }
    virtual ~ComplicatedLeptonOp () { }

    virtual bool on() { return !analysis.genOnly && analysis.complicatedLeptons; }

  protected:
    void do_readData(TString dirPath);
    void do_execute();
    void do_init(Registry& registry) {
      SimpleLeptonOp::do_init(registry);
      if (analysis.hbb) 
        pfCandsMap = registry.access<EtaPhiMap<panda::PFParticle>>("pfCandsMap"); 
    }
  private:
    std::unique_ptr<RoccoR> rochesterCorrection{nullptr};
    std::shared_ptr<EtaPhiMap<panda::PFParticle>> pfCandsMap{nullptr}; 
  };

  class SimplePhotonOp : public AnalysisOp {
  public:
    SimplePhotonOp(panda::Event& event_,
                    Config& cfg_,
                    Utils& utils_,
                    GeneralTree& gt_,
                    int level_=0) :
      AnalysisOp("simplepho", event_, cfg_, utils_, gt_, level_),
      veryLoosePhos(std::v_make_shared<panda::Photon*>()),
      loosePhos(std::v_make_shared<panda::Photon*>()),
      tightPhos(std::v_make_shared<panda::Photon*>()) { }
    virtual ~SimplePhotonOp () { }

    virtual bool on() { return !analysis.genOnly && !analysis.complicatedPhotons; }

  protected:
    virtual void do_init(Registry& registry) {
      registry.publishConst("veryLoosePhos", veryLoosePhos);
      registry.publishConst("loosePhos", loosePhos);
      registry.publishConst("tightPhos", tightPhos);
      matchLeps = registry.accessConst<std::vector<panda::Particle*>>("matchLeps");
    }
    virtual void do_execute();
    virtual void do_reset() {
      veryLoosePhos->clear();
      loosePhos->clear();
      tightPhos->clear();
    }
    void scaleFactors();
    std::shared_ptr<std::vector<panda::Photon*>> veryLoosePhos, loosePhos, tightPhos;
  private:
    std::shared_ptr<const std::vector<panda::Particle*>> matchLeps{nullptr};

  };


  class ComplicatedPhotonOp : public SimplePhotonOp {
  public:
    ComplicatedPhotonOp(panda::Event& event_,
                         Config& cfg_,
                         Utils& utils_,
                         GeneralTree& gt_,
                         int level_=0) :
      SimplePhotonOp(event_, cfg_, utils_, gt_, level_) { name="comppho"; }
    virtual ~ComplicatedPhotonOp () { }

    virtual bool on() { return !analysis.genOnly && analysis.complicatedPhotons; }

  protected:
    void do_init(Registry& registry) {
      SimplePhotonOp::do_init(registry);
//      if (!analysis.darkg)
    }
    void do_execute();
  private:
    bool pfChargedPhotonMatch(const panda::Photon& photon);
    std::shared_ptr<const std::vector<panda::Particle*>> matchLeps{nullptr};
  };

  class TauOp : public AnalysisOp {
  public:
    TauOp(panda::Event& event_,
           Config& cfg_,
           Utils& utils_,
           GeneralTree& gt_,
           int level_=0) :
      AnalysisOp("tau", event_, cfg_, utils_, gt_, level_) { }
    virtual ~TauOp () { }

    virtual bool on() { return !analysis.genOnly; }

  protected:
    void do_init(Registry& registry) {
      matchLeps = registry.accessConst<std::vector<panda::Particle*>>("matchLeps");
    }
    void do_execute();
  private:
    std::shared_ptr<const std::vector<panda::Particle*>> matchLeps{nullptr};
};

  class GenLepOp : public AnalysisOp {
  public:
    GenLepOp(panda::Event& event_,
              Config& cfg_,
              Utils& utils_,
              GeneralTree& gt_,
              int level_=0) :
      AnalysisOp("genlep", event_, cfg_, utils_, gt_, level_) { }
    virtual ~GenLepOp () { }

    virtual bool on() { return analysis.vbf && !analysis.isData; }

  protected:
    void do_execute();
    void do_init(Registry& registry) {
      genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
    }
  private:
    std::shared_ptr<const std::vector<panda::Particle*>> genP{nullptr};
};

}

#endif
