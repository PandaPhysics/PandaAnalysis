#ifndef LEPPHOMODS
#define LEPPHOMODS

#include "Module.h"
#include "array"

namespace pa {
  class SimpleLeptonMod : public AnalysisMod {
  public: 
    SimpleLeptonMod(panda::EventAnalysis& event_, 
                    const Config& cfg_,
                    const Utils& utils_,
                    GeneralTree& gt_) : 
      AnalysisMod("simplep", event_, cfg_, utils_, gt_) { }
    virtual ~SimpleLeptonMod () { }

    virtual bool on() { return !analysis.genOnly && !analysis.complicatedLeptons; }
    
  protected:
    virtual void do_init(Registry& registry) {
      registry.publishConst("looseLeps", &looseLeps);
      registry.publishConst("tightLeps", &tightLeps);
      registry.publishConst("lepPdgId", &lepPdgId);
    }
    virtual void do_execute(); 
    virtual void do_reset() { 
      looseLeps.clear();
      tightLeps.clear();
      for (auto& id : lepPdgId)
        id = 0;
    }
    void scaleFactors();

    std::vector<panda::Lepton*> looseLeps, tightLeps;
    std::array<int,4> lepPdgId;
    
  };


  class ComplicatedLeptonMod : public SimpleLeptonMod {
  public: 
    ComplicatedLeptonMod(panda::EventAnalysis& event_, 
                         const Config& cfg_,
                         const Utils& utils_,
                         GeneralTree& gt_) : 
      AnalysisMod("complep", event_, cfg_, utils_, gt_) { }
    ~ComplicatedLeptonMod () { delete rochesterCorrection; }

    virtual bool on() { return !analysis.genOnly && analysis.complicatedLeptons; }
    
  protected:
    void do_readData(TString dirPath);
    void do_execute(); 
    void do_initalize(Registry& registry) {
      SimpleLeptonMod::do_initalize(registry);
      genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
    }
  private:
    RoccoR *rochesterCorrection{nullptr};
    const std::vector<panda::Particle*> *genP{nullptr};
  };

  class InclusiveLeptonMod : public AnalysisMod {
  public: 
    InclusiveLeptonMod(panda::EventAnalysis& event_, 
                       const Config& cfg_,
                       const Utils& utils_,
                       GeneralTree& gt_) : 
      AnalysisMod("inclep", event_, cfg_, utils_, gt_) { }
    virtual ~InclusiveLeptonMod () { }

    virtual bool on() { return !analysis.genOnly && analysis.hbb; }
    
  protected:
    virtual void do_init(Registry& registry) {
      registry.publishConst("inclusiveLeps", &inclusiveLeps);
    }
    virtual void do_execute(); 
    virtual void do_reset() { 
      inclusiveLeps.clear();
    }
  private:
    std::vector<panda::Lepton*> inclusiveLeps;
  };


  class SimplePhotonMod : public AnalysisMod {
  public: 
    SimplePhotonMod(panda::EventAnalysis& event_, 
                    const Config& cfg_,
                    const Utils& utils_,
                    GeneralTree& gt_) : 
      AnalysisMod("simplepho", event_, cfg_, utils_, gt_) { }
    virtual ~SimplePhotonMod () { }

    virtual bool on() { return !analysis.genOnly && !analysis.complicatedPhotons; }
    
  protected:
    virtual void do_init(Registry& registry) {
      registry.publishConst("loosePhos", &loosePhos);
    }
    virtual void do_execute(); 
    virtual void do_reset() { 
      loosePhos.clear();
    }
    void scaleFactors();
    std::vector<panda::Photon*> loosePhos;
  };


  class ComplicatedPhotonMod : public SimplePhotonMod {
  public: 
    ComplicatedPhotonMod(panda::EventAnalysis& event_, 
                         const Config& cfg_,
                         const Utils& utils_,
                         GeneralTree& gt_) : 
      AnalysisMod("comppho", event_, cfg_, utils_, gt_) { }
    ~ComplicatedPhotonMod () { }

    virtual bool on() { return !analysis.genOnly && analysis.complicatedPhotons; }
    
  protected:
    void do_init(Registry& registry) {
      SimplePhotonMod::do_init(registry);
      matchLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps");
    }
    void do_execute(); 
  private:
    const std::vector<panda::Lepton*> *matchLeps; 
    bool pfChargedPhotonMatch(const panda::Photon& photon);
  };


  class TauMod : public AnalysisMod {
  public: 
    TauMod(panda::EventAnalysis& event_, 
           const Config& cfg_,
           const Utils& utils_,
           GeneralTree& gt_) : 
      AnalysisMod("tau", event_, cfg_, utils_, gt_) { }
    ~TauMod () { }

    virtual bool on() { return !analysis.genOnly; }
    
  protected:
    void do_init(Registry& registry) {
      matchLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps");
    }
    void do_execute(); 
  private:
    const std::vector<panda::Lepton*> *matchLeps; 
  };

  class GenLepMod : public AnalysisMod {
  public: 
    GenLepMod(panda::EventAnalysis& event_, 
              const Config& cfg_,
              const Utils& utils_,
              GeneralTree& gt_) : 
      AnalysisMod("genlep", event_, cfg_, utils_, gt_) { }
    ~GenLepMod () { }

    virtual bool on() { return analysis.vbf && !analysis.isData; }
    
  protected:
    void do_execute(); 
    void do_initalize(Registry& registry) {
      genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
    }
  private:
    const std::vector<panda::Particle*> *genP{nullptr};
  };

}

#endif
