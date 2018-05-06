#ifndef THEORYMODS
#define THEORYMODS

#include "Module.h"

namespace pa {
  class GenStudyEWKMod : public AnalysisMod {
  public: 
    GenStudyEWKMod(panda::EventAnalysis& event_, 
                    const Config& cfg_,
                    const Utils& utils_,
                    GeneralTree& gt_) : 
      AnalysisMod("genstudyewk", event_, cfg_, utils_, gt_) { }
    virtual ~GenStudyEWKMod () { }

    virtual bool on() { return !analysis.isData && (analysis.complicatedLeptons || analysis.complicatedPhotons); }
    
  protected:
    virtual void do_init(Registry& registry) {
      genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
      looseLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps");
      loosePhos = registry.accessConst<std::vector<panda::Lepton*>>("loosePhos");
    }
    virtual void do_execute(); 
  private:
    const std::vector<panda::Particle*> *genP{nullptr};
    const std::vector<panda::Lepton*> *looseLeps{nullptr}; 
    const std::vector<panda::Lepton*> *loosePhos{nullptr}; 
    
  };


  class QCDUncMod : public AnalysisMod {
  public: 
    QCDUncMod(panda::EventAnalysis& event_, 
                    const Config& cfg_,
                    const Utils& utils_,
                    GeneralTree& gt_) : 
      AnalysisMod("qcdunc", event_, cfg_, utils_, gt_) { }
    virtual ~QCDUncMod () { }

    virtual bool on() { return !analysis.isData; }
    
  protected:
    void do_execute(); 
  };


  class SignalGenMod : public AnalysisMod {
  public: 
    SignalGenMod(panda::EventAnalysis& event_, 
                    const Config& cfg_,
                    const Utils& utils_,
                    GeneralTree& gt_) : 
      AnalysisMod("signalweight", event_, cfg_, utils_, gt_) { }
    virtual ~SignalGenMod () { }

    virtual bool on() { return !analysis.isData && analysis.processType==kSignal; }
    
  protected:
    void do_execute(); 
    void do_init(Registry& registry) {
      wIDs = registry.accessConst<std::vector<TString>>("wIDs"); 
      genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
    }
  private:
    const std::vector<TString> *wIDs {nullptr};
    const std::vector<panda::Particle*> *genP{nullptR}
  };


  class HFCountingMod : public AnalysisMod {
  public: 
    HFCountingMod(panda::EventAnalysis& event_, 
                    const Config& cfg_,
                    const Utils& utils_,
                    GeneralTree& gt_) : 
      AnalysisMod("hfcounting", event_, cfg_, utils_, gt_) { }
    virtual ~HFCountingMod () { }

    virtual bool on() { return !analysis.isData; }
    
  protected:
    void do_execute(); 
    void do_init(Registry& registry) {
      genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
    }
  private:
    const std::vector<panda::Particle*> *genP{nullptR}
  };

  class KFactorMod : public AnalysisMod {
  public: 
    KFactorMod(panda::EventAnalysis& event_, 
                    const Config& cfg_,
                    const Utils& utils_,
                    GeneralTree& gt_) : 
      AnalysisMod("kfactor", event_, cfg_, utils_, gt_) { }
    virtual ~KFactorMod () { }

    virtual bool on() { return !analysis.isData; }
    
  protected:
    void do_execute(); 
    void do_init(Registry& registry) {
      genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
      matchLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps");
    }
  private:
    void toppt(); 
    void vpt(); 

    const std::vector<panda::Particle*> *genP{nullptr}
    const std::vector<panda::Lepton*> *matchLeps{nullptr}; 
  };
}

#endif
