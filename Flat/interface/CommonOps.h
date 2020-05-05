#ifndef COMMONOPS
#define COMMONOPS

#include "Operator.h"
#include "AnalyzerUtilities.h"
#include "PandaAnalysis/Utilities/interface/EtaPhiMap.h"

namespace pa {

  class RecoilOp : public AnalysisOp {
  public:
    RecoilOp(panda::Event& event_,
              Config& cfg_,
              Utils& utils_,
              GeneralTree& gt_,
              int level_=0) :
      AnalysisOp("recoil", event_, cfg_, utils_, gt_, level_) { }
    virtual ~RecoilOp () { }

    bool on() { return !analysis.genOnly && analysis.recoil; }

  protected:
    void do_init(Registry& registry) {
      looseLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps");
      cout <<"Recoil init" << endl;
      loosePhos = registry.accessConst<std::vector<panda::Photon*>>("loosePhos");
      lepPdgId = registry.accessConst<std::array<int,4>>("lepPdgId");
      jesShifts = registry.access<std::vector<JESHandler>>("jesShifts");
    }
    void do_execute();

  private:
    std::shared_ptr<const std::vector<panda::Lepton*>> looseLeps{nullptr};
    std::shared_ptr<const std::vector<panda::Photon*>> loosePhos{nullptr};
    std::shared_ptr<const std::array<int,4>> lepPdgId {nullptr};
    std::shared_ptr<std::vector<JESHandler>> jesShifts{nullptr};
  };

  class TriggerOp : public AnalysisOp {
  public:
    TriggerOp(panda::Event& event_,
               Config& cfg_,
               Utils& utils_,
               GeneralTree& gt_,
               int level_=0) :
      AnalysisOp("trigger", event_, cfg_, utils_, gt_, level_),
      triggerHandlers(kNTrig) { }
    virtual ~TriggerOp () { }

    bool on() { return !analysis.genOnly && (analysis.isData || analysis.mcTriggers); }

  protected:
    void do_init(Registry& registry);
    void do_execute();

  private:
//    void checkEle32();
    std::vector<TriggerHandler> triggerHandlers;
  };

  class TriggerEffOp : public AnalysisOp {
  public:
    TriggerEffOp(panda::Event& event_,
               Config& cfg_,
               Utils& utils_,
               GeneralTree& gt_,
               int level_=0) :
      AnalysisOp("triggereff", event_, cfg_, utils_, gt_, level_) { }
    virtual ~TriggerEffOp () { }

    bool on() { return !analysis.genOnly && !analysis.isData; }

  protected:
    void do_init(Registry& registry) {
      looseLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps");
    }
    void do_execute();

  private:
    std::shared_ptr<const std::vector<panda::Lepton*>> looseLeps{nullptr};
};

  class GlobalOp : public AnalysisOp {
  public:
    GlobalOp(panda::Event& event_,
               Config& cfg_,
               Utils& utils_,
               GeneralTree& gt_,
               int level_=0) :
            AnalysisOp("global", event_, cfg_, utils_, gt_, level_) { }

    virtual ~GlobalOp () { }

  protected:
    void do_init(Registry& registry) {
    }
    void do_execute();
    void do_reset() {
    }
  };
}

#endif
