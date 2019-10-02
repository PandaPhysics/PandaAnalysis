#ifndef JCORROPS
#define JCORROPS

#include "Operator.h"
#include "AnalyzerUtilities.h"
#include "JetsOps.h"


namespace pa {
  class JetCorrOp : public BaseJetOp {
  public: 
    JetCorrOp(panda::EventAnalysis& event_, 
               Config& cfg_,
               Utils& utils_,
               GeneralTree& gt_,
               int level_=0) : 
      BaseJetOp("jetcorr", event_, cfg_, utils_, gt_, level_),
      jesShifts(std::v_make_shared<JESHandler>(jes2i(shiftjes::N))) {
        JESLOOP {
          (*jesShifts)[shift].shift_idx = shift;
	  if (analysis.puppiJets)
	    jetType = "AK4PFPuppi";
	  else
	    jetType = "AK4PFchs";
        }
      }
    virtual ~JetCorrOp () { }
    virtual bool on() { return !analysis.genOnly; }
    
  protected:
    void do_init(Registry& registry) {
      registry.publish("jesShifts", jesShifts);
      auto dummy = registry.access<std::vector<JESHandler>>("jesShifts");
    }
    void do_execute();
    void do_reset() {
      for (auto& s : *jesShifts){
        s.clear();
      }
    }    
  private:
    std::shared_ptr<std::vector<JESHandler>> jesShifts;
  };
}

#endif
