#ifndef JETSOPS
#define JETSOPS

#include "Operator.h"
#include <numeric>

namespace pa {

  class JetFlavorOp : public AnalysisOp {
  public:
    JetFlavorOp(panda::Event& event_,
                 Config& cfg_,
                 Utils& utils_,
                 GeneralTree& gt_,
                 int level_=0) :
      AnalysisOp("jetflavor", event_, cfg_, utils_, gt_, level_) { }
    virtual ~JetFlavorOp () {}

    bool on() { return !analysis.isData && analysis.jetFlavorJets; }
  protected:
    void do_init(Registry& registry) {
      jesShifts = registry.access<std::vector<JESHandler>>("jesShifts");
      if (!analysis.isData)
        genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
    }
    void do_execute();
  private:
    std::shared_ptr<std::vector<JESHandler>> jesShifts{nullptr};
    std::shared_ptr<const std::vector<panda::Particle*>> genP{nullptr};

    void clusteredFlavor(JetWrapper&);
  };


  class VBFSystemOp : public AnalysisOp {
  public:
    VBFSystemOp(panda::Event& event_,
                 Config& cfg_,
                 Utils& utils_,
                 GeneralTree& gt_,
                 int level_=0) :
      AnalysisOp("vbfsystem", event_, cfg_, utils_, gt_, level_) { }
    virtual ~VBFSystemOp () {}

    bool on() { return analysis.vbf; }
  protected:
    void do_init(Registry& registry) {
      currentJES = registry.access<JESHandler*>("currentJES");
    }
    void do_execute();
  private:
    std::shared_ptr<JESHandler*> currentJES{nullptr};
  };

  class BaseJetOp : public AnalysisOp {
  public:
    BaseJetOp(TString name,
               panda::Event& event_,
               Config& cfg_,
               Utils& utils_,
               GeneralTree& gt_,
               int level_=0) :
      AnalysisOp(name, event_, cfg_, utils_, gt_, level_),
      recalcJER(false) {
        if (analysis.year == 2016) {
          jecV = "V4"; jecReco = "23Sep2016";
          campaign = "Summer16";
          jerV = "Spring16_25nsV10";
          eraGroups = {"BCD","EF","G","H"};
          spacer = "";
          if (analysis.useDeepCSV) { 
            csvL = 0.2219; csvM = 0.6324; 
          } else { 
            csvL = 0.5426; csvM = 0.8484; 
          }
        } 
	else if (analysis.year == 2017) {
          jecV = "V32"; jecReco = "17Nov2017";
          campaign = "Fall17";
          jerV = "Fall17_25nsV1";
          eraGroups = {"B","C","DE","F"};
          spacer = "_";
          if (analysis.useDeepCSV) { 
            csvL = 0.1522; csvM = 0.4941; 
          } 
          else if(analysis.useDeepJet){
            csvL = 0.0521; csvM = 0.3033;
          }
          else { 
            csvL = 0.5803; csvM = 0.8838; 
          }
	}
	else if (analysis.year == 2018) {
          jecV = "V19"; jecReco = "Autumn18";
          campaign = "Winter19";
          jerV = "Autumn18_V7";
          eraGroups = {"A","B","C","D"};
          spacer = "_";
          if (analysis.useDeepCSV) { 
            csvL = 0.1241; csvM = 0.4184; 
          } 
          else if(analysis.useDeepJet){
            csvL = 0.0494; csvM = 0.2770;
          }
          else { 
            csvL = 0.5803; csvM = 0.8838; 
          }
        }
      }
    virtual ~BaseJetOp() { }
    bool csvLoose(float csv) { return csv > csvL; }
    bool csvMed(float csv) { return csv > csvM; }
    TLorentzVector met_corr, met_jer_corr;
  protected:
    bool recalcJER;
    virtual void do_execute() = 0;
    virtual void do_readData(TString path);
    JetWrapper shiftJet(const panda::Jet& jet, shiftjes shift, bool smear=false);
    JetWrapper shiftJet(const panda::FatJet& jet, shiftjes shift, bool smear=false);

    void shiftMET(const panda::RecoMet& met, TLorentzVector& v, shiftjes shift, bool metjer=false);
    void lowptshift(const panda::CorrT1METJet& jet, bool smear=false);
 
    std::map<TString,std::unique_ptr<FactorizedJetCorrector>> scales; // era/MC -> scale
    std::map<TString,std::vector<std::shared_ptr<JetCorrectionUncertainty>>> scaleUncs; // era/MC -> (src -> unc)
    std::unique_ptr<JERReader> jer{nullptr}; //!< fatjet jet energy resolution reader

    std::vector<std::shared_ptr<JetCorrectionUncertainty>> *scaleUnc  {nullptr}; // src -> unc
    FactorizedJetCorrector *scale{nullptr};

    TString jecV, jecReco, jetType, campaign, spacer, jerV;
    std::vector<TString> eraGroups;
    float csvL, csvM;

  private:
    void setScaleUnc(TString, TString);
  };

  class JetOp : public BaseJetOp {
  public:
    JetOp(panda::Event& event_,
           Config& cfg_,
           Utils& utils_,
           GeneralTree& gt_,
           int level_=0) :
      BaseJetOp("jet", event_, cfg_, utils_, gt_, level_),
      currentJet(std::make_shared<JetWrapper*>(nullptr)),
      currentJES(std::make_shared<JESHandler*>(nullptr)) {
	ak4Jets = &(event.Jet);
        recalcJER = analysis.rerunJER; 
        vbf = addSubOp<VBFSystemOp>();
	jetType = "AK4PFchs";	  
    }
    virtual ~JetOp () { }

    virtual bool on() { return !analysis.genOnly; }

  protected:
    void do_init(Registry& registry) {
      registry.publish("currentJet", currentJet);
      registry.publish("currentJES", currentJES);
      jesShifts = registry.access<std::vector<JESHandler>>("jesShifts");
      matchLeps = registry.accessConst<std::vector<panda::Particle*>>("matchLeps");
      matchPhos = registry.accessConst<std::vector<panda::Photon*>>("tightPhos");
      matchVeryLoosePhos = registry.accessConst<std::vector<panda::Photon*>>("veryLoosePhos");
    }
    void do_execute();

  private:
    VBFSystemOp *vbf{nullptr};

    std::shared_ptr<std::vector<JESHandler>> jesShifts{nullptr};

    std::shared_ptr<const std::vector<panda::Particle*>> matchLeps{nullptr};
    std::shared_ptr<const std::vector<panda::Photon*>> matchPhos{nullptr};
    std::shared_ptr<const std::vector<panda::Photon*>> matchVeryLoosePhos{nullptr};

    panda::JetCollection *ak4Jets{nullptr};

    std::shared_ptr<JetWrapper*> currentJet; // shared ptr to a bare address
    std::shared_ptr<JESHandler*> currentJES;

    void setupJES();
    void varyJES();
  };
}

#endif
