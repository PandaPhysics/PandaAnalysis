#ifndef DEEPOPS
#define DEEPOPS

#include "Operator.h"
#include "set"
#include "PandaAnalysis/Utilities/interface/EnergyCorrelations.h"
#include "fastjet/contrib/Njettiness.hh"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "TMVA/Reader.h"
#include <algorithm>

namespace pa {

  class TFInferOp : public AnalysisOp {
  public:
    TFInferOp(TString name_,
               panda::EventAnalysis& event_,
               Config& cfg_,
               Utils& utils_,
               GeneralTree& gt_,
               int level_=0) :
      AnalysisOp(name_, event_, cfg_, utils_,  gt_, level_),
      p_inputs(std::v_make_shared<float>()),
      p_outputs(std::v_make_shared<float>()),
      inputs(*p_inputs),
      outputs(*p_outputs),
      transpose(false),
      t_i(1)    { }
    ~TFInferOp() { }

  protected:
    virtual void do_init(Registry& registry) {
      // just in case you want to set these in other ops
      registry.publishConst(name+"_outputs", p_outputs);
      registry.publish(name+"_inputs", p_inputs);
    }
    virtual void do_readData(TString dirPath) = 0; // must call build(weightpath)
    virtual void do_execute() = 0; // must call eval()
    virtual void do_reset() {
      fill(inputs.begin(), inputs.end(), -99);
      fill(outputs.begin(), outputs.end(), -99);
    }

    virtual void build(TString weightpath) final;
    virtual void eval() final;

    TString inputName{0};
    std::vector<std::string> outputNames;
    int n_inputs{0}, n_outputs{0};
    std::shared_ptr<std::vector<float>> p_inputs, p_outputs; // keep this around for publication
    std::vector<float>& inputs;
    std::vector<float>& outputs;
    bool transpose;

  private:
    std::unique_ptr<tensorflow::GraphDef> graph{nullptr}; 
    std::unique_ptr<tensorflow::Session> sess{nullptr};
    tensorflow::NamedTensorList t_i;
  };


  class BRegDeepOp : public TFInferOp {
  public:
    BRegDeepOp(panda::EventAnalysis& event_,
                Config& cfg_,
                Utils& utils_,
                GeneralTree& gt_,
                int level_=0) :
      TFInferOp("bregdeep", event_, cfg_, utils_, gt_, level_) {
      //n_inputs = 43;
      n_inputs = 51;
      if (analysis.bjetDeepReg_withPFs)
	n_inputs = 76;
      if (analysis.bjetDeepReg_withPFs_3)
	n_inputs = 66;
      if (analysis.bjetDeepReg_withPFs_10)
	n_inputs = 101;
      if (analysis.bjetDeepReg_withD0)
	n_inputs = 86;
      n_outputs = 3;
      if (analysis.bjetDeepReg_EtaPhi)
	n_outputs = 5;

      //inputName = "ffwd_inp";
      inputName = "input";
      outputNames.reserve(n_outputs);
      for (int i = 0; i != n_outputs; ++i)
        outputNames.push_back(Form("output_%i/BiasAdd", i));
      //outputNames.push_back("ffwd_out/BiasAdd");
    }

    virtual bool on() { return !analysis.genOnly && analysis.hbb && analysis.bjetDeepReg; }
  protected:
    void do_readData(TString dirPath) {
      TString modelfile = dirPath+"/trainings/breg_graph.pb";
      //downloadData("http://t3serv001.mit.edu/~snarayan/pandadata/trainings/breg/v2/quantiles/graph.pb",
      //downloadData("http://t3serv001.mit.edu/~snarayan/pandadata/trainings/breg_training_2017_updated.pb",

      if (analysis.year == 2016)
	downloadData("http://t3serv001.mit.edu/~snarayan/pandadata/trainings/sidbreg_v0/graph.pb",
		     modelfile, true);
      else if (analysis.year == 2017){
	if (analysis.bjetDeepReg_withPFs)
	  downloadData("http://t3serv001.mit.edu/~bmaier/figs/hbb/regression/graph_v2017_withPFs.pb",
		     modelfile, true);
	else if (analysis.bjetDeepReg_withPFs_3)
	  downloadData("http://t3serv001.mit.edu/~bmaier/figs/hbb/regression/graph_v2017_withPFs_3.pb",
		     modelfile, true);
	else if (analysis.bjetDeepReg_withD0)
	  downloadData("http://t3serv001.mit.edu/~bmaier/figs/hbb/regression/graph_v2017_withD0.pb",
		     modelfile, true);
	else{
	  downloadData("http://t3serv001.mit.edu/~bmaier/figs/hbb/regression/graph_v2017_2017.pb",
		       modelfile, true);
	  //downloadData("http://t3serv001.mit.edu/~snarayan/pandadata/trainings/sidbreg_v0/graph.pb",
	  //	     modelfile, true);
	}
      }
      else{
	if (analysis.bjetDeepReg_EtaPhi)
	  downloadData("http://t3serv001.mit.edu/~bmaier/figs/hbb/regression/graph_v2018_withPFs_EtaPhi.pb",
		       modelfile, true);
	else if (analysis.bjetDeepReg_withPFs){
	  downloadData("http://t3serv001.mit.edu/~bmaier/figs/hbb/regression/graph_v2018_withPFs_defaultfix_0.pb",
	  	       modelfile, true);
	}
	else if (analysis.bjetDeepReg_withD0)
	  downloadData("http://t3serv001.mit.edu/~bmaier/figs/hbb/regression/graph_v2018_withD0_defaultfix_0.pb",
		       modelfile, true);
	else {
	  downloadData("http://t3serv001.mit.edu/~bmaier/figs/hbb/regression/graph_v2018_defaultfix_0.pb",
	  	     modelfile, true);
	  
	}
      }
      build(modelfile);
    }
    virtual void do_init(Registry& registry) {
      TFInferOp::do_init(registry);
      currentJet = registry.access<JetWrapper*>("higgsDaughterJet");
    }
    void do_execute();
  private:
    std::shared_ptr<JetWrapper*> currentJet{nullptr};
  };

  /* if someone wants to implement inferences for other VH classifiers,
   * please do not just copy-paste the code for this op. we should create 
   * a base classifier op and subclass it for each channel
   */
  class ZvvHClassOp : public TFInferOp {
  public:
    ZvvHClassOp(panda::EventAnalysis& event_,
                Config& cfg_,
                Utils& utils_,
                GeneralTree& gt_,
                int level_=0) :
      TFInferOp("zvvhclass", event_, cfg_, utils_, gt_, level_) {
      n_inputs = 12;
      n_outputs = 1;
      inputName = "input";
      outputNames = {"output/Softmax"};
      transpose = true; 
    }

    virtual bool on() { return !analysis.genOnly && analysis.hbb; }
  protected:
    void do_readData(TString dirPath) {
      TString modelfile = dirPath+"/trainings/zh_graph.pb";
      downloadData("http://t3serv001.mit.edu/~snarayan/pandadata/trainings/sidzh_v2/graph.pb",
                   modelfile, true);
      build(modelfile);
    }
    virtual void do_init(Registry& registry) {
      TFInferOp::do_init(registry);
    }
    void do_execute();
  };

  class PrepareLSTMOp : public AnalysisOp {
  public:
    PrepareLSTMOp(panda::EventAnalysis& event_,
	    Config& cfg_,
            Utils& utils_,
            GeneralTree& gt_,
            int level_=0) :
    AnalysisOp("lstm", event_, cfg_, utils_, gt_, level_) { 
      jetDef.reset(new fastjet::JetDefinition(fastjet::antikt_algorithm,0.4));
    };

    virtual ~PrepareLSTMOp() { }
    bool on() { return analysis.hbb; }
  protected:
    void do_init(Registry& registry) {
      currentJet = registry.access<JetWrapper*>("higgsDaughterJet");
      fOut = registry.access<TFile>("fOut");
      incrementAux(false);
    }
    void do_execute();
    void do_reset() {
      recoJetInfo.reset();
    }
    void do_terminate() {
      incrementAux(true);
    }
    void incrementAux(bool close = false);
  private:
    std::shared_ptr<TFile> fOut{nullptr};
    std::shared_ptr<JetWrapper*> currentJet{nullptr};
    std::unique_ptr<fastjet::JetDefinition> jetDef{nullptr};
    std::unique_ptr<TFile>fAux{nullptr};
    std::unique_ptr<TTree>tAux{nullptr};
    int auxCounter{0};
    RecoJetInfo recoJetInfo;
  };

  template <typename GENP>
  class DeepGenOp : public AnalysisOp {
  public:
    DeepGenOp(panda::EventAnalysis& event_,
               Config& cfg_,
               Utils& utils_,
               GeneralTree& gt_,
               int level_=0);
    virtual ~DeepGenOp () { }

    virtual bool on() { return !analysis.isData && analysis.deepGen; }

  protected:
    void do_init(Registry& registry) {
      genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
      fOut = registry.access<TFile>("fOut");
      incrementAux(false);
    }
    void do_execute();
    void do_reset() {
      genJetInfo.reset();
      if (grid != nullptr)
        grid->clear();
    }
    void do_terminate() {
      incrementAux(true);
    }
    void countGenPartons(std::unordered_set<const GENP*>&);
    void incrementAux(bool close = false);
  private:
    std::shared_ptr<const std::vector<panda::Particle*>> genP{nullptr};
    std::shared_ptr<TFile> fOut{nullptr};

    std::unique_ptr<fastjet::JetDefinition> jetDef{nullptr};
    std::unique_ptr<fastjet::contrib::Njettiness>tauN{nullptr};
    std::unique_ptr<ECFCalculator>ecfcalc{nullptr};
    std::unique_ptr<ParticleGridder>grid{nullptr};
    GenJetInfo genJetInfo;

    std::unique_ptr<TFile>fAux{nullptr};
    std::unique_ptr<TTree>tAux{nullptr};
    int auxCounter{0};
  };

  typedef DeepGenOp<panda::GenParticle> DeepPGenOp;
  typedef DeepGenOp<panda::UnpackedGenParticle> DeepUGenOp;


  class BRegBDTOp : public AnalysisOp {
  public:
    BRegBDTOp(panda::EventAnalysis& event_,
                Config& cfg_,
                Utils& utils_,
                GeneralTree& gt_,
                int level_=0) :
      AnalysisOp("bregbdt", event_, cfg_, utils_, gt_, level_) { }
    ~BRegBDTOp() { }

    virtual bool on() { return !analysis.genOnly && analysis.hbb && analysis.bjetBDTReg; }
  protected:
    void do_readData(TString dirPath);
    virtual void do_init(Registry& registry) {
      currentJet = registry.access<JetWrapper*>("higgsDaughterJet");
    }
    void do_execute();
  private:
    std::shared_ptr<JetWrapper*> currentJet{nullptr};
    std::unique_ptr<TMVA::Reader> bjetregReader;
    std::unique_ptr<float[]> bjetreg_vars; 
  };

}


#endif
