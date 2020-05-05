#include "PandaCore/Tools/interface/Common.h"
#include "AnalyzerUtilities.h"
#include "TString.h"
#include "GeneralTree.h"
#include <algorithm>

#ifndef PANDA_SELECTION
#define PANDA_SELECTION

#define __ACCPFUNC(x) [=](const GeneralTree* gt) { return x; }

namespace pa {
  template <typename T>
  class BaseSelection {
  public:
    enum Stage {
      sGen, sReco
    };

    BaseSelection(Stage stage_, TString n=""): stage(stage_), name(n) { }
    virtual ~BaseSelection() { }
    
    virtual void report() const final { 
      logger.debug("Selection::" + name, Form("Accepted %i/%i events", nPassed, nTotal)); 
    }
    virtual void set_gt(const T* gt_) final { gt = gt_; }
    // if called at a different stage, just return true
    virtual bool accept(Stage stage_) final { 
      if (stage_ != stage)
        return true;
      bool good = do_accept(); 
      ++nTotal; if (good) ++nPassed;
      return good;
    };
    virtual bool anded() const final { return is_anded; }
    virtual TString get_name() const final { return name; }
  protected:
    virtual bool do_accept() const = 0;
    Stage stage;
    int nTotal{0}, nPassed{0};
    const T* gt{nullptr};
    TString name;
    bool is_anded{false};
  };

  typedef BaseSelection<GeneralTree> Selection;

  // Lambda class that can be bound at runtime or subclassed (see below)
  class LambdaSel : public Selection {
  public:
    typedef std::function<bool(const GeneralTree*)> accept_func;
    LambdaSel(Stage stage_, TString n, accept_func f_, bool anded = false): 
      Selection(stage_,n), 
      f(f_) { is_anded = anded; }
    ~LambdaSel() { }

  protected:
    virtual bool do_accept() const final { return f(gt); }
  private:
    accept_func f;
  };

  class CategorySel : public LambdaSel {
  public:
  //  CategorySel(bool tightOnly=false):
  //    LambdaSel(Selection::sReco, "Category",
  //              __ACCPFUNC((gt->category>0) || (gt->category!=0 && !tightOnly))) { }
  };

  class VBFGamma : public Selection {
  public:
  VBFGamma(): Selection(Selection::sReco, "vbfgamma") { }
  protected:
    virtual bool do_accept() const;
  };
}

#endif
