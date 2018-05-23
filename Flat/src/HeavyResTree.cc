// THIS FILE IS AUTOGENERATED //
#include "../interface/HeavyResTree.h"
// STARTCUSTOM INCLUDE

#include <iostream>

//ENDCUSTOM
HeavyResTree::HeavyResTree() {
// STARTCUSTOM CONSTRUCTOR
    for (auto ibeta : ibetas) {
        for (auto N : Ns) {
            for (auto order : orders) {
                ECFParams p;
                p.ibeta = ibeta;
                p.N = N;
                p.order = order;
                ecfParams.push_back(p);
                clf_ECFNs[p] = -1;
            }
        }
    }
// ENDCUSTOM
  Reset();
}
HeavyResTree::~HeavyResTree() {
// STARTCUSTOM DESTRUCTOR
// ENDCUSTOM
}
void HeavyResTree::SetAuxTree(TTree *t) {
// STARTCUSTOM AUX
// ENDCUSTOM
}
void HeavyResTree::Reset() {
// STARTCUSTOM RESET
    for (auto p : ecfParams) {
        clf_ECFNs[p] = -1;
    }
// ENDCUSTOM
  runNumber = 0;
  lumiNumber = 0;
  i_evt = -99;
  npv = 0;
  rho = -99;
  mcWeight = -99;
  sampleType = 0;
  gen_pt = -99;
  gen_eta = -99;
  gen_phi = -99;
  gen_size = -99;
  gen_pdgid = 0;
  clf_Tau32 = -99;
  clf_Tau21 = -99;
  clf_Tau32SD = -99;
  clf_Tau21SD = -99;
  clf_MSD = -99;
  clf_MSD_corr = -99;
  clf_Pt = -99;
  clf_Phi = -99;
  clf_Eta = -99;
  clf_M = -99;
  clf_MaxCSV = -99;
  clf_IsMatched = 0;
  clf_HTTFRec = -99;
}
void HeavyResTree::WriteTree(TTree *t) {
  treePtr = t;
// STARTCUSTOM WRITE
    for (auto p : ecfParams) {
        TString ecfn(makeECFString(p));
        Book("clf_"+ecfn,&(clf_ECFNs[p]),"clf_"+ecfn+"/F");
    }
// ENDCUSTOM
  Book("runNumber",&runNumber,"runNumber/I");
  Book("lumiNumber",&lumiNumber,"lumiNumber/I");
  Book("i_evt",&i_evt,"i_evt/l");
  Book("npv",&npv,"npv/I");
  Book("rho",&rho,"rho/F");
  Book("mcWeight",&mcWeight,"mcWeight/F");
  Book("sampleType",&sampleType,"sampleType/I");
  Book("gen_pt",&gen_pt,"gen_pt/F");
  Book("gen_eta",&gen_eta,"gen_eta/F");
  Book("gen_phi",&gen_phi,"gen_phi/F");
  Book("gen_size",&gen_size,"gen_size/F");
  Book("gen_pdgid",&gen_pdgid,"gen_pdgid/I");
  Book("clf_Tau32",&clf_Tau32,"clf_Tau32/F");
  Book("clf_Tau21",&clf_Tau21,"clf_Tau21/F");
  Book("clf_Tau32SD",&clf_Tau32SD,"clf_Tau32SD/F");
  Book("clf_Tau21SD",&clf_Tau21SD,"clf_Tau21SD/F");
  Book("clf_MSD",&clf_MSD,"clf_MSD/F");
  Book("clf_MSD_corr",&clf_MSD_corr,"clf_MSD_corr/F");
  Book("clf_Pt",&clf_Pt,"clf_Pt/F");
  Book("clf_Phi",&clf_Phi,"clf_Phi/F");
  Book("clf_Eta",&clf_Eta,"clf_Eta/F");
  Book("clf_M",&clf_M,"clf_M/F");
  Book("clf_MaxCSV",&clf_MaxCSV,"clf_MaxCSV/F");
  Book("clf_IsMatched",&clf_IsMatched,"clf_IsMatched/I");
  Book("clf_HTTFRec",&clf_HTTFRec,"clf_HTTFRec/F");
}