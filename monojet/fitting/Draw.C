#include <iostream>
#include <TH1.h>
#include <TTree.h>
#include <math.h>

void Draw(){


  TString file1 = "fittingForest_singlemuon.root";
  TString file2 = "new/fittingForest_singlemuon.root";

  TFile *f1 = new TFile(file1);
  TFile *f2 = new TFile(file2);

//  TString hname[6]={"Zmm_zll","Zee_zll","Wen_zll","Wmn_zll","Zee_qcd","Zmm_qcd"};

  f2->cd();
  TH1F *h1=(TH1F*)f2->Get("hGJLOMonoX");
  TH1F *h2=(TH1F*)f2->Get("hGJNLOMonoX");
  h1->SetDirectory(0);
  h2->SetDirectory(0);
  f1->cd();
  TH1F *h3=(TH1F*)f1->Get("hGJLOMonoX");
  TH1F *h4=(TH1F*)f1->Get("hGJLO_DRMonoX");
  h3->SetDirectory(0);
  h4->SetDirectory(0); 

//  h1->Rebin(2);
//  h2->Rebin(2);
//  h3->Rebin(2); 
//  h4->Rebin(2);
 

  h1->Scale(1.156); //scale lumi
  h2->Scale(1.156);

  TPaveText *ll = new TPaveText(0.10, 0.90, 0.95, 0.95, "NDC");
  ll->SetTextSize(0.03);
  ll->SetTextFont(42);
//  ll->SetFillColor(0);
  ll->SetBorderSize(0);
  ll->SetMargin(0.01);
  ll->SetTextAlign(12); // align left
  TString text = "CMS Preliminary";
  ll->AddText(0.01,0.5,text);
  text = "         #sqrt{s} =  13TeV" ;
  //text = "#sqrt{s} = 13 TeV, L = 14.77 fb^{-1}" ;
   ll->AddText(0.65, 0.6, text);

  TCanvas *c1 = new TCanvas("c1","c1",600,800);
  c1->cd();
//  c1->SetLogy(1);
  c1->SetTicks(1,1);


  TH1F *hd = NULL, *hn=NULL,*hn1=NULL,*hn2=NULL;
  hn = (TH1F*)h4->Clone();
  hd = (TH1F*)h3->Clone();
  hn1 = (TH1F*)h1->Clone();
  hn2 = (TH1F*)h2->Clone();
  hn->Divide(hd);
  hn1->Divide(hd);
  hn2->Divide(hd);
  

  hn->SetStats(kFALSE);
  hd->SetStats(kFALSE); 
  hn1->SetStats(kFALSE);
  hn2->SetStats(kFALSE);

//  h1->SetMinimum(0.1);
  h1->SetFillStyle(0);
  h2->SetFillStyle(0);
  h3->SetFillStyle(0);
  h4->SetFillStyle(0); 
  h2->SetLineWidth(2);
  h2->SetLineColor(2);
  h1->SetLineWidth(2);
  h1->SetLineColor(3);
  h3->SetLineWidth(2);
  h4->SetLineWidth(2);
  h4->SetLineColor(4);
  h2->GetYaxis()->SetLabelSize(0.03);
  h2->GetXaxis()->SetLabelSize(0.03);  
  h2->GetYaxis()->SetTitleSize(0.05);
  h2->GetXaxis()->SetTitleSize(0.05);
  h2->SetTitle("");
  h2->Draw("HIST");
  h1->Draw("HIST same");
  h3->Draw("HIST same");
  h4->Draw("HIST same");

  TLegend *legend = new TLegend(0.6,0.70,0.85,0.87);
  legend->AddEntry(h1,"2016 LO","l");
  legend->AddEntry(h2,"2016 NLO","l");
  legend->AddEntry(h3,"2017 LO","l");
  legend->AddEntry(h4,"2017 DR-LO","l");

  legend->Draw();
  ll->Draw("same");
  double canvasratio = 0.3;
  c1->SetBottomMargin(canvasratio + (1-canvasratio)*c1->GetBottomMargin()-canvasratio*c1->GetTopMargin());

  canvasratio = 0.16;
  TPad *ratioPad = new TPad("BottomPad","",0,0,1,1);
  ratioPad->SetTopMargin((1-canvasratio) - (1-canvasratio)*ratioPad->GetBottomMargin()+canvasratio*ratioPad->GetTopMargin());
  ratioPad->SetFillStyle(4000);
//  ratioPad->SetFillColor(4000);
  ratioPad->SetFrameFillColor(4000);
  ratioPad->SetFrameFillStyle(4000);
  ratioPad->SetFrameBorderMode(0);
  ratioPad->SetTicks(1,1);
  hn->GetYaxis()->SetLabelSize(0.03);
  hn->GetXaxis()->SetLabelSize(0.03);
  hn->GetYaxis()->SetTitleSize(0.05);
  hn->GetXaxis()->SetTitleSize(0.05);

  ratioPad->Draw();
  ratioPad->cd();

  hn->SetMarkerStyle(20);
  hn->SetMarkerSize(0.8);
  hn->SetMarkerColor(4);
  hn1->SetMarkerStyle(20);
  hn1->SetMarkerSize(0.8);
  hn1->SetMarkerColor(3);
  hn2->SetMarkerStyle(20);
  hn2->SetMarkerSize(0.8);
  hn2->SetMarkerColor(2);

  hn->SetMinimum(0.5);
  hn->SetMaximum(1.5);
  hn->SetTitle("");
  hn->Draw("Psame");
  hn1->Draw("Psame");
  hn2->Draw("Psame");

  c1->Update();
  c1->SaveAs("plot/kfactor.png");
  f1->Close();
  f2->Close();

}
