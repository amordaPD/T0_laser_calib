#include <iostream>
using namespace std;

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TGraph.h>

void Analyse_ScanTune(){
  ///////////////////////////////////////////////
  /////// channel 4 behaves quite well  /////////
  ///////////////////////////////////////////////
  
  TFile *f_77 = new TFile("fit_results/runScan-T77-pos0_large_scale_out_FR.root");
  TFile *f_82 = new TFile("fit_results/runScan-T82-pos0_large_scale_out_FR.root");
  TFile *f_84 = new TFile("fit_results/runScan-T84-pos0_large_scale_out_FR.root");
  TFile *f_86 = new TFile("fit_results/runScan-T86-pos0_large_scale_out_FR.root");
  TFile *f_88 = new TFile("fit_results/runScan-T88-pos0_large_scale_out_FR.root");
  TFile *f_90 = new TFile("fit_results/runScan-T90-pos0_large_scale_out_FR.root");

  TGraphErrors* FRAC_SIG_T77 = (TGraphErrors*)f_77->Get("FRAC_SIG");FRAC_SIG_T77->SetName("FRAC_SIG_T77");FRAC_SIG_T77->SetTitle("FRAC_SIG_T77");
  TGraphErrors* FRAC_SIG_T82 = (TGraphErrors*)f_82->Get("FRAC_SIG");FRAC_SIG_T82->SetName("FRAC_SIG_T82");FRAC_SIG_T82->SetTitle("FRAC_SIG_T82");
  TGraphErrors* FRAC_SIG_T84 = (TGraphErrors*)f_84->Get("FRAC_SIG");FRAC_SIG_T84->SetName("FRAC_SIG_T84");FRAC_SIG_T84->SetTitle("FRAC_SIG_T84");
  TGraphErrors* FRAC_SIG_T86 = (TGraphErrors*)f_86->Get("FRAC_SIG");FRAC_SIG_T86->SetName("FRAC_SIG_T86");FRAC_SIG_T86->SetTitle("FRAC_SIG_T86");
  TGraphErrors* FRAC_SIG_T88 = (TGraphErrors*)f_88->Get("FRAC_SIG");FRAC_SIG_T88->SetName("FRAC_SIG_T88");FRAC_SIG_T88->SetTitle("FRAC_SIG_T88");
  TGraphErrors* FRAC_SIG_T90 = (TGraphErrors*)f_90->Get("FRAC_SIG");FRAC_SIG_T90->SetName("FRAC_SIG_T90");FRAC_SIG_T90->SetTitle("FRAC_SIG_T90");

  FRAC_SIG_T77->SetMarkerStyle(20);FRAC_SIG_T77->SetMarkerColor(1);         FRAC_SIG_T77->SetMarkerSize(2);
  FRAC_SIG_T82->SetMarkerStyle(20);FRAC_SIG_T82->SetMarkerColor(2);         FRAC_SIG_T82->SetMarkerSize(2);
  FRAC_SIG_T84->SetMarkerStyle(20);FRAC_SIG_T84->SetMarkerColor(4);         FRAC_SIG_T84->SetMarkerSize(2);
  FRAC_SIG_T86->SetMarkerStyle(20);FRAC_SIG_T86->SetMarkerColor(8);         FRAC_SIG_T86->SetMarkerSize(2);
  FRAC_SIG_T88->SetMarkerStyle(20);FRAC_SIG_T88->SetMarkerColor(kOrange-2); FRAC_SIG_T88->SetMarkerSize(2);
  FRAC_SIG_T90->SetMarkerStyle(20);FRAC_SIG_T90->SetMarkerColor(11);        FRAC_SIG_T90->SetMarkerSize(2);

  
  TMultiGraph *mg_FRAC_SIG = new TMultiGraph();mg_FRAC_SIG->SetName("mg_FRAC_SIG");
  mg_FRAC_SIG->SetTitle("F_{sig} vs Channel");
  mg_FRAC_SIG->Add(FRAC_SIG_T77);
  mg_FRAC_SIG->Add(FRAC_SIG_T82);
  mg_FRAC_SIG->Add(FRAC_SIG_T84);
  mg_FRAC_SIG->Add(FRAC_SIG_T86);
  mg_FRAC_SIG->Add(FRAC_SIG_T88);
  mg_FRAC_SIG->Add(FRAC_SIG_T90);
  TCanvas *c_FRAC_SIG = new TCanvas("c_FRAC_SIG","c_FRAC_SIG");
  mg_FRAC_SIG->Draw("AP");
  mg_FRAC_SIG->GetXaxis()->SetTitle("Channel");
  mg_FRAC_SIG->GetYaxis()->SetTitle("F_{Sig} ");
  gPad->BuildLegend();
  gPad->Update();
  gPad->SetGridy();

  ///////////////////////////////////////////
  /////// T0 ////////////////////////////////
  ///////////////////////////////////////////
  
  TGraphErrors* T0_T77 = (TGraphErrors*)f_77->Get("MEAN_L_B_0");T0_T77->SetName("T0_T77");T0_T77->SetTitle("T0_T77");
  TGraphErrors* T0_T82 = (TGraphErrors*)f_82->Get("MEAN_L_B_0");T0_T82->SetName("T0_T82");T0_T82->SetTitle("T0_T82");
  TGraphErrors* T0_T84 = (TGraphErrors*)f_84->Get("MEAN_L_B_0");T0_T84->SetName("T0_T84");T0_T84->SetTitle("T0_T84");
  TGraphErrors* T0_T86 = (TGraphErrors*)f_86->Get("MEAN_L_B_0");T0_T86->SetName("T0_T86");T0_T86->SetTitle("T0_T86");
  TGraphErrors* T0_T88 = (TGraphErrors*)f_88->Get("MEAN_L_B_0");T0_T88->SetName("T0_T88");T0_T88->SetTitle("T0_T88");
  TGraphErrors* T0_T90 = (TGraphErrors*)f_90->Get("MEAN_L_B_0");T0_T90->SetName("T0_T90");T0_T90->SetTitle("T0_T90");

  T0_T77->SetMarkerStyle(20);T0_T77->SetMarkerColor(1);         T0_T77->SetMarkerSize(2);
  T0_T82->SetMarkerStyle(20);T0_T82->SetMarkerColor(2);         T0_T82->SetMarkerSize(2);
  T0_T84->SetMarkerStyle(20);T0_T84->SetMarkerColor(4);         T0_T84->SetMarkerSize(2);
  T0_T86->SetMarkerStyle(20);T0_T86->SetMarkerColor(8);         T0_T86->SetMarkerSize(2);
  T0_T88->SetMarkerStyle(20);T0_T88->SetMarkerColor(kOrange-2); T0_T88->SetMarkerSize(2);
  T0_T90->SetMarkerStyle(20);T0_T90->SetMarkerColor(11);        T0_T90->SetMarkerSize(2);

  
  TMultiGraph *mg_T0 = new TMultiGraph();mg_T0->SetName("mg_T0");
  mg_T0->SetTitle("T_{0} vs Channel");
  mg_T0->Add(T0_T77);
  mg_T0->Add(T0_T82);
  mg_T0->Add(T0_T84);
  mg_T0->Add(T0_T86);
  mg_T0->Add(T0_T88);
  mg_T0->Add(T0_T90);
  TCanvas *c_T0 = new TCanvas("c_T0","c_T0");
  mg_T0->Draw("AP");
  mg_T0->GetXaxis()->SetTitle("Channel");
  mg_T0->GetYaxis()->SetTitle("T_{0} [ns]");
  gPad->BuildLegend();
  gPad->Update();
  gPad->SetGridy();
  




  ////////////////////////////////////////////////
  /// Single pixel graphs //////
  //////////////////////////////


  Int_t n=14;
  
  float TUNES[6];
  TUNES[0]=77;
  TUNES[1]=82;
  TUNES[2]=84;
  TUNES[3]=86;
  TUNES[4]=88;
  TUNES[5]=90;
  float err_TUNES[6];
  err_TUNES[0]=0;
  err_TUNES[1]=0;
  err_TUNES[2]=0;
  err_TUNES[3]=0;
  err_TUNES[4]=0;
  err_TUNES[5]=0;


  Double_t x_au;
  Double_t FRAC_CH_au;
  float FRAC_CH[6];
  float err_FRAC_CH[6];
  
  FRAC_SIG_T77->GetPoint(n,x_au,FRAC_CH_au); FRAC_CH[0]=FRAC_CH_au;
  FRAC_SIG_T82->GetPoint(n,x_au,FRAC_CH_au); FRAC_CH[1]=FRAC_CH_au;
  FRAC_SIG_T84->GetPoint(n,x_au,FRAC_CH_au); FRAC_CH[2]=FRAC_CH_au;
  FRAC_SIG_T86->GetPoint(n,x_au,FRAC_CH_au); FRAC_CH[3]=FRAC_CH_au;
  FRAC_SIG_T88->GetPoint(n,x_au,FRAC_CH_au); FRAC_CH[4]=FRAC_CH_au;
  FRAC_SIG_T90->GetPoint(n,x_au,FRAC_CH_au); FRAC_CH[5]=FRAC_CH_au;

  err_FRAC_CH[0]=(float)FRAC_SIG_T77->GetErrorXhigh(n);
  err_FRAC_CH[1]=(float)FRAC_SIG_T82->GetErrorXhigh(n);
  err_FRAC_CH[2]=(float)FRAC_SIG_T84->GetErrorXhigh(n);
  err_FRAC_CH[3]=(float)FRAC_SIG_T86->GetErrorXhigh(n);
  err_FRAC_CH[4]=(float)FRAC_SIG_T88->GetErrorXhigh(n);
  err_FRAC_CH[5]=(float)FRAC_SIG_T90->GetErrorXhigh(n);

  
  TGraphErrors *TGE_FRAC_CH = new TGraphErrors(6,TUNES,FRAC_CH,err_TUNES,err_FRAC_CH);TGE_FRAC_CH->SetName("FRAC_CH 4"); TGE_FRAC_CH->SetTitle("FRAC_CH 4"); 
  TGE_FRAC_CH->SetMarkerStyle(20);TGE_FRAC_CH->SetMarkerColor(4);         TGE_FRAC_CH->SetMarkerSize(2);
  TCanvas* c_FRAC_CH = new TCanvas("c_FRAC_CH","c_FRAC_CH");
  TGE_FRAC_CH->Draw("AP");
  TGE_FRAC_CH->GetXaxis()->SetTitle("Laser Tune");
  TGE_FRAC_CH->GetYaxis()->SetTitle("F_{Sig}");
  gPad->BuildLegend();
  gPad->Update();
  gPad->SetGridy();











  
  Double_t T0_au;
  float T0[6];
  float err_T0[6];
  
  T0_T77->GetPoint(n,x_au,T0_au); T0[0]=T0_au;
  T0_T82->GetPoint(n,x_au,T0_au); T0[1]=T0_au;
  T0_T84->GetPoint(n,x_au,T0_au); T0[2]=T0_au;
  T0_T86->GetPoint(n,x_au,T0_au); T0[3]=T0_au;
  T0_T88->GetPoint(n,x_au,T0_au); T0[4]=T0_au;
  T0_T90->GetPoint(n,x_au,T0_au); T0[5]=T0_au;

  err_T0[0]=T0_T77->GetErrorX(n);
  err_T0[1]=T0_T82->GetErrorX(n);
  err_T0[2]=T0_T84->GetErrorX(n);
  err_T0[3]=T0_T86->GetErrorX(n);
  err_T0[4]=T0_T88->GetErrorX(n);
  err_T0[5]=T0_T90->GetErrorX(n);

  
  TGraphErrors *TGE_T0 = new TGraphErrors(6,TUNES,T0,err_TUNES,err_T0);TGE_T0->SetName("T0 che 4"); TGE_T0->SetTitle("T0 ch 4"); 
  TGE_T0->SetMarkerStyle(20);TGE_T0->SetMarkerColor(4);         TGE_T0->SetMarkerSize(2);
  TCanvas* c_T0_CH = new TCanvas("c_T0_CH","c_T0_CH");
  TGE_T0->Draw("AP");
  TGE_T0->GetXaxis()->SetTitle("Laser Tune");
  TGE_T0->GetYaxis()->SetTitle("T_{0} [ns]");
  gPad->BuildLegend();
  gPad->Update();
  gPad->SetGridy();

  


  //////////////////////







  
}
