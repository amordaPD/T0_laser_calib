#include <iostream>
using namespace std;

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include "TAxis.h"
#include <TGraph.h>

void Analyse_ScanTune(){
  ///////////////////////////////////////////////
  /////// channel 4 behaves quite well  /////////
  ///////////////////////////////////////////////

  TString file_path="/";
  //TFile *f_77 = new TFile("fit_results/runScan-T77-pos0_large_scale_out_FR.root");
  TFile *f_80 = new TFile("fit_results"+file_path+"T80_FR.root");
  TFile *f_82 = new TFile("fit_results"+file_path+"T82_FR.root");
  TFile *f_84 = new TFile("fit_results"+file_path+"T84_FR.root");
  TFile *f_86 = new TFile("fit_results"+file_path+"T86_FR.root");
  TFile *f_88 = new TFile("fit_results"+file_path+"T88_FR.root");
  TFile *f_90 = new TFile("fit_results"+file_path+"T90_FR.root");
  TFile *f_78 = new TFile("fit_results"+file_path+"T78_FR.root");
  TFile *f_76 = new TFile("fit_results"+file_path+"T76_FR.root");
  TFile *f_74 = new TFile("fit_results"+file_path+"T74_FR.root");
  TFile *f_72 = new TFile("fit_results"+file_path+"T72_FR.root");
  TFile *f_73 = new TFile("fit_results"+file_path+"T73_FR.root");
  TFile *f_75 = new TFile("fit_results"+file_path+"T75_FR.root");
  TFile *f_77 = new TFile("fit_results"+file_path+"T77_FR.root");
  TFile *f_79 = new TFile("fit_results"+file_path+"T79_FR.root");
  TFile *f_81 = new TFile("fit_results"+file_path+"T81_FR.root");
  TFile *f_795 = new TFile("fit_results"+file_path+"T79_FR.root");

  //TGraphErrors* FRAC_SIG_T77 = (TGraphErrors*)f_77->Get("FRAC_SIG");FRAC_SIG_T77->SetName("FRAC_SIG_T77");FRAC_SIG_T77->SetTitle("FRAC_SIG_T77");
  TGraphErrors* FRAC_SIG_T80 = (TGraphErrors*)f_80->Get("FRAC_SIG");FRAC_SIG_T80->SetName("FRAC_SIG_T80");FRAC_SIG_T80->SetTitle("FRAC_SIG_T80");
  TGraphErrors* FRAC_SIG_T82 = (TGraphErrors*)f_82->Get("FRAC_SIG");FRAC_SIG_T82->SetName("FRAC_SIG_T82");FRAC_SIG_T82->SetTitle("FRAC_SIG_T82");
  TGraphErrors* FRAC_SIG_T84 = (TGraphErrors*)f_84->Get("FRAC_SIG");FRAC_SIG_T84->SetName("FRAC_SIG_T84");FRAC_SIG_T84->SetTitle("FRAC_SIG_T84");
  TGraphErrors* FRAC_SIG_T86 = (TGraphErrors*)f_86->Get("FRAC_SIG");FRAC_SIG_T86->SetName("FRAC_SIG_T86");FRAC_SIG_T86->SetTitle("FRAC_SIG_T86");
  TGraphErrors* FRAC_SIG_T88 = (TGraphErrors*)f_88->Get("FRAC_SIG");FRAC_SIG_T88->SetName("FRAC_SIG_T88");FRAC_SIG_T88->SetTitle("FRAC_SIG_T88");
  TGraphErrors* FRAC_SIG_T90 = (TGraphErrors*)f_90->Get("FRAC_SIG");FRAC_SIG_T90->SetName("FRAC_SIG_T90");FRAC_SIG_T90->SetTitle("FRAC_SIG_T90");
  TGraphErrors* FRAC_SIG_T78 = (TGraphErrors*)f_78->Get("FRAC_SIG");FRAC_SIG_T78->SetName("FRAC_SIG_T78");FRAC_SIG_T78->SetTitle("FRAC_SIG_T78");
  TGraphErrors* FRAC_SIG_T76 = (TGraphErrors*)f_76->Get("FRAC_SIG");FRAC_SIG_T76->SetName("FRAC_SIG_T76");FRAC_SIG_T76->SetTitle("FRAC_SIG_T76");
  TGraphErrors* FRAC_SIG_T74 = (TGraphErrors*)f_74->Get("FRAC_SIG");FRAC_SIG_T74->SetName("FRAC_SIG_T74");FRAC_SIG_T74->SetTitle("FRAC_SIG_T74");
  TGraphErrors* FRAC_SIG_T72 = (TGraphErrors*)f_72->Get("FRAC_SIG");FRAC_SIG_T72->SetName("FRAC_SIG_T72");FRAC_SIG_T72->SetTitle("FRAC_SIG_T72");
  TGraphErrors* FRAC_SIG_T73 = (TGraphErrors*)f_73->Get("FRAC_SIG");FRAC_SIG_T73->SetName("FRAC_SIG_T73");FRAC_SIG_T73->SetTitle("FRAC_SIG_T73");
  TGraphErrors* FRAC_SIG_T75 = (TGraphErrors*)f_75->Get("FRAC_SIG");FRAC_SIG_T75->SetName("FRAC_SIG_T75");FRAC_SIG_T75->SetTitle("FRAC_SIG_T75");
  TGraphErrors* FRAC_SIG_T77 = (TGraphErrors*)f_77->Get("FRAC_SIG");FRAC_SIG_T77->SetName("FRAC_SIG_T77");FRAC_SIG_T77->SetTitle("FRAC_SIG_T77");
  TGraphErrors* FRAC_SIG_T79 = (TGraphErrors*)f_79->Get("FRAC_SIG");FRAC_SIG_T79->SetName("FRAC_SIG_T79");FRAC_SIG_T79->SetTitle("FRAC_SIG_T79");
  TGraphErrors* FRAC_SIG_T81 = (TGraphErrors*)f_81->Get("FRAC_SIG");FRAC_SIG_T81->SetName("FRAC_SIG_T81");FRAC_SIG_T81->SetTitle("FRAC_SIG_T81");
  TGraphErrors* FRAC_SIG_T795 = (TGraphErrors*)f_795->Get("FRAC_SIG");FRAC_SIG_T795->SetName("FRAC_SIG_T795");FRAC_SIG_T795->SetTitle("FRAC_SIG_T795");

  
  FRAC_SIG_T80->SetMarkerStyle(20);FRAC_SIG_T80->SetMarkerColor(9);         FRAC_SIG_T80->SetMarkerSize(2);
  FRAC_SIG_T82->SetMarkerStyle(20);FRAC_SIG_T82->SetMarkerColor(2);         FRAC_SIG_T82->SetMarkerSize(2);
  FRAC_SIG_T84->SetMarkerStyle(20);FRAC_SIG_T84->SetMarkerColor(4);         FRAC_SIG_T84->SetMarkerSize(2);
  FRAC_SIG_T86->SetMarkerStyle(20);FRAC_SIG_T86->SetMarkerColor(8);         FRAC_SIG_T86->SetMarkerSize(2);
  FRAC_SIG_T88->SetMarkerStyle(20);FRAC_SIG_T88->SetMarkerColor(kOrange-2); FRAC_SIG_T88->SetMarkerSize(2);
  FRAC_SIG_T90->SetMarkerStyle(20);FRAC_SIG_T90->SetMarkerColor(11);        FRAC_SIG_T90->SetMarkerSize(2);
  FRAC_SIG_T78->SetMarkerStyle(20);FRAC_SIG_T78->SetMarkerColor(7);         FRAC_SIG_T78->SetMarkerSize(2);
  FRAC_SIG_T76->SetMarkerStyle(20);FRAC_SIG_T76->SetMarkerColor(6);         FRAC_SIG_T76->SetMarkerSize(2);
  FRAC_SIG_T74->SetMarkerStyle(20);FRAC_SIG_T74->SetMarkerColor(5);         FRAC_SIG_T74->SetMarkerSize(2);
  FRAC_SIG_T73->SetMarkerStyle(20);FRAC_SIG_T73->SetMarkerColor(12);        FRAC_SIG_T73->SetMarkerSize(2);
  FRAC_SIG_T72->SetMarkerStyle(20);FRAC_SIG_T72->SetMarkerColor(3);         FRAC_SIG_T72->SetMarkerSize(2);
  FRAC_SIG_T75->SetMarkerStyle(20);FRAC_SIG_T75->SetMarkerColor(13);        FRAC_SIG_T75->SetMarkerSize(2);
  FRAC_SIG_T77->SetMarkerStyle(20);FRAC_SIG_T77->SetMarkerColor(1);         FRAC_SIG_T77->SetMarkerSize(2);
  FRAC_SIG_T79->SetMarkerStyle(20);FRAC_SIG_T79->SetMarkerColor(15);        FRAC_SIG_T79->SetMarkerSize(2);
  FRAC_SIG_T81->SetMarkerStyle(20);FRAC_SIG_T81->SetMarkerColor(16);         FRAC_SIG_T81->SetMarkerSize(2);
  FRAC_SIG_T795->SetMarkerStyle(20);FRAC_SIG_T795->SetMarkerColor(17);      FRAC_SIG_T795->SetMarkerSize(2);

  
  TMultiGraph *mg_FRAC_SIG = new TMultiGraph();mg_FRAC_SIG->SetName("mg_FRAC_SIG");
  mg_FRAC_SIG->SetTitle("F_{sig} vs Channel");
  //mg_FRAC_SIG->Add(FRAC_SIG_T77);
  mg_FRAC_SIG->Add(FRAC_SIG_T80);
  mg_FRAC_SIG->Add(FRAC_SIG_T82);
  mg_FRAC_SIG->Add(FRAC_SIG_T84);
  mg_FRAC_SIG->Add(FRAC_SIG_T86);
  mg_FRAC_SIG->Add(FRAC_SIG_T88);
  mg_FRAC_SIG->Add(FRAC_SIG_T90);
  mg_FRAC_SIG->Add(FRAC_SIG_T78);
  mg_FRAC_SIG->Add(FRAC_SIG_T76);
  mg_FRAC_SIG->Add(FRAC_SIG_T74);
  mg_FRAC_SIG->Add(FRAC_SIG_T72);
  mg_FRAC_SIG->Add(FRAC_SIG_T73);
  mg_FRAC_SIG->Add(FRAC_SIG_T75);
  mg_FRAC_SIG->Add(FRAC_SIG_T77);
  mg_FRAC_SIG->Add(FRAC_SIG_T79);
  mg_FRAC_SIG->Add(FRAC_SIG_T81);
  mg_FRAC_SIG->Add(FRAC_SIG_T795);
  
  TCanvas *c_FRAC_SIG = new TCanvas("c_FRAC_SIG","c_FRAC_SIG");
  mg_FRAC_SIG->Draw("AP");
  mg_FRAC_SIG->GetYaxis()->SetRangeUser(0.65,1);
  mg_FRAC_SIG->GetXaxis()->SetTitle("Channel");
  mg_FRAC_SIG->GetYaxis()->SetTitle("F_{Sig} ");
  gPad->BuildLegend();
  gPad->Update();
  gPad->SetGridy();

  ///////////////////////////////////////////
  /////// T0 ////////////////////////////////
  ///////////////////////////////////////////
  
  //TGraphErrors* T0_T77 = (TGraphErrors*)f_77->Get("MEAN_L_B_0");T0_T77->SetName("T0_T77");T0_T77->SetTitle("T0_T77");
  TGraphErrors* T0_T80 = (TGraphErrors*)f_80->Get("MEAN_L_B_0");T0_T80->SetName("T0_T80");T0_T80->SetTitle("T0_T80");
  TGraphErrors* T0_T82 = (TGraphErrors*)f_82->Get("MEAN_L_B_0");T0_T82->SetName("T0_T82");T0_T82->SetTitle("T0_T82");
  TGraphErrors* T0_T84 = (TGraphErrors*)f_84->Get("MEAN_L_B_0");T0_T84->SetName("T0_T84");T0_T84->SetTitle("T0_T84");
  TGraphErrors* T0_T86 = (TGraphErrors*)f_86->Get("MEAN_L_B_0");T0_T86->SetName("T0_T86");T0_T86->SetTitle("T0_T86");
  TGraphErrors* T0_T88 = (TGraphErrors*)f_88->Get("MEAN_L_B_0");T0_T88->SetName("T0_T88");T0_T88->SetTitle("T0_T88");
  TGraphErrors* T0_T90 = (TGraphErrors*)f_90->Get("MEAN_L_B_0");T0_T90->SetName("T0_T90");T0_T90->SetTitle("T0_T90");
  TGraphErrors* T0_T78 = (TGraphErrors*)f_78->Get("MEAN_L_B_0");T0_T78->SetName("T0_T78");T0_T88->SetTitle("T0_T78");
  TGraphErrors* T0_T76 = (TGraphErrors*)f_76->Get("MEAN_L_B_0");T0_T76->SetName("T0_T76");T0_T76->SetTitle("T0_T76");
  TGraphErrors* T0_T74 = (TGraphErrors*)f_74->Get("MEAN_L_B_0");T0_T74->SetName("T0_T74");T0_T74->SetTitle("T0_T74");
  TGraphErrors* T0_T72 = (TGraphErrors*)f_72->Get("MEAN_L_B_0");T0_T72->SetName("T0_T72");T0_T72->SetTitle("T0_T72");
  TGraphErrors* T0_T73 = (TGraphErrors*)f_73->Get("MEAN_L_B_0");T0_T73->SetName("T0_T73");T0_T73->SetTitle("T0_T73");
  TGraphErrors* T0_T75 = (TGraphErrors*)f_75->Get("MEAN_L_B_0");T0_T75->SetName("T0_T75");T0_T75->SetTitle("T0_T75");
  TGraphErrors* T0_T77 = (TGraphErrors*)f_77->Get("MEAN_L_B_0");T0_T77->SetName("T0_T77");T0_T77->SetTitle("T0_T77");
  TGraphErrors* T0_T79 = (TGraphErrors*)f_79->Get("MEAN_L_B_0");T0_T79->SetName("T0_T79");T0_T79->SetTitle("T0_T79");
  TGraphErrors* T0_T81 = (TGraphErrors*)f_81->Get("MEAN_L_B_0");T0_T81->SetName("T0_T81");T0_T81->SetTitle("T0_T81");
  TGraphErrors* T0_T795 = (TGraphErrors*)f_795->Get("MEAN_L_B_0");T0_T795->SetName("T0_T795");T0_T795->SetTitle("T0_T795");


  T0_T80->SetMarkerStyle(20);T0_T80->SetMarkerColor(9);         T0_T80->SetMarkerSize(2);
  T0_T82->SetMarkerStyle(20);T0_T82->SetMarkerColor(2);         T0_T82->SetMarkerSize(2);
  T0_T84->SetMarkerStyle(20);T0_T84->SetMarkerColor(4);         T0_T84->SetMarkerSize(2);
  T0_T86->SetMarkerStyle(20);T0_T86->SetMarkerColor(8);         T0_T86->SetMarkerSize(2);
  T0_T88->SetMarkerStyle(20);T0_T88->SetMarkerColor(kOrange-2); T0_T88->SetMarkerSize(2);
  T0_T90->SetMarkerStyle(20);T0_T90->SetMarkerColor(11);        T0_T90->SetMarkerSize(2);
  T0_T78->SetMarkerStyle(20);T0_T78->SetMarkerColor(7);         T0_T78->SetMarkerSize(2);
  T0_T76->SetMarkerStyle(20);T0_T76->SetMarkerColor(6);         T0_T76->SetMarkerSize(2);
  T0_T74->SetMarkerStyle(20);T0_T74->SetMarkerColor(5);         T0_T74->SetMarkerSize(2);
  T0_T72->SetMarkerStyle(20);T0_T72->SetMarkerColor(12);        T0_T72->SetMarkerSize(2);
  T0_T73->SetMarkerStyle(20);T0_T73->SetMarkerColor(3);         T0_T73->SetMarkerSize(2);
  T0_T75->SetMarkerStyle(20);T0_T75->SetMarkerColor(13);        T0_T75->SetMarkerSize(2);
  T0_T77->SetMarkerStyle(20);T0_T77->SetMarkerColor(3);         T0_T77->SetMarkerSize(2);
  T0_T79->SetMarkerStyle(20);T0_T79->SetMarkerColor(13);        T0_T79->SetMarkerSize(2);
  T0_T81->SetMarkerStyle(20);T0_T81->SetMarkerColor(3);         T0_T81->SetMarkerSize(2);
  T0_T795->SetMarkerStyle(20);T0_T795->SetMarkerColor(13);      T0_T795->SetMarkerSize(2);

  
  TMultiGraph *mg_T0 = new TMultiGraph();mg_T0->SetName("mg_T0");
  mg_T0->SetTitle("T_{0} vs Channel");
  //mg_T0->Add(T0_T77);
  mg_T0->Add(T0_T80);
  mg_T0->Add(T0_T82);
  mg_T0->Add(T0_T84);
  mg_T0->Add(T0_T86);
  mg_T0->Add(T0_T88);
  mg_T0->Add(T0_T90);
  mg_T0->Add(T0_T78);
  mg_T0->Add(T0_T76);
  mg_T0->Add(T0_T74);
  mg_T0->Add(T0_T72);
  mg_T0->Add(T0_T73);
  mg_T0->Add(T0_T75);
  mg_T0->Add(T0_T77);
  mg_T0->Add(T0_T79);
  mg_T0->Add(T0_T81);
  mg_T0->Add(T0_T795);
  TCanvas *c_T0 = new TCanvas("c_T0","c_T0");
  mg_T0->Draw("AP");
  mg_T0->GetXaxis()->SetTitle("Channel");
  mg_T0->GetYaxis()->SetTitle("T_{0} [ns]");
  gPad->BuildLegend();
  gPad->Update();
  gPad->SetGridy();
  


  /////////////////////////////////////
  /// Amplitude ///////////////////////
  /////////////////////////////////////
  TGraphErrors* AMP_SIG_T80 = (TGraphErrors*)f_80->Get("AMP_SIG");AMP_SIG_T80->SetName("AMP_SIG_T80");AMP_SIG_T80->SetTitle("AMP_SIG_T80");
  TGraphErrors* AMP_SIG_T82 = (TGraphErrors*)f_82->Get("AMP_SIG");AMP_SIG_T82->SetName("AMP_SIG_T82");AMP_SIG_T82->SetTitle("AMP_SIG_T82");
  TGraphErrors* AMP_SIG_T84 = (TGraphErrors*)f_84->Get("AMP_SIG");AMP_SIG_T84->SetName("AMP_SIG_T84");AMP_SIG_T84->SetTitle("AMP_SIG_T84");
  TGraphErrors* AMP_SIG_T86 = (TGraphErrors*)f_86->Get("AMP_SIG");AMP_SIG_T86->SetName("AMP_SIG_T86");AMP_SIG_T86->SetTitle("AMP_SIG_T86");
  TGraphErrors* AMP_SIG_T88 = (TGraphErrors*)f_88->Get("AMP_SIG");AMP_SIG_T88->SetName("AMP_SIG_T88");AMP_SIG_T88->SetTitle("AMP_SIG_T88");
  TGraphErrors* AMP_SIG_T90 = (TGraphErrors*)f_90->Get("AMP_SIG");AMP_SIG_T90->SetName("AMP_SIG_T90");AMP_SIG_T90->SetTitle("AMP_SIG_T90");
  TGraphErrors* AMP_SIG_T78 = (TGraphErrors*)f_78->Get("AMP_SIG");AMP_SIG_T78->SetName("AMP_SIG_T78");AMP_SIG_T78->SetTitle("AMP_SIG_T78");
  TGraphErrors* AMP_SIG_T76 = (TGraphErrors*)f_76->Get("AMP_SIG");AMP_SIG_T76->SetName("AMP_SIG_T76");AMP_SIG_T76->SetTitle("AMP_SIG_T76");
  TGraphErrors* AMP_SIG_T74 = (TGraphErrors*)f_74->Get("AMP_SIG");AMP_SIG_T74->SetName("AMP_SIG_T74");AMP_SIG_T74->SetTitle("AMP_SIG_T74");
  TGraphErrors* AMP_SIG_T72 = (TGraphErrors*)f_72->Get("AMP_SIG");AMP_SIG_T72->SetName("AMP_SIG_T72");AMP_SIG_T72->SetTitle("AMP_SIG_T72");
  TGraphErrors* AMP_SIG_T73 = (TGraphErrors*)f_73->Get("AMP_SIG");AMP_SIG_T73->SetName("AMP_SIG_T73");AMP_SIG_T73->SetTitle("AMP_SIG_T73");
  TGraphErrors* AMP_SIG_T75 = (TGraphErrors*)f_75->Get("AMP_SIG");AMP_SIG_T75->SetName("AMP_SIG_T75");AMP_SIG_T75->SetTitle("AMP_SIG_T75");
  TGraphErrors* AMP_SIG_T77 = (TGraphErrors*)f_77->Get("AMP_SIG");AMP_SIG_T77->SetName("AMP_SIG_T77");AMP_SIG_T77->SetTitle("AMP_SIG_T77");
  TGraphErrors* AMP_SIG_T79 = (TGraphErrors*)f_79->Get("AMP_SIG");AMP_SIG_T79->SetName("AMP_SIG_T79");AMP_SIG_T79->SetTitle("AMP_SIG_T79");
  TGraphErrors* AMP_SIG_T81 = (TGraphErrors*)f_81->Get("AMP_SIG");AMP_SIG_T81->SetName("AMP_SIG_T81");AMP_SIG_T81->SetTitle("AMP_SIG_T81");
  TGraphErrors* AMP_SIG_T795 = (TGraphErrors*)f_795->Get("AMP_SIG");AMP_SIG_T795->SetName("AMP_SIG_T795");AMP_SIG_T795->SetTitle("AMP_SIG_T795");

  // AMP_SIG_T77->SetMarkerStyle(20);AMP_SIG_T77->SetMarkerColor(1);         AMP_SIG_T77->SetMarkerSize(2);
  AMP_SIG_T80->SetMarkerStyle(20);AMP_SIG_T80->SetMarkerColor(9);         AMP_SIG_T80->SetMarkerSize(2);
  AMP_SIG_T82->SetMarkerStyle(20);AMP_SIG_T82->SetMarkerColor(2);         AMP_SIG_T82->SetMarkerSize(2);
  AMP_SIG_T84->SetMarkerStyle(20);AMP_SIG_T84->SetMarkerColor(4);         AMP_SIG_T84->SetMarkerSize(2);
  AMP_SIG_T86->SetMarkerStyle(20);AMP_SIG_T86->SetMarkerColor(8);         AMP_SIG_T86->SetMarkerSize(2);
  AMP_SIG_T88->SetMarkerStyle(20);AMP_SIG_T88->SetMarkerColor(kOrange-2); AMP_SIG_T88->SetMarkerSize(2);
  AMP_SIG_T90->SetMarkerStyle(20);AMP_SIG_T90->SetMarkerColor(11);        AMP_SIG_T90->SetMarkerSize(2);
  AMP_SIG_T78->SetMarkerStyle(20);AMP_SIG_T78->SetMarkerColor(7);         AMP_SIG_T78->SetMarkerSize(2);
  AMP_SIG_T76->SetMarkerStyle(20);AMP_SIG_T76->SetMarkerColor(6);         AMP_SIG_T76->SetMarkerSize(2);
  AMP_SIG_T74->SetMarkerStyle(20);AMP_SIG_T74->SetMarkerColor(5);         AMP_SIG_T74->SetMarkerSize(2);
  AMP_SIG_T72->SetMarkerStyle(20);AMP_SIG_T72->SetMarkerColor(12);        AMP_SIG_T72->SetMarkerSize(2);
  AMP_SIG_T73->SetMarkerStyle(20);AMP_SIG_T73->SetMarkerColor(3);         AMP_SIG_T73->SetMarkerSize(2);
  AMP_SIG_T75->SetMarkerStyle(20);AMP_SIG_T75->SetMarkerColor(13);        AMP_SIG_T75->SetMarkerSize(2);
  AMP_SIG_T77->SetMarkerStyle(20);AMP_SIG_T77->SetMarkerColor(3);         AMP_SIG_T77->SetMarkerSize(2);
  AMP_SIG_T79->SetMarkerStyle(20);AMP_SIG_T79->SetMarkerColor(13);        AMP_SIG_T79->SetMarkerSize(2);
  AMP_SIG_T81->SetMarkerStyle(20);AMP_SIG_T81->SetMarkerColor(3);         AMP_SIG_T81->SetMarkerSize(2);
  AMP_SIG_T795->SetMarkerStyle(20);AMP_SIG_T795->SetMarkerColor(13);        AMP_SIG_T795->SetMarkerSize(2);

  
  TMultiGraph *mg_AMP_SIG = new TMultiGraph();mg_AMP_SIG->SetName("mg_AMP_SIG");
  mg_AMP_SIG->SetTitle("Amplitude peak position vs Channel");
  //mg_AMP_SIG->Add(AMP_SIG_T77);
  mg_AMP_SIG->Add(AMP_SIG_T80);
  mg_AMP_SIG->Add(AMP_SIG_T82);
  mg_AMP_SIG->Add(AMP_SIG_T84);
  mg_AMP_SIG->Add(AMP_SIG_T86);
  mg_AMP_SIG->Add(AMP_SIG_T88);
  mg_AMP_SIG->Add(AMP_SIG_T90);
  mg_AMP_SIG->Add(AMP_SIG_T78);
  mg_AMP_SIG->Add(AMP_SIG_T76);
  mg_AMP_SIG->Add(AMP_SIG_T74);
  mg_AMP_SIG->Add(AMP_SIG_T72);
  mg_AMP_SIG->Add(AMP_SIG_T73);
  mg_AMP_SIG->Add(AMP_SIG_T75);
  mg_AMP_SIG->Add(AMP_SIG_T77);
  mg_AMP_SIG->Add(AMP_SIG_T79);
  mg_AMP_SIG->Add(AMP_SIG_T81);
  mg_AMP_SIG->Add(AMP_SIG_T795);
  
  TCanvas *c_AMP_SIG = new TCanvas("c_AMP_SIG","c_AMP_SIG");
  mg_AMP_SIG->Draw("AP");
  mg_AMP_SIG->GetYaxis()->SetRangeUser(115,140);
  mg_AMP_SIG->GetXaxis()->SetTitle("Channel");
  mg_AMP_SIG->GetYaxis()->SetTitle("Amplitude peak position [ADC counts] ");
  gPad->BuildLegend();
  gPad->Update();
  gPad->SetGridy();



  ////////////////////////
  ////Time separation ////
  ////////////////////////

  TGraphErrors* DELTA_LSP_T0_T80 = (TGraphErrors*)f_80->Get("DELTA_LSP_T0");DELTA_LSP_T0_T80->SetName("DELTA_LSP_T0_T80");DELTA_LSP_T0_T80->SetTitle("DELTA_LSP_T0_T80");
  TGraphErrors* DELTA_LSP_T0_T82 = (TGraphErrors*)f_82->Get("DELTA_LSP_T0");DELTA_LSP_T0_T82->SetName("DELTA_LSP_T0_T82");DELTA_LSP_T0_T82->SetTitle("DELTA_LSP_T0_T82");
  TGraphErrors* DELTA_LSP_T0_T84 = (TGraphErrors*)f_84->Get("DELTA_LSP_T0");DELTA_LSP_T0_T84->SetName("DELTA_LSP_T0_T84");DELTA_LSP_T0_T84->SetTitle("DELTA_LSP_T0_T84");
  TGraphErrors* DELTA_LSP_T0_T86 = (TGraphErrors*)f_86->Get("DELTA_LSP_T0");DELTA_LSP_T0_T86->SetName("DELTA_LSP_T0_T86");DELTA_LSP_T0_T86->SetTitle("DELTA_LSP_T0_T86");
  TGraphErrors* DELTA_LSP_T0_T88 = (TGraphErrors*)f_88->Get("DELTA_LSP_T0");DELTA_LSP_T0_T88->SetName("DELTA_LSP_T0_T88");DELTA_LSP_T0_T88->SetTitle("DELTA_LSP_T0_T88");
  TGraphErrors* DELTA_LSP_T0_T90 = (TGraphErrors*)f_90->Get("DELTA_LSP_T0");DELTA_LSP_T0_T90->SetName("DELTA_LSP_T0_T90");DELTA_LSP_T0_T90->SetTitle("DELTA_LSP_T0_T90");
  TGraphErrors* DELTA_LSP_T0_T78 = (TGraphErrors*)f_78->Get("DELTA_LSP_T0");DELTA_LSP_T0_T78->SetName("DELTA_LSP_T0_T78");DELTA_LSP_T0_T78->SetTitle("DELTA_LSP_T0_T78");
  TGraphErrors* DELTA_LSP_T0_T76 = (TGraphErrors*)f_76->Get("DELTA_LSP_T0");DELTA_LSP_T0_T76->SetName("DELTA_LSP_T0_T76");DELTA_LSP_T0_T76->SetTitle("DELTA_LSP_T0_T76");
  TGraphErrors* DELTA_LSP_T0_T74 = (TGraphErrors*)f_74->Get("DELTA_LSP_T0");DELTA_LSP_T0_T74->SetName("DELTA_LSP_T0_T74");DELTA_LSP_T0_T74->SetTitle("DELTA_LSP_T0_T74");
  TGraphErrors* DELTA_LSP_T0_T72 = (TGraphErrors*)f_72->Get("DELTA_LSP_T0");DELTA_LSP_T0_T72->SetName("DELTA_LSP_T0_T72");DELTA_LSP_T0_T72->SetTitle("DELTA_LSP_T0_T72");
  TGraphErrors* DELTA_LSP_T0_T73 = (TGraphErrors*)f_73->Get("DELTA_LSP_T0");DELTA_LSP_T0_T73->SetName("DELTA_LSP_T0_T73");DELTA_LSP_T0_T73->SetTitle("DELTA_LSP_T0_T73");
  TGraphErrors* DELTA_LSP_T0_T75 = (TGraphErrors*)f_75->Get("DELTA_LSP_T0");DELTA_LSP_T0_T75->SetName("DELTA_LSP_T0_T75");DELTA_LSP_T0_T75->SetTitle("DELTA_LSP_T0_T75");
  TGraphErrors* DELTA_LSP_T0_T77 = (TGraphErrors*)f_77->Get("DELTA_LSP_T0");DELTA_LSP_T0_T77->SetName("DELTA_LSP_T0_T77");DELTA_LSP_T0_T77->SetTitle("DELTA_LSP_T0_T77");
  TGraphErrors* DELTA_LSP_T0_T79 = (TGraphErrors*)f_79->Get("DELTA_LSP_T0");DELTA_LSP_T0_T79->SetName("DELTA_LSP_T0_T79");DELTA_LSP_T0_T79->SetTitle("DELTA_LSP_T0_T79");
  TGraphErrors* DELTA_LSP_T0_T81 = (TGraphErrors*)f_81->Get("DELTA_LSP_T0");DELTA_LSP_T0_T81->SetName("DELTA_LSP_T0_T81");DELTA_LSP_T0_T81->SetTitle("DELTA_LSP_T0_T81");
  TGraphErrors* DELTA_LSP_T0_T795 = (TGraphErrors*)f_795->Get("DELTA_LSP_T0");DELTA_LSP_T0_T795->SetName("DELTA_LSP_T0_T795");DELTA_LSP_T0_T795->SetTitle("DELTA_LSP_T0_T795");

  // DELTA_LSP_T0_T77->SetMarkerStyle(20);DELTA_LSP_T0_T77->SetMarkerColor(1);         DELTA_LSP_T0_T77->SetMarkerSize(2);
  DELTA_LSP_T0_T80->SetMarkerStyle(20);DELTA_LSP_T0_T80->SetMarkerColor(9);         DELTA_LSP_T0_T80->SetMarkerSize(2);
  DELTA_LSP_T0_T82->SetMarkerStyle(20);DELTA_LSP_T0_T82->SetMarkerColor(2);         DELTA_LSP_T0_T82->SetMarkerSize(2);
  DELTA_LSP_T0_T84->SetMarkerStyle(20);DELTA_LSP_T0_T84->SetMarkerColor(4);         DELTA_LSP_T0_T84->SetMarkerSize(2);
  DELTA_LSP_T0_T86->SetMarkerStyle(20);DELTA_LSP_T0_T86->SetMarkerColor(8);         DELTA_LSP_T0_T86->SetMarkerSize(2);
  DELTA_LSP_T0_T88->SetMarkerStyle(20);DELTA_LSP_T0_T88->SetMarkerColor(kOrange-2); DELTA_LSP_T0_T88->SetMarkerSize(2);
  DELTA_LSP_T0_T90->SetMarkerStyle(20);DELTA_LSP_T0_T90->SetMarkerColor(11);        DELTA_LSP_T0_T90->SetMarkerSize(2);
  DELTA_LSP_T0_T78->SetMarkerStyle(20);DELTA_LSP_T0_T78->SetMarkerColor(7);         DELTA_LSP_T0_T78->SetMarkerSize(2);
  DELTA_LSP_T0_T76->SetMarkerStyle(20);DELTA_LSP_T0_T76->SetMarkerColor(6);         DELTA_LSP_T0_T76->SetMarkerSize(2);
  DELTA_LSP_T0_T74->SetMarkerStyle(20);DELTA_LSP_T0_T74->SetMarkerColor(5);         DELTA_LSP_T0_T74->SetMarkerSize(2);
  DELTA_LSP_T0_T72->SetMarkerStyle(20);DELTA_LSP_T0_T72->SetMarkerColor(12);        DELTA_LSP_T0_T72->SetMarkerSize(2);
  DELTA_LSP_T0_T73->SetMarkerStyle(20);DELTA_LSP_T0_T73->SetMarkerColor(3);         DELTA_LSP_T0_T73->SetMarkerSize(2);
  DELTA_LSP_T0_T75->SetMarkerStyle(20);DELTA_LSP_T0_T75->SetMarkerColor(13);        DELTA_LSP_T0_T75->SetMarkerSize(2);
  DELTA_LSP_T0_T77->SetMarkerStyle(20);DELTA_LSP_T0_T77->SetMarkerColor(3);         DELTA_LSP_T0_T77->SetMarkerSize(2);
  DELTA_LSP_T0_T79->SetMarkerStyle(20);DELTA_LSP_T0_T79->SetMarkerColor(13);        DELTA_LSP_T0_T79->SetMarkerSize(2);
  DELTA_LSP_T0_T81->SetMarkerStyle(20);DELTA_LSP_T0_T81->SetMarkerColor(3);         DELTA_LSP_T0_T81->SetMarkerSize(2);
  DELTA_LSP_T0_T795->SetMarkerStyle(20);DELTA_LSP_T0_T795->SetMarkerColor(13);        DELTA_LSP_T0_T795->SetMarkerSize(2);



  
  TGraphErrors* DELTA_HSP_LSP_T80 = (TGraphErrors*)f_80->Get("DELTA_HSP_LSP");DELTA_HSP_LSP_T80->SetName("DELTA_HSP_LSP_T80");DELTA_HSP_LSP_T80->SetTitle("DELTA_HSP_LSP_T80");
  TGraphErrors* DELTA_HSP_LSP_T82 = (TGraphErrors*)f_82->Get("DELTA_HSP_LSP");DELTA_HSP_LSP_T82->SetName("DELTA_HSP_LSP_T82");DELTA_HSP_LSP_T82->SetTitle("DELTA_HSP_LSP_T82");
  TGraphErrors* DELTA_HSP_LSP_T84 = (TGraphErrors*)f_84->Get("DELTA_HSP_LSP");DELTA_HSP_LSP_T84->SetName("DELTA_HSP_LSP_T84");DELTA_HSP_LSP_T84->SetTitle("DELTA_HSP_LSP_T84");
  TGraphErrors* DELTA_HSP_LSP_T86 = (TGraphErrors*)f_86->Get("DELTA_HSP_LSP");DELTA_HSP_LSP_T86->SetName("DELTA_HSP_LSP_T86");DELTA_HSP_LSP_T86->SetTitle("DELTA_HSP_LSP_T86");
  TGraphErrors* DELTA_HSP_LSP_T88 = (TGraphErrors*)f_88->Get("DELTA_HSP_LSP");DELTA_HSP_LSP_T88->SetName("DELTA_HSP_LSP_T88");DELTA_HSP_LSP_T88->SetTitle("DELTA_HSP_LSP_T88");
  TGraphErrors* DELTA_HSP_LSP_T90 = (TGraphErrors*)f_90->Get("DELTA_HSP_LSP");DELTA_HSP_LSP_T90->SetName("DELTA_HSP_LSP_T90");DELTA_HSP_LSP_T90->SetTitle("DELTA_HSP_LSP_T90");
  TGraphErrors* DELTA_HSP_LSP_T78 = (TGraphErrors*)f_78->Get("DELTA_HSP_LSP");DELTA_HSP_LSP_T78->SetName("DELTA_HSP_LSP_T78");DELTA_HSP_LSP_T78->SetTitle("DELTA_HSP_LSP_T78");
  TGraphErrors* DELTA_HSP_LSP_T76 = (TGraphErrors*)f_76->Get("DELTA_HSP_LSP");DELTA_HSP_LSP_T76->SetName("DELTA_HSP_LSP_T76");DELTA_HSP_LSP_T76->SetTitle("DELTA_HSP_LSP_T76");
  TGraphErrors* DELTA_HSP_LSP_T74 = (TGraphErrors*)f_74->Get("DELTA_HSP_LSP");DELTA_HSP_LSP_T74->SetName("DELTA_HSP_LSP_T74");DELTA_HSP_LSP_T74->SetTitle("DELTA_HSP_LSP_T74");
  TGraphErrors* DELTA_HSP_LSP_T72 = (TGraphErrors*)f_72->Get("DELTA_HSP_LSP");DELTA_HSP_LSP_T72->SetName("DELTA_HSP_LSP_T72");DELTA_HSP_LSP_T72->SetTitle("DELTA_HSP_LSP_T72");
  TGraphErrors* DELTA_HSP_LSP_T73 = (TGraphErrors*)f_73->Get("DELTA_HSP_LSP");DELTA_HSP_LSP_T73->SetName("DELTA_HSP_LSP_T73");DELTA_HSP_LSP_T73->SetTitle("DELTA_HSP_LSP_T73");
  TGraphErrors* DELTA_HSP_LSP_T75 = (TGraphErrors*)f_75->Get("DELTA_HSP_LSP");DELTA_HSP_LSP_T75->SetName("DELTA_HSP_LSP_T75");DELTA_HSP_LSP_T75->SetTitle("DELTA_HSP_LSP_T75");
  TGraphErrors* DELTA_HSP_LSP_T77 = (TGraphErrors*)f_77->Get("DELTA_HSP_LSP");DELTA_HSP_LSP_T77->SetName("DELTA_HSP_LSP_T77");DELTA_HSP_LSP_T77->SetTitle("DELTA_HSP_LSP_T77");
  TGraphErrors* DELTA_HSP_LSP_T79 = (TGraphErrors*)f_79->Get("DELTA_HSP_LSP");DELTA_HSP_LSP_T79->SetName("DELTA_HSP_LSP_T79");DELTA_HSP_LSP_T79->SetTitle("DELTA_HSP_LSP_T79");
  TGraphErrors* DELTA_HSP_LSP_T81 = (TGraphErrors*)f_81->Get("DELTA_HSP_LSP");DELTA_HSP_LSP_T81->SetName("DELTA_HSP_LSP_T81");DELTA_HSP_LSP_T81->SetTitle("DELTA_HSP_LSP_T81");
  TGraphErrors* DELTA_HSP_LSP_T795 = (TGraphErrors*)f_795->Get("DELTA_HSP_LSP");DELTA_HSP_LSP_T795->SetName("DELTA_HSP_LSP_T795");DELTA_HSP_LSP_T795->SetTitle("DELTA_HSP_LSP_T795");

  // DELTA_HSP_LSP_T77->SetMarkerStyle(20);DELTA_HSP_LSP_T77->SetMarkerColor(1);         DELTA_HSP_LSP_T77->SetMarkerSize(2);
  DELTA_HSP_LSP_T80->SetMarkerStyle(20);DELTA_HSP_LSP_T80->SetMarkerColor(9);         DELTA_HSP_LSP_T80->SetMarkerSize(2);
  DELTA_HSP_LSP_T82->SetMarkerStyle(20);DELTA_HSP_LSP_T82->SetMarkerColor(2);         DELTA_HSP_LSP_T82->SetMarkerSize(2);
  DELTA_HSP_LSP_T84->SetMarkerStyle(20);DELTA_HSP_LSP_T84->SetMarkerColor(4);         DELTA_HSP_LSP_T84->SetMarkerSize(2);
  DELTA_HSP_LSP_T86->SetMarkerStyle(20);DELTA_HSP_LSP_T86->SetMarkerColor(8);         DELTA_HSP_LSP_T86->SetMarkerSize(2);
  DELTA_HSP_LSP_T88->SetMarkerStyle(20);DELTA_HSP_LSP_T88->SetMarkerColor(kOrange-2); DELTA_HSP_LSP_T88->SetMarkerSize(2);
  DELTA_HSP_LSP_T90->SetMarkerStyle(20);DELTA_HSP_LSP_T90->SetMarkerColor(11);        DELTA_HSP_LSP_T90->SetMarkerSize(2);
  DELTA_HSP_LSP_T78->SetMarkerStyle(20);DELTA_HSP_LSP_T78->SetMarkerColor(7);         DELTA_HSP_LSP_T78->SetMarkerSize(2);
  DELTA_HSP_LSP_T76->SetMarkerStyle(20);DELTA_HSP_LSP_T76->SetMarkerColor(6);         DELTA_HSP_LSP_T76->SetMarkerSize(2);
  DELTA_HSP_LSP_T74->SetMarkerStyle(20);DELTA_HSP_LSP_T74->SetMarkerColor(5);         DELTA_HSP_LSP_T74->SetMarkerSize(2);
  DELTA_HSP_LSP_T72->SetMarkerStyle(20);DELTA_HSP_LSP_T72->SetMarkerColor(12);        DELTA_HSP_LSP_T72->SetMarkerSize(2);
  DELTA_HSP_LSP_T73->SetMarkerStyle(20);DELTA_HSP_LSP_T73->SetMarkerColor(3);         DELTA_HSP_LSP_T73->SetMarkerSize(2);
  DELTA_HSP_LSP_T75->SetMarkerStyle(20);DELTA_HSP_LSP_T75->SetMarkerColor(13);        DELTA_HSP_LSP_T75->SetMarkerSize(2);
  DELTA_HSP_LSP_T77->SetMarkerStyle(20);DELTA_HSP_LSP_T77->SetMarkerColor(3);         DELTA_HSP_LSP_T77->SetMarkerSize(2);
  DELTA_HSP_LSP_T79->SetMarkerStyle(20);DELTA_HSP_LSP_T79->SetMarkerColor(13);        DELTA_HSP_LSP_T79->SetMarkerSize(2);
  DELTA_HSP_LSP_T81->SetMarkerStyle(20);DELTA_HSP_LSP_T81->SetMarkerColor(3);         DELTA_HSP_LSP_T81->SetMarkerSize(2);
  DELTA_HSP_LSP_T795->SetMarkerStyle(20);DELTA_HSP_LSP_T795->SetMarkerColor(13);        DELTA_HSP_LSP_T795->SetMarkerSize(2);





  
  TGraphErrors* SIGMA_H_B_0_T80 = (TGraphErrors*)f_80->Get("SIGMA_H_B_0");SIGMA_H_B_0_T80->SetName("SIGMA_H_B_0_T80");SIGMA_H_B_0_T80->SetTitle("SIGMA_H_B_0_T80");
  TGraphErrors* SIGMA_H_B_0_T82 = (TGraphErrors*)f_82->Get("SIGMA_H_B_0");SIGMA_H_B_0_T82->SetName("SIGMA_H_B_0_T82");SIGMA_H_B_0_T82->SetTitle("SIGMA_H_B_0_T82");
  TGraphErrors* SIGMA_H_B_0_T84 = (TGraphErrors*)f_84->Get("SIGMA_H_B_0");SIGMA_H_B_0_T84->SetName("SIGMA_H_B_0_T84");SIGMA_H_B_0_T84->SetTitle("SIGMA_H_B_0_T84");
  TGraphErrors* SIGMA_H_B_0_T86 = (TGraphErrors*)f_86->Get("SIGMA_H_B_0");SIGMA_H_B_0_T86->SetName("SIGMA_H_B_0_T86");SIGMA_H_B_0_T86->SetTitle("SIGMA_H_B_0_T86");
  TGraphErrors* SIGMA_H_B_0_T88 = (TGraphErrors*)f_88->Get("SIGMA_H_B_0");SIGMA_H_B_0_T88->SetName("SIGMA_H_B_0_T88");SIGMA_H_B_0_T88->SetTitle("SIGMA_H_B_0_T88");
  TGraphErrors* SIGMA_H_B_0_T90 = (TGraphErrors*)f_90->Get("SIGMA_H_B_0");SIGMA_H_B_0_T90->SetName("SIGMA_H_B_0_T90");SIGMA_H_B_0_T90->SetTitle("SIGMA_H_B_0_T90");
  TGraphErrors* SIGMA_H_B_0_T78 = (TGraphErrors*)f_78->Get("SIGMA_H_B_0");SIGMA_H_B_0_T78->SetName("SIGMA_H_B_0_T78");SIGMA_H_B_0_T78->SetTitle("SIGMA_H_B_0_T78");
  TGraphErrors* SIGMA_H_B_0_T76 = (TGraphErrors*)f_76->Get("SIGMA_H_B_0");SIGMA_H_B_0_T76->SetName("SIGMA_H_B_0_T76");SIGMA_H_B_0_T76->SetTitle("SIGMA_H_B_0_T76");
  TGraphErrors* SIGMA_H_B_0_T74 = (TGraphErrors*)f_74->Get("SIGMA_H_B_0");SIGMA_H_B_0_T74->SetName("SIGMA_H_B_0_T74");SIGMA_H_B_0_T74->SetTitle("SIGMA_H_B_0_T74");
  TGraphErrors* SIGMA_H_B_0_T72 = (TGraphErrors*)f_72->Get("SIGMA_H_B_0");SIGMA_H_B_0_T72->SetName("SIGMA_H_B_0_T72");SIGMA_H_B_0_T72->SetTitle("SIGMA_H_B_0_T72");
  TGraphErrors* SIGMA_H_B_0_T73 = (TGraphErrors*)f_73->Get("SIGMA_H_B_0");SIGMA_H_B_0_T73->SetName("SIGMA_H_B_0_T73");SIGMA_H_B_0_T73->SetTitle("SIGMA_H_B_0_T73");
  TGraphErrors* SIGMA_H_B_0_T75 = (TGraphErrors*)f_75->Get("SIGMA_H_B_0");SIGMA_H_B_0_T75->SetName("SIGMA_H_B_0_T75");SIGMA_H_B_0_T75->SetTitle("SIGMA_H_B_0_T75");
  TGraphErrors* SIGMA_H_B_0_T77 = (TGraphErrors*)f_77->Get("SIGMA_H_B_0");SIGMA_H_B_0_T77->SetName("SIGMA_H_B_0_T77");SIGMA_H_B_0_T77->SetTitle("SIGMA_H_B_0_T77");
  TGraphErrors* SIGMA_H_B_0_T79 = (TGraphErrors*)f_79->Get("SIGMA_H_B_0");SIGMA_H_B_0_T79->SetName("SIGMA_H_B_0_T79");SIGMA_H_B_0_T79->SetTitle("SIGMA_H_B_0_T79");
  TGraphErrors* SIGMA_H_B_0_T81 = (TGraphErrors*)f_81->Get("SIGMA_H_B_0");SIGMA_H_B_0_T81->SetName("SIGMA_H_B_0_T81");SIGMA_H_B_0_T81->SetTitle("SIGMA_H_B_0_T81");
  TGraphErrors* SIGMA_H_B_0_T795 = (TGraphErrors*)f_795->Get("SIGMA_H_B_0");SIGMA_H_B_0_T795->SetName("SIGMA_H_B_0_T795");SIGMA_H_B_0_T795->SetTitle("SIGMA_H_B_0_T795");

  // SIGMA_H_B_0_T77->SetMarkerStyle(20);SIGMA_H_B_0_T77->SetMarkerColor(1);         SIGMA_H_B_0_T77->SetMarkerSize(2);
  SIGMA_H_B_0_T80->SetMarkerStyle(20);SIGMA_H_B_0_T80->SetMarkerColor(9);         SIGMA_H_B_0_T80->SetMarkerSize(2);
  SIGMA_H_B_0_T82->SetMarkerStyle(20);SIGMA_H_B_0_T82->SetMarkerColor(2);         SIGMA_H_B_0_T82->SetMarkerSize(2);
  SIGMA_H_B_0_T84->SetMarkerStyle(20);SIGMA_H_B_0_T84->SetMarkerColor(4);         SIGMA_H_B_0_T84->SetMarkerSize(2);
  SIGMA_H_B_0_T86->SetMarkerStyle(20);SIGMA_H_B_0_T86->SetMarkerColor(8);         SIGMA_H_B_0_T86->SetMarkerSize(2);
  SIGMA_H_B_0_T88->SetMarkerStyle(20);SIGMA_H_B_0_T88->SetMarkerColor(kOrange-2); SIGMA_H_B_0_T88->SetMarkerSize(2);
  SIGMA_H_B_0_T90->SetMarkerStyle(20);SIGMA_H_B_0_T90->SetMarkerColor(11);        SIGMA_H_B_0_T90->SetMarkerSize(2);
  SIGMA_H_B_0_T78->SetMarkerStyle(20);SIGMA_H_B_0_T78->SetMarkerColor(7);         SIGMA_H_B_0_T78->SetMarkerSize(2);
  SIGMA_H_B_0_T76->SetMarkerStyle(20);SIGMA_H_B_0_T76->SetMarkerColor(6);         SIGMA_H_B_0_T76->SetMarkerSize(2);
  SIGMA_H_B_0_T74->SetMarkerStyle(20);SIGMA_H_B_0_T74->SetMarkerColor(5);         SIGMA_H_B_0_T74->SetMarkerSize(2);
  SIGMA_H_B_0_T72->SetMarkerStyle(20);SIGMA_H_B_0_T72->SetMarkerColor(12);        SIGMA_H_B_0_T72->SetMarkerSize(2);
  SIGMA_H_B_0_T73->SetMarkerStyle(20);SIGMA_H_B_0_T73->SetMarkerColor(3);         SIGMA_H_B_0_T73->SetMarkerSize(2);
  SIGMA_H_B_0_T75->SetMarkerStyle(20);SIGMA_H_B_0_T75->SetMarkerColor(13);        SIGMA_H_B_0_T75->SetMarkerSize(2);
  SIGMA_H_B_0_T77->SetMarkerStyle(20);SIGMA_H_B_0_T77->SetMarkerColor(3);         SIGMA_H_B_0_T77->SetMarkerSize(2);
  SIGMA_H_B_0_T79->SetMarkerStyle(20);SIGMA_H_B_0_T79->SetMarkerColor(13);        SIGMA_H_B_0_T79->SetMarkerSize(2);
  SIGMA_H_B_0_T81->SetMarkerStyle(20);SIGMA_H_B_0_T81->SetMarkerColor(3);         SIGMA_H_B_0_T81->SetMarkerSize(2);
  SIGMA_H_B_0_T795->SetMarkerStyle(20);SIGMA_H_B_0_T795->SetMarkerColor(13);        SIGMA_H_B_0_T795->SetMarkerSize(2);


  

  
  TGraphErrors* DELTA_H_B_0_T80 = (TGraphErrors*)f_80->Get("DELTA_H_B_0");DELTA_H_B_0_T80->SetName("DELTA_H_B_0_T80");DELTA_H_B_0_T80->SetTitle("DELTA_H_B_0_T80");
  TGraphErrors* DELTA_H_B_0_T82 = (TGraphErrors*)f_82->Get("DELTA_H_B_0");DELTA_H_B_0_T82->SetName("DELTA_H_B_0_T82");DELTA_H_B_0_T82->SetTitle("DELTA_H_B_0_T82");
  TGraphErrors* DELTA_H_B_0_T84 = (TGraphErrors*)f_84->Get("DELTA_H_B_0");DELTA_H_B_0_T84->SetName("DELTA_H_B_0_T84");DELTA_H_B_0_T84->SetTitle("DELTA_H_B_0_T84");
  TGraphErrors* DELTA_H_B_0_T86 = (TGraphErrors*)f_86->Get("DELTA_H_B_0");DELTA_H_B_0_T86->SetName("DELTA_H_B_0_T86");DELTA_H_B_0_T86->SetTitle("DELTA_H_B_0_T86");
  TGraphErrors* DELTA_H_B_0_T88 = (TGraphErrors*)f_88->Get("DELTA_H_B_0");DELTA_H_B_0_T88->SetName("DELTA_H_B_0_T88");DELTA_H_B_0_T88->SetTitle("DELTA_H_B_0_T88");
  TGraphErrors* DELTA_H_B_0_T90 = (TGraphErrors*)f_90->Get("DELTA_H_B_0");DELTA_H_B_0_T90->SetName("DELTA_H_B_0_T90");DELTA_H_B_0_T90->SetTitle("DELTA_H_B_0_T90");
  TGraphErrors* DELTA_H_B_0_T78 = (TGraphErrors*)f_78->Get("DELTA_H_B_0");DELTA_H_B_0_T78->SetName("DELTA_H_B_0_T78");DELTA_H_B_0_T78->SetTitle("DELTA_H_B_0_T78");
  TGraphErrors* DELTA_H_B_0_T76 = (TGraphErrors*)f_76->Get("DELTA_H_B_0");DELTA_H_B_0_T76->SetName("DELTA_H_B_0_T76");DELTA_H_B_0_T76->SetTitle("DELTA_H_B_0_T76");
  TGraphErrors* DELTA_H_B_0_T74 = (TGraphErrors*)f_74->Get("DELTA_H_B_0");DELTA_H_B_0_T74->SetName("DELTA_H_B_0_T74");DELTA_H_B_0_T74->SetTitle("DELTA_H_B_0_T74");
  TGraphErrors* DELTA_H_B_0_T72 = (TGraphErrors*)f_72->Get("DELTA_H_B_0");DELTA_H_B_0_T72->SetName("DELTA_H_B_0_T72");DELTA_H_B_0_T72->SetTitle("DELTA_H_B_0_T72");
  TGraphErrors* DELTA_H_B_0_T73 = (TGraphErrors*)f_73->Get("DELTA_H_B_0");DELTA_H_B_0_T73->SetName("DELTA_H_B_0_T73");DELTA_H_B_0_T73->SetTitle("DELTA_H_B_0_T73");
  TGraphErrors* DELTA_H_B_0_T75 = (TGraphErrors*)f_75->Get("DELTA_H_B_0");DELTA_H_B_0_T75->SetName("DELTA_H_B_0_T75");DELTA_H_B_0_T75->SetTitle("DELTA_H_B_0_T75");
  TGraphErrors* DELTA_H_B_0_T77 = (TGraphErrors*)f_77->Get("DELTA_H_B_0");DELTA_H_B_0_T77->SetName("DELTA_H_B_0_T77");DELTA_H_B_0_T77->SetTitle("DELTA_H_B_0_T77");
  TGraphErrors* DELTA_H_B_0_T79 = (TGraphErrors*)f_79->Get("DELTA_H_B_0");DELTA_H_B_0_T79->SetName("DELTA_H_B_0_T79");DELTA_H_B_0_T79->SetTitle("DELTA_H_B_0_T79");
  TGraphErrors* DELTA_H_B_0_T81 = (TGraphErrors*)f_81->Get("DELTA_H_B_0");DELTA_H_B_0_T81->SetName("DELTA_H_B_0_T81");DELTA_H_B_0_T81->SetTitle("DELTA_H_B_0_T81");
  TGraphErrors* DELTA_H_B_0_T795 = (TGraphErrors*)f_795->Get("DELTA_H_B_0");DELTA_H_B_0_T795->SetName("DELTA_H_B_0_T795");DELTA_H_B_0_T795->SetTitle("DELTA_H_B_0_T795");

  // DELTA_H_B_0_T77->SetMarkerStyle(20);DELTA_H_B_0_T77->SetMarkerColor(1);         DELTA_H_B_0_T77->SetMarkerSize(2);
  DELTA_H_B_0_T80->SetMarkerStyle(20);DELTA_H_B_0_T80->SetMarkerColor(9);         DELTA_H_B_0_T80->SetMarkerSize(2);
  DELTA_H_B_0_T82->SetMarkerStyle(20);DELTA_H_B_0_T82->SetMarkerColor(2);         DELTA_H_B_0_T82->SetMarkerSize(2);
  DELTA_H_B_0_T84->SetMarkerStyle(20);DELTA_H_B_0_T84->SetMarkerColor(4);         DELTA_H_B_0_T84->SetMarkerSize(2);
  DELTA_H_B_0_T86->SetMarkerStyle(20);DELTA_H_B_0_T86->SetMarkerColor(8);         DELTA_H_B_0_T86->SetMarkerSize(2);
  DELTA_H_B_0_T88->SetMarkerStyle(20);DELTA_H_B_0_T88->SetMarkerColor(kOrange-2); DELTA_H_B_0_T88->SetMarkerSize(2);
  DELTA_H_B_0_T90->SetMarkerStyle(20);DELTA_H_B_0_T90->SetMarkerColor(11);        DELTA_H_B_0_T90->SetMarkerSize(2);
  DELTA_H_B_0_T78->SetMarkerStyle(20);DELTA_H_B_0_T78->SetMarkerColor(7);         DELTA_H_B_0_T78->SetMarkerSize(2);
  DELTA_H_B_0_T76->SetMarkerStyle(20);DELTA_H_B_0_T76->SetMarkerColor(6);         DELTA_H_B_0_T76->SetMarkerSize(2);
  DELTA_H_B_0_T74->SetMarkerStyle(20);DELTA_H_B_0_T74->SetMarkerColor(5);         DELTA_H_B_0_T74->SetMarkerSize(2);
  DELTA_H_B_0_T72->SetMarkerStyle(20);DELTA_H_B_0_T72->SetMarkerColor(12);        DELTA_H_B_0_T72->SetMarkerSize(2);
  DELTA_H_B_0_T73->SetMarkerStyle(20);DELTA_H_B_0_T73->SetMarkerColor(3);         DELTA_H_B_0_T73->SetMarkerSize(2);
  DELTA_H_B_0_T75->SetMarkerStyle(20);DELTA_H_B_0_T75->SetMarkerColor(13);        DELTA_H_B_0_T75->SetMarkerSize(2);
  DELTA_H_B_0_T77->SetMarkerStyle(20);DELTA_H_B_0_T77->SetMarkerColor(3);         DELTA_H_B_0_T77->SetMarkerSize(2);
  DELTA_H_B_0_T79->SetMarkerStyle(20);DELTA_H_B_0_T79->SetMarkerColor(13);        DELTA_H_B_0_T79->SetMarkerSize(2);
  DELTA_H_B_0_T81->SetMarkerStyle(20);DELTA_H_B_0_T81->SetMarkerColor(3);         DELTA_H_B_0_T81->SetMarkerSize(2);
  DELTA_H_B_0_T795->SetMarkerStyle(20);DELTA_H_B_0_T795->SetMarkerColor(13);        DELTA_H_B_0_T795->SetMarkerSize(2);
  ////////////////////////////////////////////////
  /// Single pixel graphs //////
  //////////////////////////////


  Int_t n=4;
  
  float TUNES[16];
   
  //TUNES[8]=77;
  TUNES[15]=81;
  TUNES[14]=79.5;
  TUNES[13]=79;
  TUNES[12]=77;
  TUNES[11]=75;
  TUNES[10]=72;
  TUNES[9]=73;
  TUNES[8]=74;
  TUNES[7]=76;
  TUNES[6]=78;
  TUNES[5]=80;
  TUNES[4]=82;
  TUNES[3]=84;
  TUNES[2]=86;
  TUNES[1]=88;
  TUNES[0]=90;
  float err_TUNES[16];
  err_TUNES[0]=0;
  err_TUNES[1]=0;
  err_TUNES[2]=0;
  err_TUNES[3]=0;
  err_TUNES[4]=0;
  err_TUNES[5]=0;
  err_TUNES[6]=0;
  err_TUNES[7]=0;
  err_TUNES[8]=0;
  err_TUNES[9]=0;
  err_TUNES[10]=0;
  err_TUNES[11]=0;
  err_TUNES[12]=0;
  err_TUNES[13]=0;
  err_TUNES[14]=0;
  err_TUNES[15]=0;


  Double_t x_au;
  Double_t FRAC_CH_au;
  float FRAC_CH[16];
  float err_FRAC_CH[16];
  
  //FRAC_SIG_T77->GetPoint(n,x_au,FRAC_CH_au); FRAC_CH[8]=FRAC_CH_au;
  FRAC_SIG_T81->GetPoint(n,x_au,FRAC_CH_au); FRAC_CH[15]=FRAC_CH_au;
  FRAC_SIG_T795->GetPoint(n,x_au,FRAC_CH_au); FRAC_CH[14]=FRAC_CH_au;
  FRAC_SIG_T79->GetPoint(n,x_au,FRAC_CH_au); FRAC_CH[13]=FRAC_CH_au;
  FRAC_SIG_T77->GetPoint(n,x_au,FRAC_CH_au); FRAC_CH[12]=FRAC_CH_au;
  FRAC_SIG_T75->GetPoint(n,x_au,FRAC_CH_au); FRAC_CH[11]=FRAC_CH_au;
  FRAC_SIG_T72->GetPoint(n,x_au,FRAC_CH_au); FRAC_CH[10]=FRAC_CH_au;
  FRAC_SIG_T73->GetPoint(n,x_au,FRAC_CH_au); FRAC_CH[9]=FRAC_CH_au;
  FRAC_SIG_T74->GetPoint(n,x_au,FRAC_CH_au); FRAC_CH[8]=FRAC_CH_au;
  FRAC_SIG_T76->GetPoint(n,x_au,FRAC_CH_au); FRAC_CH[7]=FRAC_CH_au;
  FRAC_SIG_T78->GetPoint(n,x_au,FRAC_CH_au); FRAC_CH[6]=FRAC_CH_au;
  FRAC_SIG_T80->GetPoint(n,x_au,FRAC_CH_au); FRAC_CH[5]=FRAC_CH_au;
  FRAC_SIG_T82->GetPoint(n,x_au,FRAC_CH_au); FRAC_CH[4]=FRAC_CH_au;
  FRAC_SIG_T84->GetPoint(n,x_au,FRAC_CH_au); FRAC_CH[3]=FRAC_CH_au;
  FRAC_SIG_T86->GetPoint(n,x_au,FRAC_CH_au); FRAC_CH[2]=FRAC_CH_au;
  FRAC_SIG_T88->GetPoint(n,x_au,FRAC_CH_au); FRAC_CH[1]=FRAC_CH_au;
  FRAC_SIG_T90->GetPoint(n,x_au,FRAC_CH_au); FRAC_CH[0]=FRAC_CH_au;

  
  err_FRAC_CH[15]=(float)FRAC_SIG_T81->GetErrorXhigh(n);
  err_FRAC_CH[14]=(float)FRAC_SIG_T795->GetErrorXhigh(n);
  err_FRAC_CH[13]=(float)FRAC_SIG_T79->GetErrorXhigh(n);
  err_FRAC_CH[12]=(float)FRAC_SIG_T77->GetErrorXhigh(n);
  err_FRAC_CH[11]=(float)FRAC_SIG_T75->GetErrorXhigh(n);
  err_FRAC_CH[10]=(float)FRAC_SIG_T72->GetErrorXhigh(n);
  err_FRAC_CH[9]=(float)FRAC_SIG_T73->GetErrorXhigh(n);
  err_FRAC_CH[8]=(float)FRAC_SIG_T74->GetErrorXhigh(n);
  err_FRAC_CH[7]=(float)FRAC_SIG_T76->GetErrorXhigh(n);
  err_FRAC_CH[6]=(float)FRAC_SIG_T78->GetErrorXhigh(n);
  err_FRAC_CH[5]=(float)FRAC_SIG_T80->GetErrorXhigh(n);
  err_FRAC_CH[4]=(float)FRAC_SIG_T82->GetErrorXhigh(n);
  err_FRAC_CH[3]=(float)FRAC_SIG_T84->GetErrorXhigh(n);
  err_FRAC_CH[2]=(float)FRAC_SIG_T86->GetErrorXhigh(n);
  err_FRAC_CH[1]=(float)FRAC_SIG_T88->GetErrorXhigh(n);
  err_FRAC_CH[0]=(float)FRAC_SIG_T90->GetErrorXhigh(n);

  
  TGraphErrors *TGE_FRAC_CH = new TGraphErrors(16,TUNES,FRAC_CH,err_TUNES,err_FRAC_CH);TGE_FRAC_CH->SetName("FRAC_CH 4"); TGE_FRAC_CH->SetTitle(Form("FRAC_CH %d",n)); 
  TGE_FRAC_CH->SetMarkerStyle(20);TGE_FRAC_CH->SetMarkerColor(4);         TGE_FRAC_CH->SetMarkerSize(2);
  TCanvas* c_FRAC_CH = new TCanvas("c_FRAC_CH","c_FRAC_CH");
  TGE_FRAC_CH->Draw("AP");
  TGE_FRAC_CH->GetXaxis()->SetTitle("Laser Tune");
  TGE_FRAC_CH->GetYaxis()->SetTitle("F_{Sig}");
  gPad->BuildLegend();
  gPad->Update();
  gPad->SetGridy();
  //TGE_FRAC_CH->Fit("pol2");










  
  Double_t T0_au;
  float T0[16];
  float err_T0[16];

  
  T0_T81->GetPoint(n,x_au,T0_au); T0[15]=T0_au;
  T0_T795->GetPoint(n,x_au,T0_au); T0[14]=T0_au;   
  T0_T79->GetPoint(n,x_au,T0_au); T0[13]=T0_au;
  T0_T77->GetPoint(n,x_au,T0_au); T0[12]=T0_au;
  T0_T75->GetPoint(n,x_au,T0_au); T0[11]=T0_au;
  T0_T72->GetPoint(n,x_au,T0_au); T0[10]=T0_au;
  T0_T73->GetPoint(n,x_au,T0_au); T0[9]=T0_au;
  T0_T74->GetPoint(n,x_au,T0_au); T0[8]=T0_au;
  T0_T76->GetPoint(n,x_au,T0_au); T0[7]=T0_au;
  T0_T78->GetPoint(n,x_au,T0_au); T0[6]=T0_au;
  T0_T80->GetPoint(n,x_au,T0_au); T0[5]=T0_au;
  T0_T82->GetPoint(n,x_au,T0_au); T0[4]=T0_au;
  T0_T84->GetPoint(n,x_au,T0_au); T0[3]=T0_au;
  T0_T86->GetPoint(n,x_au,T0_au); T0[2]=T0_au;
  T0_T88->GetPoint(n,x_au,T0_au); T0[1]=T0_au;
  T0_T90->GetPoint(n,x_au,T0_au); T0[0]=T0_au;


  err_T0[15]=T0_T81->GetErrorX(n);
  err_T0[14]=T0_T795->GetErrorX(n);
  err_T0[13]=T0_T79->GetErrorX(n);
  err_T0[12]=T0_T77->GetErrorX(n);
  err_T0[11]=T0_T75->GetErrorX(n);
  err_T0[10]=T0_T72->GetErrorX(n);
  err_T0[9]=T0_T73->GetErrorX(n);
  err_T0[8]=T0_T74->GetErrorX(n);
  err_T0[7]=T0_T76->GetErrorX(n);
  err_T0[6]=T0_T78->GetErrorX(n);
  err_T0[5]=T0_T80->GetErrorX(n);
  err_T0[4]=T0_T82->GetErrorX(n);
  err_T0[3]=T0_T84->GetErrorX(n);
  err_T0[2]=T0_T86->GetErrorX(n);
  err_T0[1]=T0_T88->GetErrorX(n);
  err_T0[0]=T0_T90->GetErrorX(n);

  
  TGraphErrors *TGE_T0 = new TGraphErrors(16,TUNES,T0,err_TUNES,err_T0);TGE_T0->SetName("T0 che 4"); TGE_T0->SetTitle(Form("T_{0} %d",n)); 
  TGE_T0->SetMarkerStyle(20);TGE_T0->SetMarkerColor(4);         TGE_T0->SetMarkerSize(2);
  TCanvas* c_T0_CH = new TCanvas("c_T0_CH","c_T0_CH");
  TGE_T0->Draw("AP");
  TGE_T0->GetXaxis()->SetTitle("Laser Tune");
  TGE_T0->GetYaxis()->SetTitle("T_{0} [ns]");
  gPad->BuildLegend();
  gPad->Update();
  gPad->SetGridy();

  


  //////////////////////




  
  Double_t DELTA_LSP_T0_au;
  float DELTA_LSP_T0[16];
  float err_DELTA_LSP_T0[16];

  
  DELTA_LSP_T0_T81->GetPoint(n,x_au,DELTA_LSP_T0_au); DELTA_LSP_T0[15]=DELTA_LSP_T0_au;
  DELTA_LSP_T0_T795->GetPoint(n,x_au,DELTA_LSP_T0_au); DELTA_LSP_T0[14]=DELTA_LSP_T0_au;   
  DELTA_LSP_T0_T79->GetPoint(n,x_au,DELTA_LSP_T0_au); DELTA_LSP_T0[13]=DELTA_LSP_T0_au;
  DELTA_LSP_T0_T77->GetPoint(n,x_au,DELTA_LSP_T0_au); DELTA_LSP_T0[12]=DELTA_LSP_T0_au;
  DELTA_LSP_T0_T75->GetPoint(n,x_au,DELTA_LSP_T0_au); DELTA_LSP_T0[11]=DELTA_LSP_T0_au;
  DELTA_LSP_T0_T72->GetPoint(n,x_au,DELTA_LSP_T0_au); DELTA_LSP_T0[10]=DELTA_LSP_T0_au;
  DELTA_LSP_T0_T73->GetPoint(n,x_au,DELTA_LSP_T0_au); DELTA_LSP_T0[9]=DELTA_LSP_T0_au;
  DELTA_LSP_T0_T74->GetPoint(n,x_au,DELTA_LSP_T0_au); DELTA_LSP_T0[8]=DELTA_LSP_T0_au;
  DELTA_LSP_T0_T76->GetPoint(n,x_au,DELTA_LSP_T0_au); DELTA_LSP_T0[7]=DELTA_LSP_T0_au;
  DELTA_LSP_T0_T78->GetPoint(n,x_au,DELTA_LSP_T0_au); DELTA_LSP_T0[6]=DELTA_LSP_T0_au;
  DELTA_LSP_T0_T80->GetPoint(n,x_au,DELTA_LSP_T0_au); DELTA_LSP_T0[5]=DELTA_LSP_T0_au;
  DELTA_LSP_T0_T82->GetPoint(n,x_au,DELTA_LSP_T0_au); DELTA_LSP_T0[4]=DELTA_LSP_T0_au;
  DELTA_LSP_T0_T84->GetPoint(n,x_au,DELTA_LSP_T0_au); DELTA_LSP_T0[3]=DELTA_LSP_T0_au;
  DELTA_LSP_T0_T86->GetPoint(n,x_au,DELTA_LSP_T0_au); DELTA_LSP_T0[2]=DELTA_LSP_T0_au;
  DELTA_LSP_T0_T88->GetPoint(n,x_au,DELTA_LSP_T0_au); DELTA_LSP_T0[1]=DELTA_LSP_T0_au;
  DELTA_LSP_T0_T90->GetPoint(n,x_au,DELTA_LSP_T0_au); DELTA_LSP_T0[0]=DELTA_LSP_T0_au;


  err_DELTA_LSP_T0[15]=DELTA_LSP_T0_T81->GetErrorX(n);
  err_DELTA_LSP_T0[14]=DELTA_LSP_T0_T795->GetErrorX(n);
  err_DELTA_LSP_T0[13]=DELTA_LSP_T0_T79->GetErrorX(n);
  err_DELTA_LSP_T0[12]=DELTA_LSP_T0_T77->GetErrorX(n);
  err_DELTA_LSP_T0[11]=DELTA_LSP_T0_T75->GetErrorX(n);
  err_DELTA_LSP_T0[10]=DELTA_LSP_T0_T72->GetErrorX(n);
  err_DELTA_LSP_T0[9]=DELTA_LSP_T0_T73->GetErrorX(n);
  err_DELTA_LSP_T0[8]=DELTA_LSP_T0_T74->GetErrorX(n);
  err_DELTA_LSP_T0[7]=DELTA_LSP_T0_T76->GetErrorX(n);
  err_DELTA_LSP_T0[6]=DELTA_LSP_T0_T78->GetErrorX(n);
  err_DELTA_LSP_T0[5]=DELTA_LSP_T0_T80->GetErrorX(n);
  err_DELTA_LSP_T0[4]=DELTA_LSP_T0_T82->GetErrorX(n);
  err_DELTA_LSP_T0[3]=DELTA_LSP_T0_T84->GetErrorX(n);
  err_DELTA_LSP_T0[2]=DELTA_LSP_T0_T86->GetErrorX(n);
  err_DELTA_LSP_T0[1]=DELTA_LSP_T0_T88->GetErrorX(n);
  err_DELTA_LSP_T0[0]=DELTA_LSP_T0_T90->GetErrorX(n);

  
  TGraphErrors *TGE_DELTA_LSP_T0 = new TGraphErrors(16,TUNES,DELTA_LSP_T0,err_TUNES,err_DELTA_LSP_T0);TGE_DELTA_LSP_T0->SetName("DELTA_LSP_T0 che 4"); TGE_DELTA_LSP_T0->SetTitle(Form("T_{L}^{SP}-T_{0} ch %d",n)); 
  TGE_DELTA_LSP_T0->SetMarkerStyle(20);TGE_DELTA_LSP_T0->SetMarkerColor(4);         TGE_DELTA_LSP_T0->SetMarkerSize(2);
  TCanvas* c_DELTA_LSP_T0_CH = new TCanvas("c_DELTA_LSP_T0_CH","c_DELTA_LSP_T0_CH");
  TGE_DELTA_LSP_T0->Draw("AP");
  TGE_DELTA_LSP_T0->GetXaxis()->SetTitle("Laser Tune");
  TGE_DELTA_LSP_T0->GetYaxis()->SetTitle("T_{L}^{SP}-T_{0} [ns]");
  gPad->BuildLegend();
  gPad->Update();
  gPad->SetGridy();






  
  Double_t DELTA_HSP_LSP_au;
  float DELTA_HSP_LSP[16];
  float err_DELTA_HSP_LSP[16];

  
  DELTA_HSP_LSP_T81->GetPoint(n,x_au,DELTA_HSP_LSP_au); DELTA_HSP_LSP[15]=DELTA_HSP_LSP_au;
  DELTA_HSP_LSP_T795->GetPoint(n,x_au,DELTA_HSP_LSP_au); DELTA_HSP_LSP[14]=DELTA_HSP_LSP_au;   
  DELTA_HSP_LSP_T79->GetPoint(n,x_au,DELTA_HSP_LSP_au); DELTA_HSP_LSP[13]=DELTA_HSP_LSP_au;
  DELTA_HSP_LSP_T77->GetPoint(n,x_au,DELTA_HSP_LSP_au); DELTA_HSP_LSP[12]=DELTA_HSP_LSP_au;
  DELTA_HSP_LSP_T75->GetPoint(n,x_au,DELTA_HSP_LSP_au); DELTA_HSP_LSP[11]=DELTA_HSP_LSP_au;
  DELTA_HSP_LSP_T72->GetPoint(n,x_au,DELTA_HSP_LSP_au); DELTA_HSP_LSP[10]=DELTA_HSP_LSP_au;
  DELTA_HSP_LSP_T73->GetPoint(n,x_au,DELTA_HSP_LSP_au); DELTA_HSP_LSP[9]=DELTA_HSP_LSP_au;
  DELTA_HSP_LSP_T74->GetPoint(n,x_au,DELTA_HSP_LSP_au); DELTA_HSP_LSP[8]=DELTA_HSP_LSP_au;
  DELTA_HSP_LSP_T76->GetPoint(n,x_au,DELTA_HSP_LSP_au); DELTA_HSP_LSP[7]=DELTA_HSP_LSP_au;
  DELTA_HSP_LSP_T78->GetPoint(n,x_au,DELTA_HSP_LSP_au); DELTA_HSP_LSP[6]=DELTA_HSP_LSP_au;
  DELTA_HSP_LSP_T80->GetPoint(n,x_au,DELTA_HSP_LSP_au); DELTA_HSP_LSP[5]=DELTA_HSP_LSP_au;
  DELTA_HSP_LSP_T82->GetPoint(n,x_au,DELTA_HSP_LSP_au); DELTA_HSP_LSP[4]=DELTA_HSP_LSP_au;
  DELTA_HSP_LSP_T84->GetPoint(n,x_au,DELTA_HSP_LSP_au); DELTA_HSP_LSP[3]=DELTA_HSP_LSP_au;
  DELTA_HSP_LSP_T86->GetPoint(n,x_au,DELTA_HSP_LSP_au); DELTA_HSP_LSP[2]=DELTA_HSP_LSP_au;
  DELTA_HSP_LSP_T88->GetPoint(n,x_au,DELTA_HSP_LSP_au); DELTA_HSP_LSP[1]=DELTA_HSP_LSP_au;
  DELTA_HSP_LSP_T90->GetPoint(n,x_au,DELTA_HSP_LSP_au); DELTA_HSP_LSP[0]=DELTA_HSP_LSP_au;


  err_DELTA_HSP_LSP[15]=DELTA_HSP_LSP_T81->GetErrorX(n);
  err_DELTA_HSP_LSP[14]=DELTA_HSP_LSP_T795->GetErrorX(n);
  err_DELTA_HSP_LSP[13]=DELTA_HSP_LSP_T79->GetErrorX(n);
  err_DELTA_HSP_LSP[12]=DELTA_HSP_LSP_T77->GetErrorX(n);
  err_DELTA_HSP_LSP[11]=DELTA_HSP_LSP_T75->GetErrorX(n);
  err_DELTA_HSP_LSP[10]=DELTA_HSP_LSP_T72->GetErrorX(n);
  err_DELTA_HSP_LSP[9]=DELTA_HSP_LSP_T73->GetErrorX(n);
  err_DELTA_HSP_LSP[8]=DELTA_HSP_LSP_T74->GetErrorX(n);
  err_DELTA_HSP_LSP[7]=DELTA_HSP_LSP_T76->GetErrorX(n);
  err_DELTA_HSP_LSP[6]=DELTA_HSP_LSP_T78->GetErrorX(n);
  err_DELTA_HSP_LSP[5]=DELTA_HSP_LSP_T80->GetErrorX(n);
  err_DELTA_HSP_LSP[4]=DELTA_HSP_LSP_T82->GetErrorX(n);
  err_DELTA_HSP_LSP[3]=DELTA_HSP_LSP_T84->GetErrorX(n);
  err_DELTA_HSP_LSP[2]=DELTA_HSP_LSP_T86->GetErrorX(n);
  err_DELTA_HSP_LSP[1]=DELTA_HSP_LSP_T88->GetErrorX(n);
  err_DELTA_HSP_LSP[0]=DELTA_HSP_LSP_T90->GetErrorX(n);

  
  TGraphErrors *TGE_DELTA_HSP_LSP = new TGraphErrors(16,TUNES,DELTA_HSP_LSP,err_TUNES,err_DELTA_HSP_LSP);TGE_DELTA_HSP_LSP->SetName("DELTA_HSP_LSP che 4"); TGE_DELTA_HSP_LSP->SetTitle(Form("T_{H}^{SP}-T_{L}^{SP} ch %d",n)); 
  TGE_DELTA_HSP_LSP->SetMarkerStyle(20);TGE_DELTA_HSP_LSP->SetMarkerColor(4);         TGE_DELTA_HSP_LSP->SetMarkerSize(2);
  TCanvas* c_DELTA_HSP_LSP_CH = new TCanvas("c_DELTA_HSP_LSP_CH","c_DELTA_HSP_LSP_CH");
  TGE_DELTA_HSP_LSP->Draw("AP");
  TGE_DELTA_HSP_LSP->GetXaxis()->SetTitle("Laser Tune");
  TGE_DELTA_HSP_LSP->GetYaxis()->SetTitle("T_{H}^{SP}-T_{L}^{SP} [ns]");
  gPad->BuildLegend();
  gPad->Update();
  gPad->SetGridy();


  
  
  Double_t SIGMA_H_B_0_au;
  float SIGMA_H_B_0[16];
  float err_SIGMA_H_B_0[16];

  
  SIGMA_H_B_0_T81->GetPoint(n,x_au,SIGMA_H_B_0_au); SIGMA_H_B_0[15]=SIGMA_H_B_0_au;
  SIGMA_H_B_0_T795->GetPoint(n,x_au,SIGMA_H_B_0_au); SIGMA_H_B_0[14]=SIGMA_H_B_0_au;   
  SIGMA_H_B_0_T79->GetPoint(n,x_au,SIGMA_H_B_0_au); SIGMA_H_B_0[13]=SIGMA_H_B_0_au;
  SIGMA_H_B_0_T77->GetPoint(n,x_au,SIGMA_H_B_0_au); SIGMA_H_B_0[12]=SIGMA_H_B_0_au;
  SIGMA_H_B_0_T75->GetPoint(n,x_au,SIGMA_H_B_0_au); SIGMA_H_B_0[11]=SIGMA_H_B_0_au;
  SIGMA_H_B_0_T72->GetPoint(n,x_au,SIGMA_H_B_0_au); SIGMA_H_B_0[10]=SIGMA_H_B_0_au;
  SIGMA_H_B_0_T73->GetPoint(n,x_au,SIGMA_H_B_0_au); SIGMA_H_B_0[9]=SIGMA_H_B_0_au;
  SIGMA_H_B_0_T74->GetPoint(n,x_au,SIGMA_H_B_0_au); SIGMA_H_B_0[8]=SIGMA_H_B_0_au;
  SIGMA_H_B_0_T76->GetPoint(n,x_au,SIGMA_H_B_0_au); SIGMA_H_B_0[7]=SIGMA_H_B_0_au;
  SIGMA_H_B_0_T78->GetPoint(n,x_au,SIGMA_H_B_0_au); SIGMA_H_B_0[6]=SIGMA_H_B_0_au;
  SIGMA_H_B_0_T80->GetPoint(n,x_au,SIGMA_H_B_0_au); SIGMA_H_B_0[5]=SIGMA_H_B_0_au;
  SIGMA_H_B_0_T82->GetPoint(n,x_au,SIGMA_H_B_0_au); SIGMA_H_B_0[4]=SIGMA_H_B_0_au;
  SIGMA_H_B_0_T84->GetPoint(n,x_au,SIGMA_H_B_0_au); SIGMA_H_B_0[3]=SIGMA_H_B_0_au;
  SIGMA_H_B_0_T86->GetPoint(n,x_au,SIGMA_H_B_0_au); SIGMA_H_B_0[2]=SIGMA_H_B_0_au;
  SIGMA_H_B_0_T88->GetPoint(n,x_au,SIGMA_H_B_0_au); SIGMA_H_B_0[1]=SIGMA_H_B_0_au;
  SIGMA_H_B_0_T90->GetPoint(n,x_au,SIGMA_H_B_0_au); SIGMA_H_B_0[0]=SIGMA_H_B_0_au;


  err_SIGMA_H_B_0[15]=SIGMA_H_B_0_T81->GetErrorX(n);
  err_SIGMA_H_B_0[14]=SIGMA_H_B_0_T795->GetErrorX(n);
  err_SIGMA_H_B_0[13]=SIGMA_H_B_0_T79->GetErrorX(n);
  err_SIGMA_H_B_0[12]=SIGMA_H_B_0_T77->GetErrorX(n);
  err_SIGMA_H_B_0[11]=SIGMA_H_B_0_T75->GetErrorX(n);
  err_SIGMA_H_B_0[10]=SIGMA_H_B_0_T72->GetErrorX(n);
  err_SIGMA_H_B_0[9]=SIGMA_H_B_0_T73->GetErrorX(n);
  err_SIGMA_H_B_0[8]=SIGMA_H_B_0_T74->GetErrorX(n);
  err_SIGMA_H_B_0[7]=SIGMA_H_B_0_T76->GetErrorX(n);
  err_SIGMA_H_B_0[6]=SIGMA_H_B_0_T78->GetErrorX(n);
  err_SIGMA_H_B_0[5]=SIGMA_H_B_0_T80->GetErrorX(n);
  err_SIGMA_H_B_0[4]=SIGMA_H_B_0_T82->GetErrorX(n);
  err_SIGMA_H_B_0[3]=SIGMA_H_B_0_T84->GetErrorX(n);
  err_SIGMA_H_B_0[2]=SIGMA_H_B_0_T86->GetErrorX(n);
  err_SIGMA_H_B_0[1]=SIGMA_H_B_0_T88->GetErrorX(n);
  err_SIGMA_H_B_0[0]=SIGMA_H_B_0_T90->GetErrorX(n);

  
  TGraphErrors *TGE_SIGMA_H_B_0 = new TGraphErrors(16,TUNES,SIGMA_H_B_0,err_TUNES,err_SIGMA_H_B_0);TGE_SIGMA_H_B_0->SetName("SIGMA_H_B_0 che 4"); TGE_SIGMA_H_B_0->SetTitle(Form("#delta T_{0} ch %d",n)); 
  TGE_SIGMA_H_B_0->SetMarkerStyle(20);TGE_SIGMA_H_B_0->SetMarkerColor(4);         TGE_SIGMA_H_B_0->SetMarkerSize(2);
  TCanvas* c_SIGMA_H_B_0_CH = new TCanvas("c_SIGMA_H_B_0_CH","c_SIGMA_H_B_0_CH");
  TGE_SIGMA_H_B_0->Draw("AP");
  TGE_SIGMA_H_B_0->GetXaxis()->SetTitle("Laser Tune");
  TGE_SIGMA_H_B_0->GetYaxis()->SetTitle("#delta T_{0} [ns]");
  gPad->BuildLegend();
  gPad->Update();
  gPad->SetGridy();



  
  
  Double_t DELTA_H_B_0_au;
  float DELTA_H_B_0[16];
  float err_DELTA_H_B_0[16];

  
  DELTA_H_B_0_T81->GetPoint(n,x_au,DELTA_H_B_0_au); DELTA_H_B_0[15]=DELTA_H_B_0_au;
  DELTA_H_B_0_T795->GetPoint(n,x_au,DELTA_H_B_0_au); DELTA_H_B_0[14]=DELTA_H_B_0_au;   
  DELTA_H_B_0_T79->GetPoint(n,x_au,DELTA_H_B_0_au); DELTA_H_B_0[13]=DELTA_H_B_0_au;
  DELTA_H_B_0_T77->GetPoint(n,x_au,DELTA_H_B_0_au); DELTA_H_B_0[12]=DELTA_H_B_0_au;
  DELTA_H_B_0_T75->GetPoint(n,x_au,DELTA_H_B_0_au); DELTA_H_B_0[11]=DELTA_H_B_0_au;
  DELTA_H_B_0_T72->GetPoint(n,x_au,DELTA_H_B_0_au); DELTA_H_B_0[10]=DELTA_H_B_0_au;
  DELTA_H_B_0_T73->GetPoint(n,x_au,DELTA_H_B_0_au); DELTA_H_B_0[9]=DELTA_H_B_0_au;
  DELTA_H_B_0_T74->GetPoint(n,x_au,DELTA_H_B_0_au); DELTA_H_B_0[8]=DELTA_H_B_0_au;
  DELTA_H_B_0_T76->GetPoint(n,x_au,DELTA_H_B_0_au); DELTA_H_B_0[7]=DELTA_H_B_0_au;
  DELTA_H_B_0_T78->GetPoint(n,x_au,DELTA_H_B_0_au); DELTA_H_B_0[6]=DELTA_H_B_0_au;
  DELTA_H_B_0_T80->GetPoint(n,x_au,DELTA_H_B_0_au); DELTA_H_B_0[5]=DELTA_H_B_0_au;
  DELTA_H_B_0_T82->GetPoint(n,x_au,DELTA_H_B_0_au); DELTA_H_B_0[4]=DELTA_H_B_0_au;
  DELTA_H_B_0_T84->GetPoint(n,x_au,DELTA_H_B_0_au); DELTA_H_B_0[3]=DELTA_H_B_0_au;
  DELTA_H_B_0_T86->GetPoint(n,x_au,DELTA_H_B_0_au); DELTA_H_B_0[2]=DELTA_H_B_0_au;
  DELTA_H_B_0_T88->GetPoint(n,x_au,DELTA_H_B_0_au); DELTA_H_B_0[1]=DELTA_H_B_0_au;
  DELTA_H_B_0_T90->GetPoint(n,x_au,DELTA_H_B_0_au); DELTA_H_B_0[0]=DELTA_H_B_0_au;


  err_DELTA_H_B_0[15]=DELTA_H_B_0_T81->GetErrorX(n);
  err_DELTA_H_B_0[14]=DELTA_H_B_0_T795->GetErrorX(n);
  err_DELTA_H_B_0[13]=DELTA_H_B_0_T79->GetErrorX(n);
  err_DELTA_H_B_0[12]=DELTA_H_B_0_T77->GetErrorX(n);
  err_DELTA_H_B_0[11]=DELTA_H_B_0_T75->GetErrorX(n);
  err_DELTA_H_B_0[10]=DELTA_H_B_0_T72->GetErrorX(n);
  err_DELTA_H_B_0[9]=DELTA_H_B_0_T73->GetErrorX(n);
  err_DELTA_H_B_0[8]=DELTA_H_B_0_T74->GetErrorX(n);
  err_DELTA_H_B_0[7]=DELTA_H_B_0_T76->GetErrorX(n);
  err_DELTA_H_B_0[6]=DELTA_H_B_0_T78->GetErrorX(n);
  err_DELTA_H_B_0[5]=DELTA_H_B_0_T80->GetErrorX(n);
  err_DELTA_H_B_0[4]=DELTA_H_B_0_T82->GetErrorX(n);
  err_DELTA_H_B_0[3]=DELTA_H_B_0_T84->GetErrorX(n);
  err_DELTA_H_B_0[2]=DELTA_H_B_0_T86->GetErrorX(n);
  err_DELTA_H_B_0[1]=DELTA_H_B_0_T88->GetErrorX(n);
  err_DELTA_H_B_0[0]=DELTA_H_B_0_T90->GetErrorX(n);

  
  TGraphErrors *TGE_DELTA_H_B_0 = new TGraphErrors(16,TUNES,DELTA_H_B_0,err_TUNES,err_DELTA_H_B_0);TGE_DELTA_H_B_0->SetName("DELTA_H_B_0 che 4"); TGE_DELTA_H_B_0->SetTitle(Form("#Delta T ch %d",n)); 
  TGE_DELTA_H_B_0->SetMarkerStyle(20);TGE_DELTA_H_B_0->SetMarkerColor(4);         TGE_DELTA_H_B_0->SetMarkerSize(2);
  TCanvas* c_DELTA_H_B_0_CH = new TCanvas("c_DELTA_H_B_0_CH","c_DELTA_H_B_0_CH");
  TGE_DELTA_H_B_0->Draw("AP");
  TGE_DELTA_H_B_0->GetXaxis()->SetTitle("Laser Tune");
  TGE_DELTA_H_B_0->GetYaxis()->SetTitle("#Delta T [ns]");
  gPad->BuildLegend();
  gPad->Update();
  gPad->SetGridy();

}
