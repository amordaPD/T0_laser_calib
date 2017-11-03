#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
using namespace std;

#include <assert.h>
#include "TMath.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH1.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "TF1.h"
#include <TFitResultPtr.h>
#include <TFile.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include "TKDE.h"


void make_plots_scan() {


  TFile *_file0 = new TFile("dati/impulsatore-stabilita-lungo-thr--20-T50-F1982-1.root");
    
  TTree *t = (TTree*)_file0->Get("times");
  int n_bin =150;
  float t_min=51.2;
  float t_max=51.75;
  TH1D *h_1 = new TH1D("h_1","h_1",n_bin,t_min,t_max);;
  TH1D *h_2 = new TH1D("h_2","h_2",n_bin,t_min,t_max);;
  TH1D *h_3 = new TH1D("h_3","h_3",n_bin,t_min,t_max);;
  TH1D *h_4 = new TH1D("h_4","h_4",n_bin,t_min,t_max);;
  TH1D *h_5 = new TH1D("h_5","h_5",n_bin,t_min,t_max);;
  TH1D *h_6 = new TH1D("h_6","h_6",n_bin,t_min,t_max);;
  TH1D *h_7 = new TH1D("h_7","h_7",n_bin,t_min,t_max);;
  /*  TH1D *h_8 = new TH1D("h_8","h_8",150,7.2,7.6);;
      TH1D *h_9 = new TH1D("h_all","h_all",150,7.2,7.6);;
  */

  TCut single_photon="20<amps&&amps<40";
  TCut double_photon="50<amps&&amps<68";
  TCut triple_photon="75<amps&&amps<100";
  TCut cut_photon=double_photon;
  float temps[8];
  for(int gg=0;gg<8;gg++){
    temps[gg]=23.0+(4/7.0)*gg;
  }
  TCut temp_cut[7];
  for(int tt=0;tt<7;tt++){
      temp_cut[tt]=Form("%f<=temp&&temp<%f",temps[tt],temps[tt+1]);
  }
  
  t->Project("h_1","times",cut_photon&&temp_cut[0]);
  t->Project("h_2","times",cut_photon&&temp_cut[1]);
  t->Project("h_3","times",cut_photon&&temp_cut[2]);
  t->Project("h_4","times",cut_photon&&temp_cut[3]);
  t->Project("h_5","times",cut_photon&&temp_cut[4]);
  t->Project("h_6","times",cut_photon&&temp_cut[5]);
  t->Project("h_7","times",cut_photon&&temp_cut[6]);
  ;

  
  h_1->SetLineColor(1);
  h_2->SetLineColor(2);
  h_3->SetLineColor(3);
  h_4->SetLineColor(4);
  h_5->SetLineColor(5);
  h_6->SetLineColor(6);
  h_7->SetLineColor(7);/*
			 h_8->SetLineColor(8);
			 h_9->SetLineColor(1);
			 h_9->SetLineWidth(2);*/
  ;

  TCanvas *c_temp = new TCanvas("c_temp","c_temp");
  c_temp->Divide(2,2);
  c_temp->cd(1);
  
  h_1->DrawNormalized("E");
  h_2->DrawNormalized("same:E");
  h_3->DrawNormalized("same:E");
  h_4->DrawNormalized("same:E");
  h_5->DrawNormalized("same:E");
  h_6->DrawNormalized("same:E");
  h_7->DrawNormalized("same:E");/*
  h_8->DrawNormalized("same");
  h_9->DrawNormalized("same");*/
  gPad->BuildLegend();


  c_temp->cd(2);

  float T[7];
  float err_T[7];
  for(int h=0;h<7;h++){
    T[h]=temps[h]+0.5*(temps[h+1]-temps[h]);
    err_T[h]=0;
    cout<<T[h]<<endl;
  }
  
  float means[7];
  float resolution[7];

  means[0]=h_1->GetMean();
  means[1]=h_2->GetMean();
  means[2]=h_3->GetMean();
  means[3]=h_4->GetMean();
  means[4]=h_5->GetMean();
  means[5]=h_6->GetMean();
  means[6]=h_7->GetMean();

  resolution[0]=h_1->GetRMS();
  resolution[1]=h_2->GetRMS();
  resolution[2]=h_3->GetRMS();
  resolution[3]=h_4->GetRMS();
  resolution[4]=h_5->GetRMS();
  resolution[5]=h_6->GetRMS();
  resolution[6]=h_7->GetRMS();

  
  TGraphErrors *FRAC = new TGraphErrors(7,T,means,err_T,resolution); FRAC->SetTitle("time vs temperature");
  TF1 *f = new TF1("f", "[1] * x + [0]");
  FRAC->Fit(f);
  
  FRAC->SetMarkerStyle(20);
  FRAC->SetMarkerSize(2);
  FRAC->SetMarkerColor(2);
  FRAC->GetYaxis()->SetTitle("time [ns]");
  FRAC->GetXaxis()->SetTitle("temperature [C°]");
  FRAC->Draw("AP");

  TH1D *h_t_all_temp = new TH1D("h_t_all_temp","Detection time distribution for 2 photon events",n_bin,t_min,t_max);
  h_t_all_temp->GetXaxis()->SetTitle("time [ns]");
  h_t_all_temp->GetYaxis()->SetTitle("A.U.");
  TH1D *h_t_all_temp_corr = new TH1D("h_t_all_temp_corr","Detection time for 2 photon events",n_bin,t_min,t_max);
  h_t_all_temp_corr->GetXaxis()->SetTitle("time [ns]");
  h_t_all_temp_corr->GetYaxis()->SetTitle("A.U.");
  
  t->Project("h_t_all_temp","times",cut_photon);
  h_t_all_temp->SetLineColor(2);
  float p0 = f->GetParameter(0);
  float p1 = f->GetParameter(1);
  float t0 = 25*p1+p0;
  t->Project("h_t_all_temp_corr",Form("times-((temp-25.0)*%f)",p1),cut_photon);
  h_t_all_temp_corr->SetLineColor(4);
  c_temp->cd(3);
  h_t_all_temp_corr->DrawNormalized();
  h_t_all_temp_corr->Fit("gaus","","",51.3,51.6);
  c_temp->cd(4);
  h_t_all_temp->DrawNormalized();
  h_t_all_temp->Fit("gaus","","",51.2,51.55);
  
 
  /////////////////////////////////
  ////// Amplitude studies ////////
  /////////////////////////////////

  TH1D *h_amp_1p_1 = new TH1D("h_amp_1p_1","h_amp_1p_1",150,20,37);;
  TH1D *h_amp_1p_2 = new TH1D("h_amp_1p_2","h_amp_1p_2",150,20,37);;
  TH1D *h_amp_1p_3 = new TH1D("h_amp_1p_3","h_amp_1p_3",150,20,37);;
  TH1D *h_amp_1p_4 = new TH1D("h_amp_1p_4","h_amp_1p_4",150,20,37);;
  TH1D *h_amp_1p_5 = new TH1D("h_amp_1p_5","h_amp_1p_5",150,20,37);;
  TH1D *h_amp_1p_6 = new TH1D("h_amp_1p_6","h_amp_1p_6",150,20,37);;
  TH1D *h_amp_1p_7 = new TH1D("h_amp_1p_7","h_amp_1p_7",150,20,37);;

  t->Project("h_amp_1p_1","amps",single_photon&&temp_cut[0]);
  t->Project("h_amp_1p_2","amps",single_photon&&temp_cut[1]);
  t->Project("h_amp_1p_3","amps",single_photon&&temp_cut[2]);
  t->Project("h_amp_1p_4","amps",single_photon&&temp_cut[3]);
  t->Project("h_amp_1p_5","amps",single_photon&&temp_cut[4]);
  t->Project("h_amp_1p_6","amps",single_photon&&temp_cut[5]);
  t->Project("h_amp_1p_7","amps",single_photon&&temp_cut[6]);
  
  ;
  TH1D *h_amp_2p_1 = new TH1D("h_amp_2p_1","h_amp_2p_1",150,50,68);;
  TH1D *h_amp_2p_2 = new TH1D("h_amp_2p_2","h_amp_2p_2",150,50,68);;
  TH1D *h_amp_2p_3 = new TH1D("h_amp_2p_3","h_amp_2p_3",150,50,68);;
  TH1D *h_amp_2p_4 = new TH1D("h_amp_2p_4","h_amp_2p_4",150,50,68);;
  TH1D *h_amp_2p_5 = new TH1D("h_amp_2p_5","h_amp_2p_5",150,50,68);;
  TH1D *h_amp_2p_6 = new TH1D("h_amp_2p_6","h_amp_2p_6",150,50,68);;
  TH1D *h_amp_2p_7 = new TH1D("h_amp_2p_7","h_amp_2p_7",150,50,68);;


  t->Project("h_amp_2p_1","amps",double_photon&&temp_cut[0]);
  t->Project("h_amp_2p_2","amps",double_photon&&temp_cut[1]);
  t->Project("h_amp_2p_3","amps",double_photon&&temp_cut[2]);
  t->Project("h_amp_2p_4","amps",double_photon&&temp_cut[3]);
  t->Project("h_amp_2p_5","amps",double_photon&&temp_cut[4]);
  t->Project("h_amp_2p_6","amps",double_photon&&temp_cut[5]);
  t->Project("h_amp_2p_7","amps",double_photon&&temp_cut[6]);

  TH1D *h_amp_3p_1 = new TH1D("h_amp_3p_1","h_amp_3p_1",150,75,100);;
  TH1D *h_amp_3p_2 = new TH1D("h_amp_3p_2","h_amp_3p_2",150,75,100);;
  TH1D *h_amp_3p_3 = new TH1D("h_amp_3p_3","h_amp_3p_3",150,75,100);;
  TH1D *h_amp_3p_4 = new TH1D("h_amp_3p_4","h_amp_3p_4",150,75,100);;
  TH1D *h_amp_3p_5 = new TH1D("h_amp_3p_5","h_amp_3p_5",150,75,100);;
  TH1D *h_amp_3p_6 = new TH1D("h_amp_3p_6","h_amp_3p_6",150,75,100);;
  TH1D *h_amp_3p_7 = new TH1D("h_amp_3p_7","h_amp_3p_7",150,75,100);;


  t->Project("h_amp_3p_1","amps",triple_photon&&temp_cut[0]);
  t->Project("h_amp_3p_2","amps",triple_photon&&temp_cut[1]);
  t->Project("h_amp_3p_3","amps",triple_photon&&temp_cut[2]);
  t->Project("h_amp_3p_4","amps",triple_photon&&temp_cut[3]);
  t->Project("h_amp_3p_5","amps",triple_photon&&temp_cut[4]);
  t->Project("h_amp_3p_6","amps",triple_photon&&temp_cut[5]);
  t->Project("h_amp_3p_7","amps",triple_photon&&temp_cut[6]);

  
  h_amp_1p_1->SetLineColor(1);
  h_amp_1p_2->SetLineColor(2);
  h_amp_1p_3->SetLineColor(3);
  h_amp_1p_4->SetLineColor(4);
  h_amp_1p_5->SetLineColor(5);
  h_amp_1p_6->SetLineColor(6);
  h_amp_1p_7->SetLineColor(7);
  
  h_amp_2p_1->SetLineColor(1);
  h_amp_2p_2->SetLineColor(2);
  h_amp_2p_3->SetLineColor(3);
  h_amp_2p_4->SetLineColor(4);
  h_amp_2p_5->SetLineColor(5);
  h_amp_2p_6->SetLineColor(6);
  h_amp_2p_7->SetLineColor(7);
  
  h_amp_3p_1->SetLineColor(1);
  h_amp_3p_2->SetLineColor(2);
  h_amp_3p_3->SetLineColor(3);
  h_amp_3p_4->SetLineColor(4);
  h_amp_3p_5->SetLineColor(5);
  h_amp_3p_6->SetLineColor(6);
  h_amp_3p_7->SetLineColor(7);
  
  TCanvas *c_amp_dist = new TCanvas("c_amp_dist","c_amp_dist");
  c_amp_dist->Divide(2,2);
  c_amp_dist->cd(1);
  h_amp_1p_1->DrawNormalized("E");
  h_amp_1p_2->DrawNormalized("same:E");
  h_amp_1p_3->DrawNormalized("same:E");
  h_amp_1p_4->DrawNormalized("same:E");
  h_amp_1p_5->DrawNormalized("same:E");
  h_amp_1p_6->DrawNormalized("same:E");
  h_amp_1p_7->DrawNormalized("same:E");
  c_amp_dist->cd(2);
  h_amp_2p_1->DrawNormalized("E");
  h_amp_2p_2->DrawNormalized("same:E");
  h_amp_2p_3->DrawNormalized("same:E");
  h_amp_2p_4->DrawNormalized("same:E");
  h_amp_2p_5->DrawNormalized("same:E");
  h_amp_2p_7->DrawNormalized("same:E");
  gPad->BuildLegend();


  c_amp_dist->cd(3);

  
  float means_amp[7];
  float resolution_amp[7];
  means_amp[0]=h_amp_2p_1->GetMean()-h_amp_1p_1->GetMean();
  means_amp[1]=h_amp_2p_2->GetMean()-h_amp_1p_2->GetMean();
  means_amp[2]=h_amp_2p_3->GetMean()-h_amp_1p_3->GetMean();
  means_amp[3]=h_amp_2p_4->GetMean()-h_amp_1p_4->GetMean();
  means_amp[4]=h_amp_2p_5->GetMean()-h_amp_1p_5->GetMean();
  means_amp[5]=h_amp_2p_6->GetMean()-h_amp_1p_6->GetMean();
  means_amp[6]=h_amp_2p_7->GetMean()-h_amp_1p_7->GetMean();
							  
  resolution_amp[0]=0;//h_amp_2p_1->GetRMS();
  resolution_amp[1]=0;//h_amp_2p_2->GetRMS();
  resolution_amp[2]=0;//h_amp_2p_3->GetRMS();
  resolution_amp[3]=0;//h_amp_2p_4->GetRMS();
  resolution_amp[4]=0;//h_amp_2p_5->GetRMS();
  resolution_amp[5]=0;//h_amp_2p_6->GetRMS();
  resolution_amp[6]=0;//h_amp_2p_6->GetRMS();

  
  TGraphErrors *FRAC_amp = new TGraphErrors(7,T,means_amp,err_T,resolution_amp); FRAC_amp->SetTitle("single-double photon amplitude separation (in [ADC counts]) vs temperature");
  
  FRAC_amp->SetMarkerStyle(20);
  FRAC_amp->SetMarkerSize(2);
  FRAC_amp->SetMarkerColor(2);
  FRAC_amp->GetXaxis()->SetTitle("temperature [C°]");
  FRAC_amp->GetYaxis()->SetTitle("#Delta_{2-1 #gamma}");
  
  
  FRAC_amp->Draw("AP");
  TF1 *f_amp = new TF1("f_amp", "[1] * x + [0]");
  FRAC_amp->Fit(f_amp);


  c_amp_dist->cd(4);
  float means_amp_3[7];
  float resolution_amp[7];
  means_amp_3[0]=h_amp_3p_1->GetMean()-h_amp_2p_1->GetMean();
  means_amp_3[1]=h_amp_3p_2->GetMean()-h_amp_2p_2->GetMean();
  means_amp_3[2]=h_amp_3p_3->GetMean()-h_amp_2p_3->GetMean();
  means_amp_3[3]=h_amp_3p_4->GetMean()-h_amp_2p_4->GetMean();
  means_amp_3[4]=h_amp_3p_5->GetMean()-h_amp_2p_5->GetMean();
  means_amp_3[5]=h_amp_3p_6->GetMean()-h_amp_2p_6->GetMean();
  means_amp_3[6]=h_amp_3p_7->GetMean()-h_amp_2p_7->GetMean();

   TGraphErrors *FRAC_amp_3 = new TGraphErrors(7,T,means_amp_3,err_T,resolution_amp); FRAC_amp_3->SetTitle("single-double photon amplitude separation (in [ADC counts]) vs temperature");
  
  FRAC_amp_3->SetMarkerStyle(20);
  FRAC_amp_3->SetMarkerSize(2);
  FRAC_amp_3->SetMarkerColor(4);
  FRAC_amp_3->GetXaxis()->SetTitle("temperature [C°]");
  FRAC_amp_3->GetYaxis()->SetTitle("#Delta_{3-2 #gamma}");
  
  FRAC_amp_3->Draw("AP:same");
  TF1 *f_amp_3 = new TF1("f_amp_3", "[1] * x + [0]");
  FRAC_amp_3->Fit(f_amp_3);

}
