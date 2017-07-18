////////TIME CALIBRATION CODE///////****
////////////////////////////////////
/////To RUN: gROOT->ProcessLine(".L Fit_single_fiber.cpp"); loop_channels(2,"T74",true); > ee.log

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
//#include "RooHist.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooMinimizer.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooMinuit.h"
#include "TAxis.h"
#include <TFile.h>
#include <vector>

using namespace RooFit ;
vector<int> Fits_status;
struct Fit_results{
  float Results[80];
  float Overall_Fracs[4];
  int n_POIs;
  float most_probable_amp[2];
  RooPlot* xframe2_amplitude_0;
  RooPlot* xframe2_fit_0;
  RooPlot* xframe2_pull_0;
  TH2 *h_correlation_0;
  float SP_separation[4];
  
};

struct DT_pars{

  float Par_0[5][8];
  float Par_1[5][8];

};

struct MC_INFO{
  int n_peaks[16];
  float T_separation[16][5];
  float err_T_separation[16][5];
  DT_pars pars;
  
};

void Fit_KEK_data(){

  TFile *file_input = new TFile("run004881_TBC4855-4858_slot01_digits.root");
  TTree *t_input = (TTree*)file_input->Get("laser");
  TH1D *h_temp = new TH1D("h_temp","h_temp",500,-50,0);
  TCut cut = "-20<time&&time<0&&quality==1&&pixel==1";
  t_input->Project("h_temp","time",cut);
  h_temp->Draw();
  Float_t max_bin = h_temp->GetMaximumBin();
  TAxis *xaxis = h_temp->GetXaxis(); 
  Double_t max_pos = xaxis->GetBinCenter(max_bin);
  //cout<< max_pos <<endl;
  TH1D *h_time = new TH1D("h_time","Time [ns]",100,-1,0.5);
  Float_t time_sc;
  Int_t pixelID=-99;
  Int_t quality=-9;
  t_input->SetBranchAddress("time",&time_sc);
  t_input->SetBranchAddress("pixel",&pixelID);
  t_input->SetBranchAddress("quality",&quality);

  Int_t n_entries = t_input->GetEntries();
  for(int i=0; i<n_entries; i++){t_input->GetEntry(i); if(-20<time_sc&&time_sc<0&&quality==1&&pixelID==1) {h_time->Fill(time_sc-max_pos);}} 
 
  h_time->Draw();
}
    

