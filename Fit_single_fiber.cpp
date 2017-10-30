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

//TString input_filename ="runScan-T74-pos0_large_scale_out";

MC_INFO Analyse_MC_inputs(TString _draw_results="draw", TString input_file){
  double a0;
  double a1;
  int col;
  int fn;
  
  TFile* f_input_MC = new TFile(input_file);
  TTree* tree_MC = (TTree*)f_input_MC->Get("inMC");
  tree_MC->SetBranchAddress("a0",&a0);
  tree_MC->SetBranchAddress("a1",&a1);
  tree_MC->SetBranchAddress("col",&col);
  tree_MC->SetBranchAddress("fn",&fn);
  DT_pars parameters_DT;
  for(int n=0;n<7*4;n++){
    tree_MC->GetEntry(n);
    parameters_DT.Par_0[col][fn]=a0;
    parameters_DT.Par_1[col][fn]=a1;
  }

  

  double t_min=0.15;
  double t_max=0.4;
  double t_min_2=-0.1;
  double t_max_2=0.1;
  TH1D *h_D_0  = new TH1D ("h_D_0","h_D_0",250,t_min_2,t_max_2);
  TH1D *h_D_1  = new TH1D ("h_D_1","h_D_1",250,t_min_2,t_max_2);
  TH1D *h_D_2  = new TH1D ("h_D_2","h_D_2",250,t_min_2,t_max_2);
  TH1D *h_D_3  = new TH1D ("h_D_3","h_D_3",250,t_min_2,t_max_2);
  TH1D *h_D_4  = new TH1D ("h_D_4","h_D_4",250,t_min_2,t_max_2);
  TH1D *h_D_5  = new TH1D ("h_D_5","h_D_5",250,t_min_2,t_max_2);
  TH1D *h_D_6  = new TH1D ("h_D_6","h_D_6",250,t_min_2,t_max_2);
  TH1D *h_D_7  = new TH1D ("h_D_7","h_D_7",250,t_min_2,t_max_2);
  TH1D *h_D_8  = new TH1D ("h_D_8","h_D_8",250,t_min_2,t_max_2);
  TH1D *h_D_9  = new TH1D ("h_D_9","h_D_9",250,t_min_2,t_max_2);
  TH1D *h_D_10 = new TH1D ("h_D_10","h_D_10",250,t_min_2,t_max_2);
  TH1D *h_D_11 = new TH1D ("h_D_11","h_D_11",250,t_min_2,t_max_2);
  TH1D *h_D_12 = new TH1D ("h_D_12","h_D_12",250,t_min_2,t_max_2);
  TH1D *h_D_13 = new TH1D ("h_D_13","h_D_13",250,t_min_2,t_max_2);
  TH1D *h_D_14 = new TH1D ("h_D_14","h_D_14",250,t_min_2,t_max_2);
  TH1D *h_D_15 = new TH1D ("h_D_15","h_D_15",250,t_min_2,t_max_2);

  
  TH1D *h_0  = new TH1D ("h_0","h_0",250,t_min,t_max);
  TH1D *h_1  = new TH1D ("h_1","h_1",250,t_min,t_max);
  TH1D *h_2  = new TH1D ("h_2","h_2",250,t_min,t_max);
  TH1D *h_3  = new TH1D ("h_3","h_3",250,t_min,t_max);
  TH1D *h_4  = new TH1D ("h_4","h_4",250,t_min,t_max);
  TH1D *h_5  = new TH1D ("h_5","h_5",250,t_min,t_max);
  TH1D *h_6  = new TH1D ("h_6","h_6",250,t_min,t_max);
  TH1D *h_7  = new TH1D ("h_7","h_7",250,t_min,t_max);
  TH1D *h_8  = new TH1D ("h_8","h_8",250,t_min,t_max);
  TH1D *h_9  = new TH1D ("h_9","h_9",250,t_min,t_max);
  TH1D *h_10 = new TH1D ("h_10","h_10",250,t_min,t_max);
  TH1D *h_11 = new TH1D ("h_11","h_11",250,t_min,t_max);
  TH1D *h_12 = new TH1D ("h_12","h_12",250,t_min,t_max);
  TH1D *h_13 = new TH1D ("h_13","h_13",250,t_min,t_max);
  TH1D *h_14 = new TH1D ("h_14","h_14",250,t_min,t_max);
  TH1D *h_15 = new TH1D ("h_15","h_15",250,t_min,t_max);
  TH1D *h_md_0  = new TH1D ("h_md_0","h_md_0",250,t_min,t_max);
  TH1D *h_md_1  = new TH1D ("h_md_1","h_md_1",250,t_min,t_max);
  TH1D *h_md_2  = new TH1D ("h_md_2","h_md_2",250,t_min,t_max);
  TH1D *h_md_3  = new TH1D ("h_md_3","h_md_3",250,t_min,t_max);
  TH1D *h_md_4  = new TH1D ("h_md_4","h_md_4",250,t_min,t_max);
  TH1D *h_md_5  = new TH1D ("h_md_5","h_md_5",250,t_min,t_max);
  TH1D *h_md_6  = new TH1D ("h_md_6","h_md_6",250,t_min,t_max);
  TH1D *h_md_7  = new TH1D ("h_md_7","h_md_7",250,t_min,t_max);
  TH1D *h_md_8  = new TH1D ("h_md_8","h_md_8",250,t_min,t_max);
  TH1D *h_md_9  = new TH1D ("h_md_9","h_md_9",250,t_min,t_max);
  TH1D *h_md_10 = new TH1D ("h_md_10","h_md_10",250,t_min,t_max);
  TH1D *h_md_11 = new TH1D ("h_md_11","h_md_11",250,t_min,t_max);
  TH1D *h_md_12 = new TH1D ("h_md_12","h_md_12",250,t_min,t_max);
  TH1D *h_md_13 = new TH1D ("h_md_13","h_md_13",250,t_min,t_max);
  TH1D *h_md_14 = new TH1D ("h_md_14","h_md_14",250,t_min,t_max);
  TH1D *h_md_15 = new TH1D ("h_md_15","h_md_15",250,t_min,t_max);
  h_md_0->SetLineColor(2);
  h_md_1->SetLineColor(2);
  h_md_2->SetLineColor(2);
  h_md_3->SetLineColor(2);
  h_md_4->SetLineColor(2);
  h_md_5->SetLineColor(2);
  h_md_6->SetLineColor(2);
  h_md_7->SetLineColor(2);
  h_md_8->SetLineColor(2);
  h_md_9->SetLineColor(2);
  h_md_10->SetLineColor(2);
  h_md_11->SetLineColor(2);
  h_md_12->SetLineColor(2);
  h_md_13->SetLineColor(2);
  h_md_14->SetLineColor(2);
  h_md_15->SetLineColor(2);

  TH1D *h_up_0  = new TH1D ("h_up_0","h_up_0",250,t_min,t_max);
  TH1D *h_up_1  = new TH1D ("h_up_1","h_up_1",250,t_min,t_max);
  TH1D *h_up_2  = new TH1D ("h_up_2","h_up_2",250,t_min,t_max);
  TH1D *h_up_3  = new TH1D ("h_up_3","h_up_3",250,t_min,t_max);
  TH1D *h_up_4  = new TH1D ("h_up_4","h_up_4",250,t_min,t_max);
  TH1D *h_up_5  = new TH1D ("h_up_5","h_up_5",250,t_min,t_max);
  TH1D *h_up_6  = new TH1D ("h_up_6","h_up_6",250,t_min,t_max);
  TH1D *h_up_7  = new TH1D ("h_up_7","h_up_7",250,t_min,t_max);
  TH1D *h_up_8  = new TH1D ("h_up_8","h_up_8",250,t_min,t_max);
  TH1D *h_up_9  = new TH1D ("h_up_9","h_up_9",250,t_min,t_max);
  TH1D *h_up_10 = new TH1D ("h_up_10","h_up_10",250,t_min,t_max);
  TH1D *h_up_11 = new TH1D ("h_up_11","h_up_11",250,t_min,t_max);
  TH1D *h_up_12 = new TH1D ("h_up_12","h_up_12",250,t_min,t_max);
  TH1D *h_up_13 = new TH1D ("h_up_13","h_up_13",250,t_min,t_max);
  TH1D *h_up_14 = new TH1D ("h_up_14","h_up_14",250,t_min,t_max);
  TH1D *h_up_15 = new TH1D ("h_up_15","h_up_15",250,t_min,t_max);
  h_up_0->SetLineColor(1);
  h_up_1->SetLineColor(1);
  h_up_2->SetLineColor(1);
  h_up_3->SetLineColor(1);
  h_up_4->SetLineColor(1);
  h_up_5->SetLineColor(1);
  h_up_6->SetLineColor(1);
  h_up_7->SetLineColor(1);
  h_up_8->SetLineColor(1);
  h_up_9->SetLineColor(1);
  h_up_10->SetLineColor(1);
  h_up_11->SetLineColor(1);
  h_up_12->SetLineColor(1);
  h_up_13->SetLineColor(1);
  h_up_14->SetLineColor(1);
  h_up_15->SetLineColor(1);


  TH1D *h_low_0  = new TH1D ("h_low_0","h_low_0",250,t_min,t_max);
  TH1D *h_low_1  = new TH1D ("h_low_1","h_low_1",250,t_min,t_max);
  TH1D *h_low_2  = new TH1D ("h_low_2","h_low_2",250,t_min,t_max);
  TH1D *h_low_3  = new TH1D ("h_low_3","h_low_3",250,t_min,t_max);
  TH1D *h_low_4  = new TH1D ("h_low_4","h_low_4",250,t_min,t_max);
  TH1D *h_low_5  = new TH1D ("h_low_5","h_low_5",250,t_min,t_max);
  TH1D *h_low_6  = new TH1D ("h_low_6","h_low_6",250,t_min,t_max);
  TH1D *h_low_7  = new TH1D ("h_low_7","h_low_7",250,t_min,t_max);
  TH1D *h_low_8  = new TH1D ("h_low_8","h_low_8",250,t_min,t_max);
  TH1D *h_low_9  = new TH1D ("h_low_9","h_low_9",250,t_min,t_max);
  TH1D *h_low_10 = new TH1D ("h_low_10","h_low_10",250,t_min,t_max);
  TH1D *h_low_11 = new TH1D ("h_low_11","h_low_11",250,t_min,t_max);
  TH1D *h_low_12 = new TH1D ("h_low_12","h_low_12",250,t_min,t_max);
  TH1D *h_low_13 = new TH1D ("h_low_13","h_low_13",250,t_min,t_max);
  TH1D *h_low_14 = new TH1D ("h_low_14","h_low_14",250,t_min,t_max);
  TH1D *h_low_15 = new TH1D ("h_low_15","h_low_15",250,t_min,t_max);
  h_low_0->SetLineColor(8);
  h_low_1->SetLineColor(8);
  h_low_2->SetLineColor(8);
  h_low_3->SetLineColor(8);
  h_low_4->SetLineColor(8);
  h_low_5->SetLineColor(8);
  h_low_6->SetLineColor(8);
  h_low_7->SetLineColor(8);
  h_low_8->SetLineColor(8);
  h_low_9->SetLineColor(8);
  h_low_10->SetLineColor(8);
  h_low_11->SetLineColor(8);
  h_low_12->SetLineColor(8);
  h_low_13->SetLineColor(8);
  h_low_14->SetLineColor(8);
  h_low_15->SetLineColor(8);

  
  
  MC_INFO mc_info;
  TCanvas *c_MC = new TCanvas("c_MC","c_MC");
  c_MC->Divide(4,4);
  TCanvas *c_MC_interval = new TCanvas("c_MC_interval","c_MC_interval");
  c_MC_interval->Divide(4,4);



  vector<float> deltas_low;
  vector<float> deltas_up;
  vector<float> deltas;
  vector<float> mean_delta;
  vector<float> err_mean_delta;
  int n_peaks;
  
  int index_channel_pixel;
  int channelID;
  TFile *f_MC_out = TFile::Open("MC_input_fit.root","RECREATE");
  TTree *t = new TTree("MC_inputs","MC_inputs");
  t->Branch("mean_delta",&mean_delta);
  t->Branch("err_mean_delta",&err_mean_delta);
  t->Branch("ch",&index_channel_pixel);
  t->Branch("peaks",&n_peaks);
  
  for(int column=1; column<=4;column++){
    for(int row=1; row<=4;row++){
      
      deltas_low.clear();
      deltas_up.clear();
      deltas.clear();
      deltas.push_back(0);
      deltas_low.push_back(0);
      deltas_up.push_back(0);
      for(int nF=1;nF<=7;nF++){
	deltas.push_back((parameters_DT.Par_0[column][nF]+(row-0.5)*parameters_DT.Par_1[column][nF])-(parameters_DT.Par_0[column][1]+(row-0.5)*parameters_DT.Par_1[column][1]));
	deltas_low.push_back((parameters_DT.Par_0[column][nF]+(row)*parameters_DT.Par_1[column][nF])-(parameters_DT.Par_0[column][1]+(row)*parameters_DT.Par_1[column][1]));
	deltas_up.push_back((parameters_DT.Par_0[column][nF]+(row-1)*parameters_DT.Par_1[column][nF])-(parameters_DT.Par_0[column][1]+(row-1)*parameters_DT.Par_1[column][1]));
      }
      //for(int ll=0; ll<deltas.size();ll++){cout<<deltas[ll]<<endl;}
      std::sort (deltas.begin(), deltas.end()); 
      //for(int ll=0; ll<deltas.size();ll++){cout<<deltas[ll]<<endl;}

      //////////////////////////////////////////////////////////
      ///////CLUSTERING OF PEAKS - BEGIN ///////////////////////
      //////////////////////////////////////////////////////////
      
      mean_delta.clear();
      err_mean_delta.clear();
      int low_nF, up_nF;
      up_nF=2;
      low_nF=2;
      float thrs_rel = 0.025;
      float thrs_range = 0.06;
      for(int nF=3;nF<=7;nF++){
	if(nF<7){
	  if(TMath::Abs(deltas[nF]-deltas[nF-1])<thrs_rel&&TMath::Abs(deltas[nF]-deltas[low_nF])<thrs_range){up_nF=nF;}
	  else{
	    mean_delta.push_back(deltas[low_nF]+0.5*(deltas[up_nF]-deltas[low_nF]));
	    err_mean_delta.push_back(0.5*(deltas[up_nF]-deltas[low_nF]));
	    low_nF=nF;}
	}else{
	  up_nF=nF;
	  mean_delta.push_back(deltas[low_nF]+0.5*(deltas[up_nF]-deltas[low_nF]));
	  err_mean_delta.push_back(0.5*(deltas[up_nF]-deltas[low_nF]));
	}
      }
      //for(int ll=0; ll<mean_delta.size();ll++){cout<<mean_delta[ll]<<endl;}
      //cout<<"%%%%%%%%%%%%%%%%%%"<<endl;
      //////////////////////////////////////////////////////////
      //////////   CLUSTERING OF PEAKS - END ///////////////////
      //////////////////////////////////////////////////////////
      
      if(column==1&&row==1){channelID=0;}//1+0; }//ch=0
      if(column==1&&row==2){channelID=1;}//1+4; }//ch=0
      if(column==1&&row==3){channelID=2;}//1+8; }//ch=0
      if(column==1&&row==4){channelID=3;}//1+12;}//ch=0
      if(column==2&&row==1){channelID=4;}//1+1;}//ch=0
      if(column==2&&row==2){channelID=5;}//1+5;}//ch=0
      if(column==2&&row==3){channelID=6;}//1+9;}//ch=0
      if(column==2&&row==4){channelID=7;}//1+13;}//ch=0
      if(column==3&&row==1){channelID=8;}//1+2;}//ch=0
      if(column==3&&row==2){channelID=9;}//1+6;}//ch=0
      if(column==3&&row==3){channelID=10;}//1+10;}//ch=0
      if(column==3&&row==4){channelID=11;}//1+14;}//ch=0
      if(column==4&&row==1){channelID=12;}//1+3;}//ch=0
      if(column==4&&row==2){channelID=13;}//1+7;}//ch=0
      if(column==4&&row==3){channelID=14;}//1+11;}//ch=0
      if(column==4&&row==4){channelID=15;}//1+15;}//ch=0
      
      if(column==1&&row==1){index_channel_pixel=1+0;  for(int ll=1; ll<=7;ll++){h_0->Fill(deltas[ll]);}   c_MC->cd(index_channel_pixel);h_0->Draw();}//ch=0
      if(column==1&&row==2){index_channel_pixel=1+4;  for(int ll=1; ll<=7;ll++){h_1->Fill(deltas[ll]);}   c_MC->cd(index_channel_pixel);h_1->Draw();}//ch=0
      if(column==1&&row==3){index_channel_pixel=1+8;  for(int ll=1; ll<=7;ll++){h_2->Fill(deltas[ll]);}   c_MC->cd(index_channel_pixel);h_2->Draw();}//ch=0
      if(column==1&&row==4){index_channel_pixel=1+12; for(int ll=1; ll<=7;ll++){h_3->Fill(deltas[ll]);}   c_MC->cd(index_channel_pixel);h_3->Draw();}//ch=0
      if(column==2&&row==1){index_channel_pixel=1+1;  for(int ll=1; ll<=7;ll++){h_4->Fill(deltas[ll]);}   c_MC->cd(index_channel_pixel);h_4->Draw();}//ch=0
      if(column==2&&row==2){index_channel_pixel=1+5;  for(int ll=1; ll<=7;ll++){h_5->Fill(deltas[ll]);}   c_MC->cd(index_channel_pixel);h_5->Draw();}//ch=0
      if(column==2&&row==3){index_channel_pixel=1+9;  for(int ll=1; ll<=7;ll++){h_6->Fill(deltas[ll]);}   c_MC->cd(index_channel_pixel);h_6->Draw();}//ch=0
      if(column==2&&row==4){index_channel_pixel=1+13; for(int ll=1; ll<=7;ll++){h_7->Fill(deltas[ll]);}   c_MC->cd(index_channel_pixel);h_7->Draw();}//ch=0
      if(column==3&&row==1){index_channel_pixel=1+2;  for(int ll=1; ll<=7;ll++){h_8->Fill(deltas[ll]);}   c_MC->cd(index_channel_pixel);h_8->Draw();}//ch=0
      if(column==3&&row==2){index_channel_pixel=1+6;  for(int ll=1; ll<=7;ll++){h_9->Fill(deltas[ll]);}   c_MC->cd(index_channel_pixel);h_9->Draw();}//ch=0
      if(column==3&&row==3){index_channel_pixel=1+10; for(int ll=1; ll<=7;ll++){h_10->Fill(deltas[ll]);}  c_MC->cd(index_channel_pixel);h_10->Draw();}//ch=0
      if(column==3&&row==4){index_channel_pixel=1+14; for(int ll=1; ll<=7;ll++){h_11->Fill(deltas[ll]);}  c_MC->cd(index_channel_pixel);h_11->Draw();}//ch=0
      if(column==4&&row==1){index_channel_pixel=1+3;  for(int ll=1; ll<=7;ll++){h_12->Fill(deltas[ll]);}  c_MC->cd(index_channel_pixel);h_12->Draw();}//ch=0
      if(column==4&&row==2){index_channel_pixel=1+7;  for(int ll=1; ll<=7;ll++){h_13->Fill(deltas[ll]);}  c_MC->cd(index_channel_pixel);h_13->Draw();}//ch=0
      if(column==4&&row==3){index_channel_pixel=1+11; for(int ll=1; ll<=7;ll++){h_14->Fill(deltas[ll]);}  c_MC->cd(index_channel_pixel);h_14->Draw();}//ch=0
      if(column==4&&row==4){index_channel_pixel=1+15; for(int ll=1; ll<=7;ll++){h_15->Fill(deltas[ll]);}  c_MC->cd(index_channel_pixel);h_15->Draw();}//ch=0


      if(column==1&&row==1){index_channel_pixel=1+0;  for(int ll=0; ll<mean_delta.size();ll++){h_md_0->Fill(mean_delta[ll]);}   c_MC->cd(index_channel_pixel);h_md_0->Draw("same");}//ch=0
      if(column==1&&row==2){index_channel_pixel=1+4;  for(int ll=0; ll<mean_delta.size();ll++){h_md_1->Fill(mean_delta[ll]);}   c_MC->cd(index_channel_pixel);h_md_1->Draw("same");}//ch=0
      if(column==1&&row==3){index_channel_pixel=1+8;  for(int ll=0; ll<mean_delta.size();ll++){h_md_2->Fill(mean_delta[ll]);}   c_MC->cd(index_channel_pixel);h_md_2->Draw("same");}//ch=0
      if(column==1&&row==4){index_channel_pixel=1+12; for(int ll=0; ll<mean_delta.size();ll++){h_md_3->Fill(mean_delta[ll]);}   c_MC->cd(index_channel_pixel);h_md_3->Draw("same");}//ch=0
      if(column==2&&row==1){index_channel_pixel=1+1;  for(int ll=0; ll<mean_delta.size();ll++){h_md_4->Fill(mean_delta[ll]);}   c_MC->cd(index_channel_pixel);h_md_4->Draw("same");}//ch=0
      if(column==2&&row==2){index_channel_pixel=1+5;  for(int ll=0; ll<mean_delta.size();ll++){h_md_5->Fill(mean_delta[ll]);}   c_MC->cd(index_channel_pixel);h_md_5->Draw("same");}//ch=0
      if(column==2&&row==3){index_channel_pixel=1+9;  for(int ll=0; ll<mean_delta.size();ll++){h_md_6->Fill(mean_delta[ll]);}   c_MC->cd(index_channel_pixel);h_md_6->Draw("same");}//ch=0
      if(column==2&&row==4){index_channel_pixel=1+13; for(int ll=0; ll<mean_delta.size();ll++){h_md_7->Fill(mean_delta[ll]);}   c_MC->cd(index_channel_pixel);h_md_7->Draw("same");}//ch=0
      if(column==3&&row==1){index_channel_pixel=1+2;  for(int ll=0; ll<mean_delta.size();ll++){h_md_8->Fill(mean_delta[ll]);}   c_MC->cd(index_channel_pixel);h_md_8->Draw("same");}//ch=0
      if(column==3&&row==2){index_channel_pixel=1+6;  for(int ll=0; ll<mean_delta.size();ll++){h_md_9->Fill(mean_delta[ll]);}   c_MC->cd(index_channel_pixel);h_md_9->Draw("same");}//ch=0
      if(column==3&&row==3){index_channel_pixel=1+10; for(int ll=0; ll<mean_delta.size();ll++){h_md_10->Fill(mean_delta[ll]);}  c_MC->cd(index_channel_pixel);h_md_10->Draw("same");}//ch=0
      if(column==3&&row==4){index_channel_pixel=1+14; for(int ll=0; ll<mean_delta.size();ll++){h_md_11->Fill(mean_delta[ll]);}  c_MC->cd(index_channel_pixel);h_md_11->Draw("same");}//ch=0
      if(column==4&&row==1){index_channel_pixel=1+3;  for(int ll=0; ll<mean_delta.size();ll++){h_md_12->Fill(mean_delta[ll]);}  c_MC->cd(index_channel_pixel);h_md_12->Draw("same");}//ch=0
      if(column==4&&row==2){index_channel_pixel=1+7;  for(int ll=0; ll<mean_delta.size();ll++){h_md_13->Fill(mean_delta[ll]);}  c_MC->cd(index_channel_pixel);h_md_13->Draw("same");}//ch=0
      if(column==4&&row==3){index_channel_pixel=1+11; for(int ll=0; ll<mean_delta.size();ll++){h_md_14->Fill(mean_delta[ll]);}  c_MC->cd(index_channel_pixel);h_md_14->Draw("same");}//ch=0
      if(column==4&&row==4){index_channel_pixel=1+15; for(int ll=0; ll<mean_delta.size();ll++){h_md_15->Fill(mean_delta[ll]);}  c_MC->cd(index_channel_pixel);h_md_15->Draw("same");}//ch=0

      
      if(column==1&&row==1){index_channel_pixel=1+0;  for(int ll=0; ll<deltas_up.size();ll++){h_up_0->Fill(deltas_up[ll]);}   c_MC->cd(index_channel_pixel);h_up_0->Draw("same");}//ch=0
      if(column==1&&row==2){index_channel_pixel=1+4;  for(int ll=0; ll<deltas_up.size();ll++){h_up_1->Fill(deltas_up[ll]);}   c_MC->cd(index_channel_pixel);h_up_1->Draw("same");}//ch=0
      if(column==1&&row==3){index_channel_pixel=1+8;  for(int ll=0; ll<deltas_up.size();ll++){h_up_2->Fill(deltas_up[ll]);}   c_MC->cd(index_channel_pixel);h_up_2->Draw("same");}//ch=0
      if(column==1&&row==4){index_channel_pixel=1+12; for(int ll=0; ll<deltas_up.size();ll++){h_up_3->Fill(deltas_up[ll]);}   c_MC->cd(index_channel_pixel);h_up_3->Draw("same");}//ch=0
      if(column==2&&row==1){index_channel_pixel=1+1;  for(int ll=0; ll<deltas_up.size();ll++){h_up_4->Fill(deltas_up[ll]);}   c_MC->cd(index_channel_pixel);h_up_4->Draw("same");}//ch=0
      if(column==2&&row==2){index_channel_pixel=1+5;  for(int ll=0; ll<deltas_up.size();ll++){h_up_5->Fill(deltas_up[ll]);}   c_MC->cd(index_channel_pixel);h_up_5->Draw("same");}//ch=0
      if(column==2&&row==3){index_channel_pixel=1+9;  for(int ll=0; ll<deltas_up.size();ll++){h_up_6->Fill(deltas_up[ll]);}   c_MC->cd(index_channel_pixel);h_up_6->Draw("same");}//ch=0
      if(column==2&&row==4){index_channel_pixel=1+13; for(int ll=0; ll<deltas_up.size();ll++){h_up_7->Fill(deltas_up[ll]);}   c_MC->cd(index_channel_pixel);h_up_7->Draw("same");}//ch=0
      if(column==3&&row==1){index_channel_pixel=1+2;  for(int ll=0; ll<deltas_up.size();ll++){h_up_8->Fill(deltas_up[ll]);}   c_MC->cd(index_channel_pixel);h_up_8->Draw("same");}//ch=0
      if(column==3&&row==2){index_channel_pixel=1+6;  for(int ll=0; ll<deltas_up.size();ll++){h_up_9->Fill(deltas_up[ll]);}   c_MC->cd(index_channel_pixel);h_up_9->Draw("same");}//ch=0
      if(column==3&&row==3){index_channel_pixel=1+10; for(int ll=0; ll<deltas_up.size();ll++){h_up_10->Fill(deltas_up[ll]);}  c_MC->cd(index_channel_pixel);h_up_10->Draw("same");}//ch=0
      if(column==3&&row==4){index_channel_pixel=1+14; for(int ll=0; ll<deltas_up.size();ll++){h_up_11->Fill(deltas_up[ll]);}  c_MC->cd(index_channel_pixel);h_up_11->Draw("same");}//ch=0
      if(column==4&&row==1){index_channel_pixel=1+3;  for(int ll=0; ll<deltas_up.size();ll++){h_up_12->Fill(deltas_up[ll]);}  c_MC->cd(index_channel_pixel);h_up_12->Draw("same");}//ch=0
      if(column==4&&row==2){index_channel_pixel=1+7;  for(int ll=0; ll<deltas_up.size();ll++){h_up_13->Fill(deltas_up[ll]);}  c_MC->cd(index_channel_pixel);h_up_13->Draw("same");}//ch=0
      if(column==4&&row==3){index_channel_pixel=1+11; for(int ll=0; ll<deltas_up.size();ll++){h_up_14->Fill(deltas_up[ll]);}  c_MC->cd(index_channel_pixel);h_up_14->Draw("same");}//ch=0
      if(column==4&&row==4){index_channel_pixel=1+15; for(int ll=0; ll<deltas_up.size();ll++){h_up_15->Fill(deltas_up[ll]);}  c_MC->cd(index_channel_pixel);h_up_15->Draw("same");}//ch=0
     
      if(column==1&&row==1){index_channel_pixel=1+0;  for(int ll=0; ll<deltas_low.size();ll++){h_low_0->Fill(deltas_low[ll]);}   c_MC->cd(index_channel_pixel);h_low_0->Draw("same");}//ch=0
      if(column==1&&row==2){index_channel_pixel=1+4;  for(int ll=0; ll<deltas_low.size();ll++){h_low_1->Fill(deltas_low[ll]);}   c_MC->cd(index_channel_pixel);h_low_1->Draw("same");}//ch=0
      if(column==1&&row==3){index_channel_pixel=1+8;  for(int ll=0; ll<deltas_low.size();ll++){h_low_2->Fill(deltas_low[ll]);}   c_MC->cd(index_channel_pixel);h_low_2->Draw("same");}//ch=0
      if(column==1&&row==4){index_channel_pixel=1+12; for(int ll=0; ll<deltas_low.size();ll++){h_low_3->Fill(deltas_low[ll]);}   c_MC->cd(index_channel_pixel);h_low_3->Draw("same");}//ch=0
      if(column==2&&row==1){index_channel_pixel=1+1;  for(int ll=0; ll<deltas_low.size();ll++){h_low_4->Fill(deltas_low[ll]);}   c_MC->cd(index_channel_pixel);h_low_4->Draw("same");}//ch=0
      if(column==2&&row==2){index_channel_pixel=1+5;  for(int ll=0; ll<deltas_low.size();ll++){h_low_5->Fill(deltas_low[ll]);}   c_MC->cd(index_channel_pixel);h_low_5->Draw("same");}//ch=0
      if(column==2&&row==3){index_channel_pixel=1+9;  for(int ll=0; ll<deltas_low.size();ll++){h_low_6->Fill(deltas_low[ll]);}   c_MC->cd(index_channel_pixel);h_low_6->Draw("same");}//ch=0
      if(column==2&&row==4){index_channel_pixel=1+13; for(int ll=0; ll<deltas_low.size();ll++){h_low_7->Fill(deltas_low[ll]);}   c_MC->cd(index_channel_pixel);h_low_7->Draw("same");}//ch=0
      if(column==3&&row==1){index_channel_pixel=1+2;  for(int ll=0; ll<deltas_low.size();ll++){h_low_8->Fill(deltas_low[ll]);}   c_MC->cd(index_channel_pixel);h_low_8->Draw("same");}//ch=0
      if(column==3&&row==2){index_channel_pixel=1+6;  for(int ll=0; ll<deltas_low.size();ll++){h_low_9->Fill(deltas_low[ll]);}   c_MC->cd(index_channel_pixel);h_low_9->Draw("same");}//ch=0
      if(column==3&&row==3){index_channel_pixel=1+10; for(int ll=0; ll<deltas_low.size();ll++){h_low_10->Fill(deltas_low[ll]);}  c_MC->cd(index_channel_pixel);h_low_10->Draw("same");}//ch=0
      if(column==3&&row==4){index_channel_pixel=1+14; for(int ll=0; ll<deltas_low.size();ll++){h_low_11->Fill(deltas_low[ll]);}  c_MC->cd(index_channel_pixel);h_low_11->Draw("same");}//ch=0
      if(column==4&&row==1){index_channel_pixel=1+3;  for(int ll=0; ll<deltas_low.size();ll++){h_low_12->Fill(deltas_low[ll]);}  c_MC->cd(index_channel_pixel);h_low_12->Draw("same");}//ch=0
      if(column==4&&row==2){index_channel_pixel=1+7;  for(int ll=0; ll<deltas_low.size();ll++){h_low_13->Fill(deltas_low[ll]);}  c_MC->cd(index_channel_pixel);h_low_13->Draw("same");}//ch=0
      if(column==4&&row==3){index_channel_pixel=1+11; for(int ll=0; ll<deltas_low.size();ll++){h_low_14->Fill(deltas_low[ll]);}  c_MC->cd(index_channel_pixel);h_low_14->Draw("same");}//ch=0
      if(column==4&&row==4){index_channel_pixel=1+15; for(int ll=0; ll<deltas_low.size();ll++){h_low_15->Fill(deltas_low[ll]);}  c_MC->cd(index_channel_pixel);h_low_15->Draw("same");}//ch=0


      
      if(column==1&&row==1){index_channel_pixel=1+0;  for(int ll=0; ll<deltas_low.size();ll++){h_D_0->Fill(deltas_up[ll]-deltas_low[ll]);}   c_MC_interval->cd(index_channel_pixel);h_D_0->Draw("same");}//ch=0
      if(column==1&&row==2){index_channel_pixel=1+4;  for(int ll=0; ll<deltas_low.size();ll++){h_D_1->Fill(deltas_up[ll]-deltas_low[ll]);}   c_MC_interval->cd(index_channel_pixel);h_D_1->Draw("same");}//ch=0
      if(column==1&&row==3){index_channel_pixel=1+8;  for(int ll=0; ll<deltas_low.size();ll++){h_D_2->Fill(deltas_up[ll]-deltas_low[ll]);}   c_MC_interval->cd(index_channel_pixel);h_D_2->Draw("same");}//ch=0
      if(column==1&&row==4){index_channel_pixel=1+12; for(int ll=0; ll<deltas_low.size();ll++){h_D_3->Fill(deltas_up[ll]-deltas_low[ll]);}   c_MC_interval->cd(index_channel_pixel);h_D_3->Draw("same");}//ch=0
      if(column==2&&row==1){index_channel_pixel=1+1;  for(int ll=0; ll<deltas_low.size();ll++){h_D_4->Fill(deltas_up[ll]-deltas_low[ll]);}   c_MC_interval->cd(index_channel_pixel);h_D_4->Draw("same");}//ch=0
      if(column==2&&row==2){index_channel_pixel=1+5;  for(int ll=0; ll<deltas_low.size();ll++){h_D_5->Fill(deltas_up[ll]-deltas_low[ll]);}   c_MC_interval->cd(index_channel_pixel);h_D_5->Draw("same");}//ch=0
      if(column==2&&row==3){index_channel_pixel=1+9;  for(int ll=0; ll<deltas_low.size();ll++){h_D_6->Fill(deltas_up[ll]-deltas_low[ll]);}   c_MC_interval->cd(index_channel_pixel);h_D_6->Draw("same");}//ch=0
      if(column==2&&row==4){index_channel_pixel=1+13; for(int ll=0; ll<deltas_low.size();ll++){h_D_7->Fill(deltas_up[ll]-deltas_low[ll]);}   c_MC_interval->cd(index_channel_pixel);h_D_7->Draw("same");}//ch=0
      if(column==3&&row==1){index_channel_pixel=1+2;  for(int ll=0; ll<deltas_low.size();ll++){h_D_8->Fill(deltas_up[ll]-deltas_low[ll]);}   c_MC_interval->cd(index_channel_pixel);h_D_8->Draw("same");}//ch=0
      if(column==3&&row==2){index_channel_pixel=1+6;  for(int ll=0; ll<deltas_low.size();ll++){h_D_9->Fill(deltas_up[ll]-deltas_low[ll]);}   c_MC_interval->cd(index_channel_pixel);h_D_9->Draw("same");}//ch=0
      if(column==3&&row==3){index_channel_pixel=1+10; for(int ll=0; ll<deltas_low.size();ll++){h_D_10->Fill(deltas_up[ll]-deltas_low[ll]);}   c_MC_interval->cd(index_channel_pixel);h_D_10->Draw("same");}//ch=0
      if(column==3&&row==4){index_channel_pixel=1+14; for(int ll=0; ll<deltas_low.size();ll++){h_D_11->Fill(deltas_up[ll]-deltas_low[ll]);}   c_MC_interval->cd(index_channel_pixel);h_D_11->Draw("same");}//ch=0
      if(column==4&&row==1){index_channel_pixel=1+3;  for(int ll=0; ll<deltas_low.size();ll++){h_D_12->Fill(deltas_up[ll]-deltas_low[ll]);}   c_MC_interval->cd(index_channel_pixel);h_D_12->Draw("same");}//ch=0
      if(column==4&&row==2){index_channel_pixel=1+7;  for(int ll=0; ll<deltas_low.size();ll++){h_D_13->Fill(deltas_up[ll]-deltas_low[ll]);}   c_MC_interval->cd(index_channel_pixel);h_D_13->Draw("same");}//ch=0
      if(column==4&&row==3){index_channel_pixel=1+11; for(int ll=0; ll<deltas_low.size();ll++){h_D_14->Fill(deltas_up[ll]-deltas_low[ll]);}   c_MC_interval->cd(index_channel_pixel);h_D_14->Draw("same");}//ch=0
      if(column==4&&row==4){index_channel_pixel=1+15; for(int ll=0; ll<deltas_low.size();ll++){h_D_15->Fill(deltas_up[ll]-deltas_low[ll]);}   c_MC_interval->cd(index_channel_pixel);h_D_15->Draw("same");}//ch=0
      n_peaks=mean_delta.size()+1;
      mc_info.n_peaks[channelID]=n_peaks;
      for(int hh=0;hh<mean_delta.size();hh++){
	mc_info.T_separation[channelID][hh]=mean_delta[hh];
	mc_info.err_T_separation[channelID][hh]=err_mean_delta[hh];
      }
      t->Fill();
      
    }
  }			


  f_MC_out->Write();
  delete f_MC_out;

  return mc_info;
}
  
Fit_results Fit_head(string _draw_results="draw", int fix_params=2, TString input_tune, int ch =0 ){


  MC_INFO mc_info =  Analyse_MC_inputs("draw","inFit.root");



  
  TFile *f_input_histogram = new TFile("flat_ntuples/runScan-"+input_tune+"-pos0_large_scale_out.root"); 
  //TFile *f_input_histogram = new TFile("flat_ntuples/runScan-"+input_tune+"-pos0_out.root"); 
  //TFile *f_input_histogram = new TFile("flat_ntuples/runScan-"+input_tune+"-pos0_new_RR_out.root"); 
  bool do_prefit=true;
  bool do_prefit_fullSpectrum = true;
  bool use_NLL=true; //to set the use of fitTo method of RooAbsPdf or the explicit construction of the nll
  int MN_output_print_level=0;
  int MN_output_print_level_prefit;
  bool print_prefit_info=false;
  const char * Type_minim="Minuit";
  const char * Algo_minim="minimize";//
  const char * Type_minim_pf="Minuit";//"Minuit2";//
  const char * Algo_minim_pf="minimize";//"scan";//
  bool do_simultaneous_fit=false;
  bool add_third_signal=false;
  bool simulate_CB_tail=false;
  bool binned_fit=true;//if false fit is unbinned
  int amplitude_cut = -40;
  bool fit_highest_peak=true;
  bool no_grease=true;
  bool fix_deltas=false;
  bool add_SP_components=false;
  bool use_marginalize_method=false;
  bool gaussian_constraints_delta = true;
  //fix_deltas=true;
  //add_third_signal=true;//if(!(ch==0||ch==4||ch==8||ch==12)){  add_third_signal=true;}
  //fix_deltas=true;
  
  if(!print_prefit_info){MN_output_print_level_prefit=-1;}else{MN_output_print_level_prefit=MN_output_print_level;}
  bool draw_results;
  if(_draw_results=="draw"){draw_results=true;}else if(_draw_results=="blind"){draw_results=false;}else{draw_results=false;}
  Fit_results my_fit_results;
  TString channel_0 = Form("fiber0-0-%d",ch);
  TString channel_A_0 = Form("fibertA0-0-%d",ch);
  //bool fix_params = true;
  vector<float> POIs;
  POIs.clear();
  // S e t u p   m o d e l 
  // ---------------------
  
  double my_low_x=22.2;//
  double my_up_x=24;//24;
  if(add_SP_components){my_up_x=my_low_x+13;}
  RooRealVar x("Time","Time [ns]",my_low_x,my_up_x) ;
  RooRealVar amp("Amplitude","Amplitude [ADC counts]",-200,0) ;
  RooRealVar CH("Channel","PMT Channel",0,16) ;
  ///////fixing starting values and boundaries for two positions


  
  TTree *tree = (TTree*)f_input_histogram->Get("tree_input");
  TH2D *h_amp_histogram_0 = (TH2D*)f_input_histogram->Get(channel_A_0);
  RooDataSet ds_0_amp("ds_0_amp","ds_0_amp", RooArgSet(x,amp,CH),Import(*tree),Cut(Form("Amplitude<%d &&Channel==%d",amplitude_cut,ch)));
  RooPlot* xframe2_0_amp = amp.frame(Title(Form("amplitude pos. 0, channel %d",ch))) ;
  ds_0_amp.plotOn(xframe2_0_amp);//,DataError(RooAbsData::SumW2)) ;
  
  if(binned_fit){
    TH1 *h_input_histogram_0 = (TH1*)f_input_histogram->Get(channel_0);
    RooDataHist ds_0("ds_0","ds_0",RooArgSet(x),Import(*h_input_histogram_0)) ;
    cout<<Form("dataset 0 ch %d info :",ch)<<ds_0.Print("v")<<endl;
  }else{
    TTree *tree = (TTree*)f_input_histogram->Get("tree_input");
    RooDataSet ds_0("ds_0","ds_0", RooArgSet(x,amp,CH),Import(*tree),Cut(Form("Amplitude<%d &&Channel==%d",amplitude_cut,ch)));
  }


  double low_x_0;
  double up_x_0;
  double starting_mean_L_0;
  double low_mean_L_0;
  double up_mean_L_0;
  double starting_mean_H_0;
  double low_mean_H_0;
  double up_mean_H_0;
  double low_delta_H_0;
  double up_delta_H_0;
  double starting_delta_H_0;
  double low_delta_T_0;
  double up_delta_T_0;
  double starting_delta_T_0;
  double starting_sigma_T_0;
  double low_sigma_T_0;
  double up_sigma_T_0;
  double starting_sigma_L_0;
  double low_sigma_L_0;
  double up_sigma_L_0;
  double starting_sigma_H_0;
  double low_sigma_H_0;
  double up_sigma_H_0; 
  double starting_Frac_H_0;
  double      low_Frac_L_0;
  double       up_Frac_L_0;
  double starting_Frac_L_0;
  double      low_Frac_H_0;
  double       up_Frac_H_0;
  double starting_Frac_T_0;
  double      low_Frac_T_0;
  double       up_Frac_T_0;
  double starting_alpha_0;
  double      low_alpha_0;
  double       up_alpha_0;
  double starting_beta_0;
  double      low_beta_0;
  double       up_beta_0;




  
  double starting_alpha_CB;
  double low_alpha_CB;
  double up_alpha_CB;
  double starting_n_CB;
  double low_n_CB;
  double up_n_CB;
  
  
  /////POS 0
   low_x_0=my_low_x; up_x_0=my_up_x;
  
  low_delta_H_0=0.01080;
   up_delta_H_0=0.5;
  
  low_delta_T_0=-0.0;
   up_delta_T_0=0.5;
  
  low_sigma_L_0=0.035;
   up_sigma_L_0=0.500;
  
  low_sigma_H_0=0.045;
   up_sigma_H_0=0.5180;
  
  low_sigma_T_0=0.050;
   up_sigma_T_0=0.200;
  
  low_alpha_0=0.01;
   up_alpha_0=0.8;
  
  low_beta_0=0.2;
   up_beta_0=1.0;

   if(no_grease){
     starting_mean_L_0=22.8;
     if(ch==14)  starting_mean_L_0=22.6;
     if(ch==12)  starting_mean_L_0=22.8;
     /* if(ch==5)  starting_mean_L_0=22.6;
     if(ch==9||ch==1)  starting_mean_L_0=22.8;
     if(ch==12||ch==14) starting_mean_L_0=22.8;
     if(ch==8) starting_mean_L_0=23.0;*/
     //if(ch==5||ch==8) starting_mean_L_0=23.0;
   }else{
     starting_mean_L_0=22.7;
     if(ch==14) starting_mean_L_0=22.5;

   }
   if(mc_info.n_peaks[ch]>2){add_third_signal=false;}
   //starting_delta_H_0=0.280;
   //starting_delta_T_0=0.1200;
   starting_delta_H_0=mc_info.T_separation[ch][0];
   starting_delta_T_0=mc_info.T_separation[ch][1]-mc_info.T_separation[ch][0];if(ch==12){cout<<" 'nnamo : "<<mc_info.T_separation[ch][1]-mc_info.T_separation[ch][0]<<endl;}
   starting_sigma_L_0=0.0810;
   starting_sigma_H_0=0.0810;
   starting_sigma_T_0=0.1000;
   starting_alpha_0=0.5;
   starting_beta_0=0.5;
   cout<<"ridaje     "<<starting_delta_H_0<<"    "<<starting_delta_T_0<<endl;
   //  low_mean_L_0=starting_mean_L_0-0.15;
   //   up_mean_L_0=starting_mean_L_0+0.15;
  low_mean_L_0=starting_mean_L_0-0.35;
   up_mean_L_0=starting_mean_L_0+0.35;
  starting_mean_H_0=starting_mean_L_0+starting_delta_H_0;
   low_mean_H_0=low_mean_L_0+starting_delta_H_0;
   up_mean_H_0=up_mean_L_0+starting_delta_H_0;
   
   /*
   if(fix_deltas){
     if(ch==0) starting_delta_H_0=0.280;
     if(ch==1) starting_delta_H_0=0.261;
     if(ch==2) starting_delta_H_0=0.221;
     if(ch==3) starting_delta_H_0=0.176;
     if(ch==4) starting_delta_H_0=0.280;
     if(ch==5) starting_delta_H_0=0.262;
     if(ch==6) starting_delta_H_0=0.222;
     if(ch==7) starting_delta_H_0=0.178;
     if(ch==8) starting_delta_H_0=0.275;
     if(ch==9) starting_delta_H_0=0.258;
     if(ch==10) starting_delta_H_0=0.221;
     if(ch==11) starting_delta_H_0=0.168;
     if(ch==12) starting_delta_H_0=0.273;
     if(ch==13) starting_delta_H_0=0.260;
     if(ch==14) starting_delta_H_0=0.221;
     if(ch==15) starting_delta_H_0=0.176;
     
     if(ch==0) starting_delta_T_0=0.315-0.280;
     if(ch==1) starting_delta_T_0=0.307-0.261;
     if(ch==2) starting_delta_T_0=0.301-0.221;
     if(ch==3) starting_delta_T_0=0.294-0.176;
     if(ch==4) starting_delta_T_0=0.315-0.280;
     if(ch==5) starting_delta_T_0=0.309-0.262;
     if(ch==6) starting_delta_T_0=0.303-0.222;
     if(ch==7) starting_delta_T_0=0.298-0.178;
     if(ch==8) starting_delta_T_0=0.318-0.275;
     if(ch==9) starting_delta_T_0=0.307-0.258;
     if(ch==10) starting_delta_T_0=0.304-0.221;
     if(ch==11) starting_delta_T_0=0.294-0.168;
     if(ch==12) starting_delta_T_0=0.321-0.273;
     if(ch==13) starting_delta_T_0=0.305-0.260;
     if(ch==14) starting_delta_T_0=0.301-0.221;
     if(ch==15) starting_delta_T_0=0.300-0.176;
   }
   */


   starting_alpha_CB=-1.5;//-0.5;
  starting_n_CB=4;

  low_alpha_CB=-5;
   up_alpha_CB=0.0;
   
  low_n_CB=0;
   up_n_CB=20;

  
  RooRealVar alpha_CB("alpha_CB","alpha parameter of CB",   starting_alpha_CB,low_alpha_CB,up_alpha_CB);
  RooRealVar     n_CB("n_CB",    "exponential decay of CB ",starting_n_CB,low_n_CB,up_n_CB);
  
  RooRealVar x_0("Time","Time [ns]",low_x_0,up_x_0);//eee
  /*
  RooRealVar mean_L_0("mean_L_0","mean of L gaussian background pos 0",starting_mean_L_0,low_mean_L_0,up_mean_L_0);
  RooRealVar Delta_H_0("Delta_H_0","Delta H pos 0",starting_delta_H_0,low_delta_H_0,up_delta_H_0);
  RooRealVar Delta_T_0("Delta_T_0","Delta T pos 0",starting_delta_T_0,low_delta_T_0,up_delta_T_0);
  RooFormulaVar mean_H_0("mean_H_0","mean_H_0","mean_L_0+Delta_H_0",RooArgList(mean_L_0,Delta_H_0));
  RooFormulaVar mean_T_0("mean_T_0","mean_T_0","mean_H_0+Delta_T_0",RooArgList(mean_H_0,Delta_T_0));
  */
  RooRealVar mean_H_0("mean_H_0","mean of L gaussian background pos 0",starting_mean_H_0,low_mean_H_0,up_mean_H_0);
  /*
  RooRealVar Delta_H_0("Delta_H_0","Delta H pos 0",starting_delta_H_0,low_delta_H_0,up_delta_H_0);
  RooRealVar Delta_T_0("Delta_T_0","Delta T pos 0",starting_delta_T_0,low_delta_T_0,up_delta_T_0);
  RooGaussian DElta_H_0("DElta_H_0","gaussian L_0",Delta_H_0,RooConst(mc_info.T_separation[ch][0]),RooConst(mc_info.err_T_separation[ch][0])) ;
  RooGaussian DElta_T_0("DElta_T_0","gaussian L_0",Delta_T_0,RooConst(mc_info.T_separation[ch][1]),RooConst(mc_info.err_T_separation[ch][1])) ;
  RooFormulaVar mean_L_0("mean_L_0","mean_l_0","mean_H_0-DElta_H_0",RooArgList(mean_H_0,DElta_H_0));
  RooFormulaVar mean_T_0("mean_T_0","mean_T_0","mean_H_0+DElta_T_0",RooArgList(mean_H_0,DElta_T_0));*/
  RooRealVar Delta_H_0("Delta_H_0","Delta H pos 0",starting_delta_H_0,low_delta_H_0,up_delta_H_0);
  RooRealVar Delta_T_0("Delta_T_0","Delta T pos 0",starting_delta_T_0,low_delta_T_0,up_delta_T_0);
  RooGaussian DElta_H_0("DElta_H_0","gaussian L_0",Delta_H_0,RooConst(mc_info.T_separation[ch][0]),RooConst(mc_info.err_T_separation[ch][0])) ;
  RooGaussian DElta_T_0("DElta_T_0","gaussian L_0",Delta_T_0,RooConst(mc_info.T_separation[ch][1]-mc_info.T_separation[ch][0]),RooConst(mc_info.err_T_separation[ch][1])) ;
  RooFormulaVar mean_L_0("mean_L_0","mean_l_0","mean_H_0-Delta_H_0",RooArgList(mean_H_0,Delta_H_0));
  RooFormulaVar mean_T_0("mean_T_0","mean_T_0","mean_H_0+Delta_T_0",RooArgList(mean_H_0,Delta_T_0));
  
  RooRealVar sigma_L_0("sigma_L_0","width of L gaussian background",starting_sigma_L_0,low_sigma_L_0,up_sigma_L_0);
  RooRealVar sigma_H_0("sigma_H_0","width of H gaussian background",starting_sigma_H_0,low_sigma_H_0,up_sigma_H_0);
  RooRealVar sigma_T_0("sigma_T_0","width of T gaussian background",starting_sigma_T_0,low_sigma_T_0,up_sigma_T_0);
  RooCBShape PDF_L_0("PDF_L_0","gaussian L_0",x,mean_L_0,sigma_L_0,alpha_CB,n_CB) ;
  RooCBShape PDF_H_0("PDF_H_0","gaussian H_0",x,mean_H_0,sigma_H_0,alpha_CB,n_CB) ;
  RooCBShape PDF_T_0("PDF_T_0","gaussian T_0",x,mean_T_0,sigma_T_0,alpha_CB,n_CB) ;  
  RooRealVar alpha_0("alpha_0","alpha_0",starting_alpha_0,low_alpha_0,up_alpha_0);
  RooRealVar beta_0("beta_0","beta_0",starting_beta_0,low_beta_0,up_beta_0);
  RooFormulaVar Frac_L_0("Frac_L_0","Frac_L_0","alpha_0",RooArgList(alpha_0));
  RooFormulaVar Frac_H_0("Frac_H_0","Frac_H_0","beta_0-alpha_0*beta_0",RooArgList(beta_0,alpha_0));
  RooFormulaVar Frac_T_0("Frac_T_0","Frac_T_0","1-alpha_0-beta_0+alpha_0*beta_0",RooArgList(alpha_0,beta_0));
  RooArgList  pdfList_sig_0(PDF_L_0,PDF_H_0); if(add_third_signal) {pdfList_sig_0.add(PDF_T_0);}//
  RooArgList  fracList_sig_0(alpha_0);if(add_third_signal) {fracList_sig_0.add(beta_0);}
  if(use_marginalize_method){
    RooAddPdf   PDF_sig_0_marginal("PDF_sig_0_marginal","PDF_sig_0_marginal",pdfList_sig_0,fracList_sig_0,kTRUE); //PDF(x|y)

    
    RooRealVar y("y","y",-2.5,2.5,"mm") ; //in mm
    // Create function w(y) = i/(bin width) //  PDF(y)
    RooRealVar a0("a0","a0",1.0/5) ;
    //RooPolyVar wy("wy","wy",y,RooArgSet(a0)) ;
    RooPolynomial wy("wy","wy",y,RooArgSet(a0)) ;

    // Create PDF(x|y)*w(y)
    RooProdPdf model("model","PDF(x|y)*w(y)",wy,Conditional(PDF_sig_0_marginal,x)) ;
  
    
    
    // M a r g i n a l i z e   m ( x , y )   t o   m ( x ) 
    // ----------------------------------------------------
    
    // modelx(x) = Int model(x,y) dy
    RooAbsPdf *PDF_sig_0 = model.createProjection(y) ;
    
  } else {
    RooAddPdf   PDF_sig_0("PDF_sig_0","PDF_sig_0",pdfList_sig_0,fracList_sig_0,kTRUE);
  }
  
  if(fix_deltas){Delta_H_0.setConstant(kTRUE); Delta_T_0.setConstant(kTRUE);}




  if(add_SP_components){
    double starting_mean_LL_SP_0;
    double low_mean_LL_SP_0;
    double up_mean_LL_SP_0;
    double starting_sigma_LL_SP_0;
    double low_sigma_LL_SP_0;
    double up_sigma_LL_SP_0;
    
    double starting_mean_L_SP_0;
    double low_mean_L_SP_0;
    double up_mean_L_SP_0;
    double starting_sigma_L_SP_0;
    double low_sigma_L_SP_0;
    double up_sigma_L_SP_0;
    double starting_mean_H_SP_0;
    double low_mean_H_SP_0;
    double up_mean_H_SP_0;
    double starting_sigma_H_SP_0;
    double low_sigma_H_SP_0;
    double up_sigma_H_SP_0;
    
    
    
    double starting_alpha_SP_0;
    double      low_alpha_SP_0;
    double       up_alpha_SP_0;
    double starting_beta_SP_0;
    double      low_beta_SP_0;
    double       up_beta_SP_0;
    
    
    
    low_sigma_LL_SP_0=0.035;
    up_sigma_LL_SP_0=.700;
    starting_sigma_LL_SP_0=0.30;
    
    low_sigma_L_SP_0=0.035;
    up_sigma_L_SP_0=2.00;
    starting_sigma_L_SP_0=0.40810;
    
    
    low_sigma_H_SP_0=0.035;
    up_sigma_H_SP_0=2.00;
    starting_sigma_H_SP_0=0.40810;
    
    starting_mean_LL_SP_0=24.5;
    low_mean_LL_SP_0=starting_mean_LL_SP_0-0.35;
    up_mean_LL_SP_0=starting_mean_LL_SP_0+0.35;
    
    starting_mean_L_SP_0=27;
    low_mean_L_SP_0=starting_mean_L_SP_0-0.35;
    up_mean_L_SP_0=starting_mean_L_SP_0+0.35;
    
    starting_mean_H_SP_0=32;
    low_mean_H_SP_0=starting_mean_H_SP_0-0.35;
    up_mean_H_SP_0=starting_mean_H_SP_0+0.35;
    
    
    
    
    low_alpha_SP_0=0.01;
    up_alpha_SP_0=0.99;
    starting_alpha_SP_0=0.8;
    
    low_beta_SP_0=0.2;
    up_beta_SP_0=0.99;
    starting_beta_SP_0=0.3;
    
    
     
    RooRealVar mean_LL_SP_0("mean_LL_SP_0","mean of L gaussian background pos 0",starting_mean_LL_SP_0,low_mean_LL_SP_0,up_mean_LL_SP_0);
    RooRealVar sigma_LL_SP_0("sigma_LL_SP_0","width of L gaussian background",starting_sigma_LL_SP_0,low_sigma_LL_SP_0,up_sigma_LL_SP_0);
    
    RooRealVar mean_L_SP_0("mean_L_SP_0","mean of L gaussian background pos 0",starting_mean_L_SP_0,low_mean_L_SP_0,up_mean_L_SP_0);
    RooRealVar sigma_L_SP_0("sigma_L_SP_0","width of L gaussian background",starting_sigma_L_SP_0,low_sigma_L_SP_0,up_sigma_L_SP_0);
    RooRealVar mean_H_SP_0("mean_H_SP_0","mean of L gaussian background pos 0",starting_mean_H_SP_0,low_mean_H_SP_0,up_mean_H_SP_0);
    RooRealVar sigma_H_SP_0("sigma_H_SP_0","width of L gaussian background",starting_sigma_H_SP_0,low_sigma_H_SP_0,up_sigma_H_SP_0);
    
    //RooCBShape PDF_LL_SP_0("PDF_LL_SP_0","gaussian L_SP_0",x,mean_LL_SP_0,sigma_LL_SP_0,alpha_CB,n_CB) ;
    //RooCBShape PDF_L_SP_0("PDF_L_SP_0","gaussian L_SP_0",x,mean_L_SP_0,sigma_L_SP_0,alpha_CB,n_CB) ;
    //RooCBShape PDF_H_SP_0("PDF_H_SP_0","gaussian L_SP_0",x,mean_H_SP_0,sigma_H_SP_0,alpha_CB,n_CB) ;

    
    RooLandau PDF_H_SP_0("PDF_H_SP_0","gaussian L_SP_0",x,mean_H_SP_0,sigma_H_SP_0);
    RooLandau PDF_L_SP_0("PDF_L_SP_0","gaussian L_SP_0",x,mean_L_SP_0,sigma_L_SP_0);
    RooLandau PDF_LL_SP_0("PDF_LL_SP_0","gaussian LL_SP_0",x,mean_LL_SP_0,sigma_LL_SP_0);
    
    RooRealVar alpha_SP_0("alpha_SP_0","alpha_SP_0",starting_alpha_SP_0,low_alpha_SP_0,up_alpha_SP_0);
    RooRealVar beta_SP_0("beta_SP_0","beta_SP_0",starting_beta_SP_0,low_beta_SP_0,up_beta_SP_0);
    
    RooArgList  pdfList_sig_SP_0(PDF_LL_SP_0,PDF_L_SP_0,PDF_H_SP_0); 
    RooArgList  fracList_sig_SP_0(beta_SP_0,alpha_SP_0);
    RooAddPdf   PDF_SP_0("PDF_SP_0","PDF_SP_0",pdfList_sig_SP_0,fracList_sig_SP_0,kTRUE);
  
  }




  
  RooRealVar a0_0("a0_0", "", 0.0, -10, 10);
  RooRealVar a1_0("a1_0", "", 0.0, -20, 20);
  RooRealVar a2_0("a2_0", "", 0.0015, -20, 20);
  RooChebychev PDF_B_0("PDF_B_0","PDF_B_0",x,RooArgList(a0_0,a1_0));//,a2_0));//
  
  RooRealVar  Frac_sig_0("Frac_sig_0","fraction of sig events", 0.9, 0.7,1.0);
  
  RooRealVar  Frac_SP_0("Frac_SP_0","fraction of secondary peaks events", 0.3, 0.1,0.99);
  if(use_marginalize_method){RooArgList  pdfList_0(*PDF_sig_0);}else{RooArgList  pdfList_0(PDF_sig_0);}
  //RooArgList  pdfList_0(PDF_sig_0);
  RooArgList  fracList_0(Frac_sig_0);
  if(add_SP_components){pdfList_0.add(PDF_SP_0);fracList_0.add(Frac_SP_0); }
  pdfList_0.add(PDF_B_0);
  if(gaussian_constraints_delta){
    //RooAddPdf   model_0("model_0","model_0",pdfList_0,fracList_0,kTRUE);
    //RooAddPdf   model_0_b("model_0_b","model_0_b",pdfList_0,fracList_0,kTRUE);
    RooAddPdf   Model_0("Model_0","Model_0",pdfList_0,fracList_0,kTRUE);
    RooArgList ListCOnstraints(Model_0,DElta_H_0); if(add_third_signal==true){ListCOnstraints.add(DElta_T_0);}
    RooProdPdf   model_0("model_0","model_0",ListCOnstraints);//RooArgList(Model_0,DElta_H_0));
    RooProdPdf   model_0_b("model_0_b","model_0_b",ListCOnstraints);//RooArgList(Model_0,DElta_H_0));
  }else{
    RooAddPdf   model_0("model_0","model_0",pdfList_0,fracList_0,kTRUE);
    RooAddPdf   model_0_b("model_0_b","model_0_b",pdfList_0,fracList_0,kTRUE);
  }
  
  
  
 
  
  
  
  
  
  
  
  
  // RooArgSet* model_params_0 = model_0.getParameters(x) ;
  //  model_params_0->Print("v") ;
  
  
  // I m p o r t    d a t a
  // ------------------------
  
  
  
  int Belle2_pixel;
  if(ch==0)  {Belle2_pixel=64 ;}
  if(ch==1)  {Belle2_pixel=128 ;}
  if(ch==2)  {Belle2_pixel=192 ;}
  if(ch==3)  {Belle2_pixel=256 ;}
  if(ch==4)  {Belle2_pixel=63 ;}
  if(ch==5)  {Belle2_pixel=127 ;}
  if(ch==6)  {Belle2_pixel=191 ;}
  if(ch==7)  {Belle2_pixel=255 ;}
  if(ch==8)  {Belle2_pixel=62 ;}
  if(ch==9)  {Belle2_pixel=126 ;}
  if(ch==10) {Belle2_pixel=190 ;}
  if(ch==11) {Belle2_pixel=254 ;}
  if(ch==12) {Belle2_pixel=61 ;}
  if(ch==13) {Belle2_pixel=125 ;}
  if(ch==14) {Belle2_pixel=189 ;}
  if(ch==15) {Belle2_pixel=253 ;}
  
  
  
  
  // D o  F i t 
  // ------------------------
  
  
  
  
  
  
  if(do_prefit){
    
    cout<<"___________________________________"<<endl;
    cout<<"                                   "<<endl;
    cout<<"-----------------------------------"<<endl;
    cout<<"***********************************"<<endl;
    cout<<"-----------------------------------"<<endl;
    cout<<"___________________________________"<<endl;
    cout<<"STATUS MIGRAD EDM 0 0"<<endl;
    
    //Delta_H_0.setConstant(kTRUE);
    sigma_L_0.setConstant(kTRUE);
    sigma_H_0.setConstant(kTRUE);
    RooFitResult* fit_results_0_b = model_0_b.fitTo(ds_0,Save(),Minimizer(Type_minim_pf,Algo_minim_pf),Strategy(2),SumW2Error(kFALSE),PrintLevel(MN_output_print_level_prefit),PrintEvalErrors(-1),Warnings(kFALSE),Verbose(kFALSE));//);
    sigma_L_0.setConstant(kFALSE);
    sigma_H_0.setConstant(kFALSE);
    //Delta_H_0.setConstant(kFALSE);
    if(print_prefit_info) fit_results_0_b->Print("v");
  }
  
  cout<<"___________________________________"<<endl;
  cout<<"                                   "<<endl;
  cout<<"-----------------------------------"<<endl;
  cout<<"***********************************"<<endl;
  cout<<"-----------------------------------"<<endl;
  cout<<"___________________________________"<<endl;
  cout<<"STATUS MIGRAD EDM 0"<<endl;
  
  
  RooArgList Ext_Constraints(DElta_H_0); if(add_third_signal){Ext_Constraints.add(DElta_T_0);}
  RooFitResult* fit_results_0;
  if(use_NLL){
    if(gaussian_constraints_delta) {
      RooAbsReal* nll_model_0 = model_0.createNLL(ds_0,Extended(kFALSE),ExternalConstraints(Ext_Constraints)) ;
    }else{
      RooAbsReal* nll_model_0 = model_0.createNLL(ds_0,Extended(kFALSE));
    }
    RooMinimizer RMN_0 = RooMinimizer(*nll_model_0);
    RMN_0.setErrorLevel(-1);
    RMN_0.setVerbose(kFALSE);
    RMN_0.setPrintEvalErrors(-1);
    RMN_0.setPrintLevel(MN_output_print_level);
    RMN_0.setStrategy(2);
    RMN_0.minimize(Type_minim,Algo_minim);
    fit_results_0=RMN_0.fit("hmr") ;
  }else{fit_results_0 = model_0.fitTo(ds_0,Save(),Strategy(2),SumW2Error(kFALSE),InitialHesse(true),PrintLevel(MN_output_print_level),PrintEvalErrors(-1),Warnings(kFALSE));//,InitialHesse(true));//,Hesse(kFALSE));//,Extended(kFALSE),Verbose(kFALSE));//Minimizer(Type_minim,Algo_minim),);
  }
  
  fit_results_0->Print("v");
  RooPlot* xframe2_0 = x.frame(Title(Form("pos. 0, channel %d",ch))) ;
  ds_0.plotOn(xframe2_0);//,DataError(RooAbsData::SumW2)) ;
  TH2 *h_correlation_0 = fit_results_0->correlationHist(); h_correlation_0->SetTitle(Form("correlation matrix 0 ch %d",ch));h_correlation_0->GetYaxis()->SetLabelSize(0.1); h_correlation_0->GetYaxis()->SetLabelFont(70); h_correlation_0->SetMarkerSize(2);
  model_0.plotOn(xframe2_0,Range("Fit_Range_0"),Components("PDF_L_0"),LineStyle(kDashed),LineColor(1)) ;
  model_0.plotOn(xframe2_0,Range("Fit_Range_0"),Components("PDF_H_0"),LineStyle(kDashed),LineColor(1)) ;
  if(add_third_signal) model_0.plotOn(xframe2_0,Range("Fit_Range_0"),Components("PDF_T_0"),LineStyle(kDashed),LineColor(1)) ;
  model_0.plotOn(xframe2_0,Range("Fit_Range_0"),Components("PDF_sig_0"),LineStyle(kDashed),LineColor(5)) ;
  model_0.plotOn(xframe2_0,Range("Fit_Range_0"),Components("PDF_B_0"),LineStyle(kDashed),LineColor(2)) ;
  if(add_SP_components){
    model_0.plotOn(xframe2_0,Range("Fit_Range_0"),Components("PDF_LL_SP_0"),LineStyle(kDashed),LineColor(kMagenta+4)) ;
    model_0.plotOn(xframe2_0,Range("Fit_Range_0"),Components("PDF_L_SP_0"),LineStyle(kDashed),LineColor(kMagenta+4)) ;
    model_0.plotOn(xframe2_0,Range("Fit_Range_0"),Components("PDF_H_SP_0"),LineStyle(kDashed),LineColor(kMagenta+4)) ;
  }
  model_0.plotOn(xframe2_0,Range("Fit_Range_0")) ;
  // Construct a histogram with the residuals of the data w.r.t. the curve
  RooHist* hresid_0 = xframe2_0->residHist() ;
  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist* hpull_0 = xframe2_0->pullHist() ;/*
  cout<<" 'aridaje"<<endl;
  xframe2_0->Print("V");
  RooCurve* curve1 = (RooCurve*) xframe2_0->getObject(6);  // 6
  RooHist* dataHist  = (RooHist*) xframe2_0->getHist("h_ds_0"); 
  RooHist* hpull_0 =  dataHist->makePullHist(*curve1,true);*/
  // Create a new frame to draw the pull distribution and add the distribution to the frame
  RooPlot* xframe4_0 = x.frame(Title(Form("Pull Distribution fiber 0 - ch %d",ch))) ;
  xframe4_0->addPlotable(hpull_0,"P") ;
  
  xframe2_0->GetXaxis()->SetLabelSize(0.05);
  xframe2_0->GetXaxis()->SetLabelFont(70);
  xframe2_0->GetXaxis()->SetTitleSize(0.05);
  xframe2_0->GetXaxis()->SetTitleFont(70);
  xframe4_0->GetXaxis()->SetLabelSize(0.05);
  xframe4_0->GetXaxis()->SetLabelFont(70);
  xframe4_0->GetXaxis()->SetTitleSize(0.05);
  xframe4_0->GetXaxis()->SetTitleFont(70);
  
  
  
  RooPlot* xframe2_0_log = x.frame(Title(Form("channel %d",ch))) ;
  ds_0.plotOn(xframe2_0_log);//,DataError(RooAbsData::SumW2)) ;
  model_0.plotOn(xframe2_0_log,Range("Fit_Range_0"),Components("PDF_L_0"),LineStyle(kDashed),LineColor(1)) ;
  model_0.plotOn(xframe2_0_log,Range("Fit_Range_0"),Components("PDF_H_0"),LineStyle(kDashed),LineColor(1)) ;
  if(add_third_signal) model_0.plotOn(xframe2_0_log,Range("Fit_Range_0"),Components("PDF_T_0"),LineStyle(kDashed),LineColor(1)) ;
  model_0.plotOn(xframe2_0_log,Range("Fit_Range_0"),Components("PDF_sig_0"),LineStyle(kDashed),LineColor(5)) ;
  model_0.plotOn(xframe2_0_log,Range("Fit_Range_0"),Components("PDF_B_0"),LineStyle(kDashed),LineColor(2)) ;
  model_0.plotOn(xframe2_0_log,Range("Fit_Range_0")) ;
  model_0.paramOn(xframe2_0_log);
  model_0.paramOn(xframe2_0);
  xframe2_0_log->GetXaxis()->SetLabelSize(0.05);
  xframe2_0_log->GetXaxis()->SetLabelFont(70);
  xframe2_0_log->GetXaxis()->SetTitleSize(0.05);
  xframe2_0_log->GetXaxis()->SetTitleFont(70);
  
  Fits_status.push_back(fit_results_0->status());
  
  
  alpha_CB.setVal(starting_alpha_CB);
  n_CB.setVal(starting_n_CB);
  
  
  
  double T0(-9);
  double T0_err(-9);
  
  if(!fit_highest_peak){
    T0=mean_L_0.getVal();
    T0_err=mean_L_0.getError();
  }else{
    T0=mean_H_0.getVal();
    T0_err=mean_H_0.getError();
    
  }
  
  ///////SINGLE FIBER FIT VALUES////////
  /*0   */   POIs.push_back(T0); 
  /*1   */   POIs.push_back(T0_err);
  /*2   */   POIs.push_back(Delta_H_0.getVal());
  /*3   */   POIs.push_back(Delta_H_0.getError());
  /*4   */   POIs.push_back(Delta_T_0.getVal());
  /*5   */   POIs.push_back(Delta_T_0.getError());
  /*6   */   POIs.push_back(sigma_L_0.getVal());
  /*7   */   POIs.push_back(sigma_L_0.getError());
  /*8   */   POIs.push_back(sigma_H_0.getVal());
  /*9   */   POIs.push_back(sigma_H_0.getError());
  /*10  */   POIs.push_back(sigma_T_0.getVal());
  /*11  */   POIs.push_back(sigma_T_0.getError());
  /*12  */   POIs.push_back(alpha_0.getVal());
  /*13  */   POIs.push_back(alpha_0.getError());
  /*14  */   POIs.push_back(beta_0.getVal());
  /*15  */   POIs.push_back(beta_0.getError());
  
  
  
  
  TString PAR_NAMES[8];
  PAR_NAMES[0]="mean_L";
  PAR_NAMES[1]="delta_H";
  PAR_NAMES[2]="delta_T";
  PAR_NAMES[3]="sigma_L";
  PAR_NAMES[4]="sigma_H";
  PAR_NAMES[5]="sigma_T";
  PAR_NAMES[6]="alpha";
  PAR_NAMES[7]="beta";
  
  
  
  
  cout<<"eeeih Delta H 0   ch : "<<ch<<"   "<<Delta_H_0.getVal()<<" +- "<<Delta_H_0.getError()<<endl;
  
  cout<<"<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
  cout<<"<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
  cout<<"<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
  cout<<" #### Fit quality checks, channel "<<ch<<" , fiber position 0 ####"<<endl;
  cout<<"------- General quality parameters ---------"<<endl;
  cout<<"Chi2/ndof  : "<<  xframe2_0->chiSquare()<<endl;
  cout<<"EDM : "<<fit_results_0->edm()<<endl;
  cout<<"Fit status : "<<fit_results_0->status()<<endl;
  cout<<"-------------  Hit boundaries --------------"<<endl;
  if(((POIs[2]-low_delta_H_0)<0.0001) || ((-POIs[2]+up_delta_H_0)<0.001))   cout<<"delta H_0 hits boundary"<<endl;
  //if(((POIs[4]-low_delta_T_0)<0.0001) || ((-POIs[4]+up_delta_T_0)<0.001))   cout<<"delta T_0 hits boundary"<<endl;
  if(((POIs[6]-low_sigma_L_0)<0.0001) || ((-POIs[6]+up_sigma_L_0)<0.001))   cout<<"sigma L_0 hits boundary"<<endl;
  if(((POIs[8]-low_sigma_H_0)<0.0001) || ((-POIs[8]+up_sigma_H_0)<0.001))   cout<<"sigma H_0 hits boundary"<<endl;
  //if(((POIs[10]-low_sigma_T_0)<0.0001)|| ((-POIs[10]+up_sigma_T_0)<0.001))  cout<<"sigma T_0 hits boundary"<<endl;
  if(((POIs[12]-low_alpha_0)<0.0001)  || ((-POIs[12]+up_alpha_0)<0.001))    cout<<"alpha_0  hits boundary"<<endl;
  //if(((POIs[14]-low_beta_0)<0.0001)   || ((-POIs[14]+up_beta_0)<0.001))     cout<<"beta_0 hits boundary"<<endl;
  cout<<"--------------- Large/small errors --------------"<<endl;
  for(int h=0;h<8;h++){
    if(h!=2&&h!=5&&h!=7){
      if(TMath::Abs(POIs[2*h+1]/POIs[2*h])>0.50)    cout<<"large error at "<<PAR_NAMES[h]<<"_0"<<endl;
      if(TMath::Abs(POIs[2*h+1]/POIs[2*h])<0.00001) cout<<"small error at "<<PAR_NAMES[h]<<"_0"<<endl;
    }
  }
  cout<<"--------------- Correlation matrix --------------"<<endl;
  for(int ff=1;ff<h_correlation_0->GetXaxis()->GetNbins()+1;ff++){
    for(int gg=1;gg<h_correlation_0->GetXaxis()->GetNbins()-ff+1;gg++){
      if(TMath::Abs(h_correlation_0->GetBinContent(ff,gg))<0.001||TMath::Abs(h_correlation_0->GetBinContent(ff,gg))>0.99){
	cout<<"low/high "<<h_correlation_0->GetXaxis()->GetBinLabel(ff)<<"-"<<h_correlation_0->GetYaxis()->GetBinLabel(gg)<<" correlation : "<<h_correlation_0->GetBinContent(ff,gg)<< endl;
      }
    }
  }
  cout<<"   "<<endl;
  cout<<"***************************************************************************"<<endl;
  
  
  
  
  TF1 *line_0s = new TF1("line_0s","0",-1,16);line_0s->SetLineColor(8);
  TF1 *line_1ps = new TF1("line_1ps","1",-1,16);line_1ps->SetLineColor(4);
  TF1 *line_1ns = new TF1("line_1ns","-1",-1,16);line_1ns->SetLineColor(4);
  TF1 *line_2ps = new TF1("line_2ps","2",-1,16);line_2ps->SetLineColor(kOrange-2);
  TF1 *line_2ns = new TF1("line_2ns","-2",-1,16);line_2ns->SetLineColor(kOrange-2);
  TF1 *line_3ps = new TF1("line_3ps","3",-1,16);line_3ps->SetLineColor(2);
  TF1 *line_3ns = new TF1("line_3ns","-3",-1,16);line_3ns->SetLineColor(2);
  if(draw_results){  
    TCanvas* c_Fit = new TCanvas("Fit results","Fit results",0,0,1124,700) ;
    c_Fit->Divide(1,3) ;
    gStyle->SetOptFit(0111); 
    c_Fit->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe2_0->GetYaxis()->SetTitleOffset(1.6) ; xframe2_0->Draw() ;
    //c_Fit->cd(4) ; gPad->SetLeftMargin(0.15) ; xframe2_0_log->GetYaxis()->SetTitleOffset(1.6) ; gPad->SetLogy(1); xframe2_0_log->Draw() ;
    c_Fit->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe4_0->GetYaxis()->SetTitleOffset(1.6) ; gPad->SetGridy(); xframe4_0->GetYaxis()->SetRangeUser(-5,5); xframe4_0->Draw() ;
    line_3ns->Draw("same");line_3ps->Draw("same");
    line_2ns->Draw("same");line_2ps->Draw("same");
    line_1ns->Draw("same");line_1ps->Draw("same");
    line_0s->Draw("same");line_0s->Draw("same");
    c_Fit->cd(3) ; gPad->SetLeftMargin(0.15) ; h_correlation_0->Draw("colz:text");
    
    
    TCanvas* c_Fit_0 = new TCanvas("c_Fit_0","c_Fit_0");
    xframe2_0->Draw() ;
  }
  
  
  
  
  //my_fit_results.Results=POIs;
  for(int g=0; g<POIs.size();g++){
    my_fit_results.Results[g]=POIs[g];
  }
  my_fit_results.Overall_Fracs[0]=Frac_sig_0.getVal();
  my_fit_results.Overall_Fracs[1]=Frac_sig_0.getError();
  if(add_SP_components) {
    my_fit_results.Overall_Fracs[2]=Frac_SP_0.getVal();
    my_fit_results.Overall_Fracs[3]=Frac_SP_0.getError();
    my_fit_results.SP_separation[0]=mean_L_SP_0.getVal()-T0;
    my_fit_results.SP_separation[1]=mean_L_SP_0.getError()+T0_err;
    my_fit_results.SP_separation[2]=mean_H_SP_0.getVal()-mean_L_SP_0.getVal();
    my_fit_results.SP_separation[3]=mean_H_SP_0.getError()+mean_L_SP_0.getError();
  }else{
    my_fit_results.Overall_Fracs[2]=-9;
    my_fit_results.Overall_Fracs[3]=-9;
    my_fit_results.SP_separation[0]=-9;
    my_fit_results.SP_separation[1]=-9;
    my_fit_results.SP_separation[2]=-9;
    my_fit_results.SP_separation[3]=-9;
  }
  my_fit_results.n_POIs=POIs.size();
  my_fit_results.xframe2_amplitude_0=xframe2_0_amp;
  my_fit_results.xframe2_fit_0=xframe2_0;
  my_fit_results.xframe2_pull_0=xframe4_0;
  my_fit_results.h_correlation_0=h_correlation_0;
  my_fit_results.most_probable_amp[0]=TMath::Abs(h_amp_histogram_0->ProjectionY()->GetMaximumBin());
  my_fit_results.most_probable_amp[1]=h_amp_histogram_0->ProjectionY()->GetRMS();
  
  
  //  delete  h_input_histogram_0;
  
  
  return my_fit_results;
  
  //return POIs;
}



vector<float> loop_channels(int deep_fixed_params,TString input_tune, bool plot_summaries){ //rel_weight is the relative weight of the 1 to 0 pois, i.e. total dataset= pos.0 + rel_weight*pos.1 - in the real case rel_weight=1.0
  
  
  TCanvas *c_pos0_AllChannels = new TCanvas("c_pos0_AllChannels","c_pos0_AllChannels",0,0,1124,700);
  TCanvas *c_pos0_AllChannels_pulls = new TCanvas("c_pos0_AllChannels_pulls","c_pos0_AllChannels_pulls",0,0,1124,700);
  TCanvas *c_pos0_AllChannels_corr = new TCanvas("c_pos0_AllChannels_corr","c_pos0_AllChannels_corr",0,0,1124,700);
  TCanvas *c_pos0_AllChannels_amp = new TCanvas("c_pos0_AllChannels_amp","c_pos0_AllChannels_amp",0,0,1124,700);
  
  
  
  vector<float> relative_weight_Frac_0; relative_weight_Frac_0.clear();
  c_pos0_AllChannels->Divide(4,4);
  
  c_pos0_AllChannels_pulls->Divide(4,4);
  c_pos0_AllChannels_corr->Divide(4,4);
  c_pos0_AllChannels_amp->Divide(4,4);
  
  
  
  
  int index_channel_pixel[16];
  index_channel_pixel[0]=1 ;
  index_channel_pixel[1]=5 ;
  index_channel_pixel[2]=9 ;
  index_channel_pixel[3]=13 ;
  index_channel_pixel[4]=2 ;
  index_channel_pixel[5]=6 ;
  index_channel_pixel[6]=10 ;
  index_channel_pixel[7]=14 ;
  index_channel_pixel[8]=3 ;
  index_channel_pixel[9]=7 ;
  index_channel_pixel[10]=11 ;
  index_channel_pixel[11]=15 ;
  index_channel_pixel[12]=4 ;
  index_channel_pixel[13]=8 ;
  index_channel_pixel[14]=12 ;
  index_channel_pixel[15]=16 ;
  
  
  
  
  
  
  vector<float> POIs;
  vector<float> FRCs;
  vector<float> AMPs;
  vector<float> SPSs;//Secondary peaks separations
  float mean_L_0_B[16];
  float err_mean_L_0_B[16];
  float Delta_H_0_B[16];
  float err_Delta_H_0_B[16];
  float Delta_T_0_B[16];
  float err_Delta_T_0_B[16];
  float sigma_L_0_B[16];
  float err_sigma_L_0_B[16];
  float sigma_H_0_B[16];
  float err_sigma_H_0_B[16];
  float sigma_T_0_B[16];
  float err_sigma_T_0_B[16];
  float mean_L_0_A[16];
  float err_mean_L_0_A[16];
  float Delta_H_0_A[16];
  float err_Delta_H_0_A[16];
  float Delta_T_0_A[16];
  float err_Delta_T_0_A[16];
  float sigma_L_0_A[16];
  float err_sigma_L_0_A[16];
  float sigma_H_0_A[16];
  float err_sigma_H_0_A[16];
  float sigma_T_0_A[16];
  float err_sigma_T_0_A[16];
  
  float mean_L_1_B[16];
  float err_mean_L_1_B[16];
  float Delta_H_1_B[16];
  float err_Delta_H_1_B[16];
  float Delta_T_1_B[16];
  float err_Delta_T_1_B[16];
  float sigma_L_1_B[16];
  float err_sigma_L_1_B[16];
  float sigma_T_1_B[16];
  float err_sigma_T_1_B[16];
  float sigma_H_1_B[16];
  float err_sigma_H_1_B[16];
  float mean_L_1_A[16];
  float err_mean_L_1_A[16];
  float Delta_H_1_A[16];
  float err_Delta_H_1_A[16];
  float Delta_T_1_A[16];
  float err_Delta_T_1_A[16];
  float sigma_L_1_A[16];
  float err_sigma_L_1_A[16];
  float sigma_H_1_A[16];
  float err_sigma_H_1_A[16];
  float sigma_T_1_A[16];
  float err_sigma_T_1_A[16];
  
  float frac_L_B_0[16];      
  float err_frac_L_B_0[16];  
  float frac_L_A_0[16];      
  float err_frac_L_A_0[16];  
  float frac_L_B_1[16];      
  float err_frac_L_B_1[16]; 
  float frac_L_A_1[16];     
  float err_frac_L_A_1[16]; 
  
  float frac_H_B_0[16];   
  float err_frac_H_B_0[16]; 
  float frac_H_A_0[16];
  float err_frac_H_A_0[16];
  float frac_H_B_1[16];
  float err_frac_H_B_1[16];
  float frac_H_A_1[16];
  float err_frac_H_A_1[16];
  
  float frac_0[16];
  float err_frac_0[16];

  float frac_Sig_0[16];
  float err_frac_Sig_0[16];
  float frac_SP_0[16];
  float err_frac_SP_0[16];

  float delta_LSP_T0[16];
  float err_delta_LSP_T0[16];
  float delta_HSP_LSP[16];
  float err_delta_HSP_LSP[16];

    
  float amp_Sig_0[16];
  float err_amp_Sig_0[16];
  
  float ref_delta_mean[16];
  ref_delta_mean[0]=0.339;
  ref_delta_mean[1]=0.304;
  ref_delta_mean[2]=0.259;
  ref_delta_mean[3]=0.198;
  ref_delta_mean[4]=0.329;
  ref_delta_mean[5]=0.295;
  ref_delta_mean[6]=0.251;
  ref_delta_mean[7]=0.194;
  ref_delta_mean[8]=0.334;
  ref_delta_mean[9]=0.304;
  ref_delta_mean[10]=0.271;
  ref_delta_mean[11]=0.185;
  ref_delta_mean[12]=0.317;
  ref_delta_mean[13]=0.308;
  ref_delta_mean[14]=0.285;
  ref_delta_mean[15]=0.23;
  float ref_sigma_L[16];
  ref_sigma_L[0]= 0.079;
  ref_sigma_L[1]= 0.084;
  ref_sigma_L[2]= 0.076;
  ref_sigma_L[3]= 0.12;
  ref_sigma_L[4]= 0.069;
  ref_sigma_L[5]= 0.080;
  ref_sigma_L[6]= 0.081;
  ref_sigma_L[7]= 0.095;
  ref_sigma_L[8]= 0.086;
  ref_sigma_L[9]= 0.090;
  ref_sigma_L[10]= 0.087;
  ref_sigma_L[11]= 0.082;
  ref_sigma_L[12]= 0.083;
  ref_sigma_L[13]= 0.093;
  ref_sigma_L[14]= 0.065;
  ref_sigma_L[15]= 0.076;
  
  float ref_sigma_H[16];
  ref_sigma_H[0]=0.094;
  ref_sigma_H[1]=  0.105;
  ref_sigma_H[2]=  0.092;
  ref_sigma_H[3]=  0.104;
  ref_sigma_H[4]=  0.083;
  ref_sigma_H[5]=  0.080;
  ref_sigma_H[6]=  0.077;
  ref_sigma_H[7]=  0.095;
  ref_sigma_H[8]=  0.096;
  ref_sigma_H[9]=  0.083;
  ref_sigma_H[10]=  0.088;
  ref_sigma_H[11]=  0.097;
  ref_sigma_H[12]=  0.082;
  ref_sigma_H[13]=  0.080;
  ref_sigma_H[14]=  0.090;
  ref_sigma_H[15]=  0.049;
  
  float ref_sigma_T[16];
  ref_sigma_T[0]=0.323;
  ref_sigma_T[1]=0.379;
  ref_sigma_T[2]=0.345;
  ref_sigma_T[3]=0.357;
  ref_sigma_T[4]=0.324;
  ref_sigma_T[5]=0.325;
  ref_sigma_T[6]=0.325;
  ref_sigma_T[7]=0.336;
  ref_sigma_T[8]=0.347;
  ref_sigma_T[9]=0.360;
  ref_sigma_T[10]=0.330;
  ref_sigma_T[11]=0.328;
  ref_sigma_T[12]=0.322;
  ref_sigma_T[13]=0.324;
  ref_sigma_T[14]=0.352;
  ref_sigma_T[15]=0.247;
  
  
  float ref_pos0_MC_B2[16];
  ref_pos0_MC_B2[0]=0.284;
  ref_pos0_MC_B2[1]=0.257;
  ref_pos0_MC_B2[2]=0.224;
  ref_pos0_MC_B2[3]=0.187;
  ref_pos0_MC_B2[4]=0.291;
  ref_pos0_MC_B2[5]=0.261;
  ref_pos0_MC_B2[6]=0.224;
  ref_pos0_MC_B2[7]=0.187;
  ref_pos0_MC_B2[8]=0.277;
  ref_pos0_MC_B2[9]=0.255;
  ref_pos0_MC_B2[10]=0.222;
  ref_pos0_MC_B2[11]=0.183;
  ref_pos0_MC_B2[12]=0.281;
  ref_pos0_MC_B2[13]=0.255;
  ref_pos0_MC_B2[14]=0.219;
  ref_pos0_MC_B2[15]=0.184;
  
  
  float ref_pos0_MC_PD[16];
  ref_pos0_MC_PD[0]=0.311;
  ref_pos0_MC_PD[1]=0.269;
  ref_pos0_MC_PD[2]=0.223;
  ref_pos0_MC_PD[3]=0.177;
  ref_pos0_MC_PD[4]=0.314;
  ref_pos0_MC_PD[5]=0.270;
  ref_pos0_MC_PD[6]=0.224;
  ref_pos0_MC_PD[7]=0.178;
  ref_pos0_MC_PD[8]=0.316;
  ref_pos0_MC_PD[9]=0.273;
  ref_pos0_MC_PD[10]=0.224;
  ref_pos0_MC_PD[11]=0.178;
  ref_pos0_MC_PD[12]=0.318;
  ref_pos0_MC_PD[13]=0.271;
  ref_pos0_MC_PD[14]=0.224;
  ref_pos0_MC_PD[15]=0.178;
  
  
  
  
  
  
  
  
  float ref_pos1_MC_B2[16];
  ref_pos1_MC_B2[0]=  0.302;
  ref_pos1_MC_B2[1]=0.280;
  ref_pos1_MC_B2[2]=0.225;
  ref_pos1_MC_B2[3]=0.187;
  ref_pos1_MC_B2[4]=0.297;
  ref_pos1_MC_B2[5]=0.273;
  ref_pos1_MC_B2[6]=0.225;
  ref_pos1_MC_B2[7]=0.189;
  ref_pos1_MC_B2[8]=0.295;
  ref_pos1_MC_B2[9]=0.266;
  ref_pos1_MC_B2[10]=0.224;
  ref_pos1_MC_B2[11]=0.188;
  ref_pos1_MC_B2[12]=0.295;
  ref_pos1_MC_B2[13]=0.263;
  ref_pos1_MC_B2[14]=0.223;
  ref_pos1_MC_B2[15]=0.186;
  
  
  
  ////////////////////////FIT RESULTS WITH 3 GAUSSIAN SIGNALS/////////////////
  
  float ref_pos0_3G[16];
  float err_ref_pos0_3G[16];
  ref_pos0_3G[0]=0.328265 ;err_ref_pos0_3G[0]=  0.0133491;
  ref_pos0_3G[1]=0.310638 ;err_ref_pos0_3G[1]=  0.00705744;
  ref_pos0_3G[2]=0.270493 ;err_ref_pos0_3G[2]=  0.010975;
  ref_pos0_3G[3]=0.216288 ;err_ref_pos0_3G[3]=  0.0158163;
  ref_pos0_3G[4]=0.337205 ;err_ref_pos0_3G[4]=  0.0106956;
  ref_pos0_3G[5]=0.277298 ;err_ref_pos0_3G[5]=  0.00802361;
  ref_pos0_3G[6]=0.254835 ;err_ref_pos0_3G[6]=  0.00799138;
  ref_pos0_3G[7]=0.2075   ;err_ref_pos0_3G[7]=  0.00316242;
  ref_pos0_3G[8]=0.34383  ;err_ref_pos0_3G[8]=  0.0280523;
  ref_pos0_3G[9]=0.295537 ;err_ref_pos0_3G[9]=  0.00905063;
  ref_pos0_3G[10]=0.275667;err_ref_pos0_3G[10]=  0.00886905;
  ref_pos0_3G[11]=0.2075  ;err_ref_pos0_3G[11]=  0.00267526;
  ref_pos0_3G[12]=0.335146;err_ref_pos0_3G[12]=  0.012047;
  ref_pos0_3G[13]=0.310555;err_ref_pos0_3G[13]=  0.010663;
  ref_pos0_3G[14]=0.284698;err_ref_pos0_3G[14]=  0.00878808;
  ref_pos0_3G[15]=0.232241;err_ref_pos0_3G[15]=  0.0150502;
  
  float ref_pos1_3G[16];
  float err_ref_pos1_3G[16];
  ref_pos1_3G[0]=0.339691  ; err_ref_pos1_3G[0]=  0.00681139;
  ref_pos1_3G[1]=0.313315  ; err_ref_pos1_3G[1]=  0.00903739;
  ref_pos1_3G[2]=0.30621   ; err_ref_pos1_3G[2]=  0.00617384;
  ref_pos1_3G[3]=0.288669  ; err_ref_pos1_3G[3]=  0.00976224;
  ref_pos1_3G[4]=0.342045  ; err_ref_pos1_3G[4]=  0.00592796;
  ref_pos1_3G[5]=0.308783  ; err_ref_pos1_3G[5]=  0.00497152;
  ref_pos1_3G[6]=0.305168  ; err_ref_pos1_3G[6]=  0.00559984;
  ref_pos1_3G[7]=0.269488  ; err_ref_pos1_3G[7]=  0.00745161;
  ref_pos1_3G[8]=0.332781  ; err_ref_pos1_3G[8]=  0.00738368;
  ref_pos1_3G[9]=0.32234   ; err_ref_pos1_3G[9]=  0.00845406;
  ref_pos1_3G[10]=0.286107 ; err_ref_pos1_3G[10]=  0.0116438;
  ref_pos1_3G[11]=0.278205 ; err_ref_pos1_3G[11]=  0.010652;
  ref_pos1_3G[12]=0.318575 ; err_ref_pos1_3G[12]=  0.00519016;
  ref_pos1_3G[13]=0.300113 ; err_ref_pos1_3G[13]=  0.0117054;
  ref_pos1_3G[14]=0.289259 ; err_ref_pos1_3G[14]=  0.0127598;
  ref_pos1_3G[15]=0.226567 ; err_ref_pos1_3G[15]=  0.0197404;
  
  
  ////////////////////////////////////////////////////////////////////////////
  
  float delta_ref_sigma_H[16];
  float delta_ref_sigma_L[16];
  
  
  float x[16];
  float err_x[16];
  
  Fits_status.clear();
  for(Int_t i=0; i<16; i++){
    x[i]=i; err_x[i]=0;
    POIs.clear();
    FRCs.clear();
    AMPs.clear();
    SPSs.clear();
    cout<<"Channel : "<<i<<endl;
    Fit_results my_fit_Results;
    
    my_fit_Results=Fit_head("blind",deep_fixed_params,input_tune,i);
    for(int hg=0;hg<my_fit_Results.n_POIs;hg++){
      POIs.push_back(my_fit_Results.Results[hg]);
    }
    SPSs.push_back(my_fit_Results.SP_separation[0]);
    SPSs.push_back(my_fit_Results.SP_separation[1]);
    SPSs.push_back(my_fit_Results.SP_separation[2]);
    SPSs.push_back(my_fit_Results.SP_separation[3]);
    FRCs.push_back(my_fit_Results.Overall_Fracs[0]); 
    FRCs.push_back(my_fit_Results.Overall_Fracs[1]); 
    FRCs.push_back(my_fit_Results.Overall_Fracs[2]); 
    FRCs.push_back(my_fit_Results.Overall_Fracs[3]);
    AMPs.push_back(my_fit_Results.most_probable_amp[0]); 
    AMPs.push_back(my_fit_Results.most_probable_amp[1]); 
    
    
    TF1 *line_0s = new TF1("line_0s","0",-1,16);line_0s->SetLineColor(8);
    TF1 *line_1ps = new TF1("line_1ps","1",-1,16);line_1ps->SetLineColor(4);
    TF1 *line_1ns = new TF1("line_1ns","-1",-1,16);line_1ns->SetLineColor(4);
    TF1 *line_2ps = new TF1("line_2ps","2",-1,16);line_2ps->SetLineColor(kOrange-2);
    TF1 *line_2ns = new TF1("line_2ns","-2",-1,16);line_2ns->SetLineColor(kOrange-2);
    TF1 *line_3ps = new TF1("line_3ps","3",-1,16);line_3ps->SetLineColor(2);
    TF1 *line_3ns = new TF1("line_3ns","-3",-1,16);line_3ns->SetLineColor(2);
    c_pos0_AllChannels->cd(index_channel_pixel[i]); my_fit_Results.xframe2_fit_0->Draw(); 
    c_pos0_AllChannels_pulls->cd(index_channel_pixel[i]); gPad->SetGridy(); my_fit_Results.xframe2_pull_0->GetYaxis()->SetRangeUser(-5,5); my_fit_Results.xframe2_pull_0->Draw();
    line_3ns->Draw("same");line_3ps->Draw("same");
    line_2ns->Draw("same");line_2ps->Draw("same");
    line_1ns->Draw("same");line_1ps->Draw("same");
    line_0s->Draw("same");line_0s->Draw("same");
    c_pos0_AllChannels_corr->cd(index_channel_pixel[i]);my_fit_Results.h_correlation_0->Draw("colz:text");
    c_pos0_AllChannels_amp->cd(index_channel_pixel[i]);my_fit_Results.xframe2_amplitude_0->Draw();
    
    
    
    
    mean_L_0_B[i]      =POIs[0];
    err_mean_L_0_B[i]  =POIs[1];
    sigma_L_0_B[i]     =POIs[6];
    err_sigma_L_0_B[i] =POIs[7];
    Delta_H_0_B[i]     =POIs[2];
    err_Delta_H_0_B[i] =POIs[3];
    sigma_H_0_B[i]     =POIs[8];
    err_sigma_H_0_B[i] =POIs[9];
    Delta_T_0_B[i]     =POIs[4];
    err_Delta_T_0_B[i] =POIs[5];
    sigma_T_0_B[i]     =POIs[10];
    err_sigma_T_0_B[i] =POIs[11];
 
    
    frac_L_B_0[i]         =POIs[12];
    err_frac_L_B_0[i]     =POIs[13];
    
    
    frac_H_B_0[i]         =POIs[14];
    err_frac_H_B_0[i]     =POIs[15];

    frac_Sig_0[i]      = FRCs[0];
    err_frac_Sig_0[i]  = FRCs[1];
    frac_SP_0[i]       = FRCs[2];
    err_frac_SP_0[i]   = FRCs[3];

    delta_LSP_T0[i]         = SPSs[0];
    err_delta_LSP_T0[i]     = SPSs[1];
    delta_HSP_LSP[i]        = SPSs[2];
    err_delta_HSP_LSP[i]    = SPSs[3];
    

    amp_Sig_0[i]     = AMPs[0];
    err_amp_Sig_0[i] = AMPs[1];
  }
  
  
  if(plot_summaries){
    
    
    
    TF1 *line = new TF1("line","0.100",-1,16);
    line->SetTitle("100 ps");
    // TLine *line = new TLine(0,100,16,100);
    line->SetLineColor(2);
    TGraphErrors *DELTA_means_ref = new TGraphErrors(16,x,ref_delta_mean,err_x,err_x);DELTA_means_ref->SetName("DELTA_means_ref"); DELTA_means_ref->SetTitle("#Delta t_{ref} pos.0 - fit pos.0");
    TGraphErrors *DELTA_means_pos0_MC_B2 = new TGraphErrors(16,x,ref_pos0_MC_B2,err_x,err_x);DELTA_means_pos0_MC_B2->SetName("DELTA_means_pos0_MC_B2"); DELTA_means_pos0_MC_B2->SetTitle("#Delta t_{ref} pos.0 - B2 MC");
    TGraphErrors *DELTA_means_pos1_MC_B2 = new TGraphErrors(16,x,ref_pos1_MC_B2,err_x,err_x);DELTA_means_pos1_MC_B2->SetName("DELTA_means_pos1_MC_B2"); DELTA_means_pos1_MC_B2->SetTitle("#Delta t_{ref} pos.1 - B2 MC");
    TGraphErrors *DELTA_means_pos0_MC_PD = new TGraphErrors(16,x,ref_pos0_MC_PD,err_x,err_x);DELTA_means_pos0_MC_PD->SetName("DELTA_means_pos0_MC_PD"); DELTA_means_pos0_MC_PD->SetTitle("#Delta t_{ref} pos.0 - PD MC");
    
    TGraphErrors *DELTA_means_pos0_3G = new TGraphErrors(16,x,ref_pos0_3G,err_x,err_ref_pos0_3G);DELTA_means_pos0_3G->SetName("DELTA_means_pos0_3G"); DELTA_means_pos0_3G->SetTitle("#Delta t_{ref} pos.0 - 3 Gaussian model");
    TGraphErrors *DELTA_means_pos1_3G = new TGraphErrors(16,x,ref_pos1_3G,err_x,err_ref_pos1_3G);DELTA_means_pos1_3G->SetName("DELTA_means_pos1_3G"); DELTA_means_pos1_3G->SetTitle("#Delta t_{ref} pos.1 - 3 Gaussian model");
    
    
    
    TGraphErrors *FRAC_SIG = new TGraphErrors(16,x,frac_Sig_0,err_x,err_frac_Sig_0);FRAC_SIG->SetName("FRAC_SIG");FRAC_SIG->SetTitle("F_{sig}");
    TGraphErrors *FRAC_SP = new TGraphErrors(16,x,frac_SP_0,err_x,err_frac_SP_0);FRAC_SP->SetName("FRAC_SP");FRAC_SP->SetTitle("F_{SP}");
    TGraphErrors *AMP_SIG = new TGraphErrors(16,x,amp_Sig_0,err_x,err_amp_Sig_0);AMP_SIG->SetName("AMP_SIG");AMP_SIG->SetTitle("Amp_{Sig}");
    
    TGraphErrors *DELTA_LSP_T0 = new TGraphErrors(16,x,delta_LSP_T0,err_x,err_delta_LSP_T0);   DELTA_LSP_T0->SetName("DELTA_LSP_T0");  DELTA_LSP_T0->SetTitle("T^{SP}_{L}-T_{0}");
    TGraphErrors *DELTA_HSP_LSP = new TGraphErrors(16,x,delta_HSP_LSP,err_x,err_delta_HSP_LSP);DELTA_HSP_LSP->SetName("DELTA_HSP_LSP");DELTA_HSP_LSP->SetTitle("T^{SP}_{H}-T^{SP}_{L}");


    
    TGraphErrors *MEAN_L_B_0 = new TGraphErrors(16,x,mean_L_0_B,err_x,err_mean_L_0_B);MEAN_L_B_0->SetName("MEAN_L_B_0"); MEAN_L_B_0->SetTitle("T_{0} - fit pos.0");
    TGraphErrors *MEAN_L_B_1 = new TGraphErrors(16,x,mean_L_1_B,err_x,err_mean_L_1_B);MEAN_L_B_1->SetName("MEAN_L_B_1"); MEAN_L_B_1->SetTitle("#Delta T_{L} #equiv T_{L}^{pos.1}-T_{L}^{pos.0} T- fit pos.1");
    TGraphErrors *MEAN_L_A_0 = new TGraphErrors(16,x,mean_L_0_A,err_x,err_mean_L_0_A);MEAN_L_A_0->SetName("MEAN_L_A_0"); MEAN_L_A_0->SetTitle("T_{L}^{pos.0} - fit pos.0 #oplus 1");
    TGraphErrors *MEAN_L_A_1 = new TGraphErrors(16,x,mean_L_1_A,err_x,err_mean_L_1_A);MEAN_L_A_1->SetName("MEAN_L_A_1"); MEAN_L_A_1->SetTitle("#Delta T_{L} #equiv T_{L}^{pos.1}-T_{L}^{pos.0} - fit pos.0 #oplus 1");
    
    TGraphErrors *DELTA_H_B_0 = new TGraphErrors(16,x,Delta_H_0_B,err_x,err_Delta_H_0_B);DELTA_H_B_0->SetName("DELTA_H_B_0"); DELTA_H_B_0->SetTitle("#Delta t pos.0 - fit pos.0");
    TGraphErrors *DELTA_H_B_1 = new TGraphErrors(16,x,Delta_H_1_B,err_x,err_Delta_H_1_B);DELTA_H_B_1->SetName("DELTA_H_B_1"); DELTA_H_B_1->SetTitle("#Delta t pos.1 - fit pos.1");
    TGraphErrors *DELTA_H_A_0 = new TGraphErrors(16,x,Delta_H_0_A,err_x,err_Delta_H_0_A);DELTA_H_A_0->SetName("DELTA_H_A_0"); DELTA_H_A_0->SetTitle("#Delta t_{LOW-HIGH} pos.0 - fit pos.0 #oplus 1");
    TGraphErrors *DELTA_H_A_1 = new TGraphErrors(16,x,Delta_H_1_A,err_x,err_Delta_H_1_A);DELTA_H_A_1->SetName("DELTA_H_A_1"); DELTA_H_A_1->SetTitle("#Delta t_{LOW-HIGH} pos.1 - fit pos.0 #oplus 1");
    TGraphErrors *DELTA_T_B_0 = new TGraphErrors(16,x,Delta_T_0_B,err_x,err_Delta_T_0_B);DELTA_T_B_0->SetName("DELTA_T_B_0"); DELTA_T_B_0->SetTitle("#Delta t_{HIGH-THIRD} pos.0 - fit pos.0");
    TGraphErrors *DELTA_T_B_1 = new TGraphErrors(16,x,Delta_T_1_B,err_x,err_Delta_T_1_B);DELTA_T_B_1->SetName("DELTA_T_B_1"); DELTA_T_B_1->SetTitle("#Delta t_{HIGH-THIRD} pos.1 - fit pos.1");
    TGraphErrors *DELTA_T_A_0 = new TGraphErrors(16,x,Delta_T_0_A,err_x,err_Delta_T_0_A);DELTA_T_A_0->SetName("DELTA_T_A_0"); DELTA_T_A_0->SetTitle("#Delta t_{HIGH-THIRD} pos.0 - fit pos.0 #oplus 1");
    TGraphErrors *DELTA_T_A_1 = new TGraphErrors(16,x,Delta_T_1_A,err_x,err_Delta_T_1_A);DELTA_T_A_1->SetName("DELTA_T_A_1"); DELTA_T_A_1->SetTitle("#Delta t_{HIGH-THIRD} pos.1 - fit pos.0 #oplus 1");
    
    TGraphErrors *SIGMA_L_B_0 = new TGraphErrors(16,x,sigma_L_0_B,err_x,err_sigma_L_0_B);SIGMA_L_B_0->SetName("SIGMA_L_B_0"); SIGMA_L_B_0->SetTitle("#delta t_{low} pos.0 - fit pos.0");
    TGraphErrors *SIGMA_L_B_1 = new TGraphErrors(16,x,sigma_L_1_B,err_x,err_sigma_L_1_B);SIGMA_L_B_1->SetName("SIGMA_L_B_1"); SIGMA_L_B_1->SetTitle("#delta t_{low} pos.1 - fit pos.1");
    TGraphErrors *SIGMA_H_B_0 = new TGraphErrors(16,x,sigma_H_0_B,err_x,err_sigma_H_0_B);SIGMA_H_B_0->SetName("SIGMA_H_B_0"); SIGMA_H_B_0->SetTitle("#delta t_{high} pos.0 - fit pos.0");
    TGraphErrors *SIGMA_H_B_1 = new TGraphErrors(16,x,sigma_H_1_B,err_x,err_sigma_H_1_B);SIGMA_H_B_1->SetName("SIGMA_H_B_1"); SIGMA_H_B_1->SetTitle("#delta t_{high} pos.1 - fit pos.1");
    TGraphErrors *SIGMA_T_B_0 = new TGraphErrors(16,x,sigma_T_0_B,err_x,err_sigma_T_0_B);SIGMA_T_B_0->SetName("SIGMA_T_B_0"); SIGMA_T_B_0->SetTitle("#delta t_{third} pos.0 - fit pos.0");
    TGraphErrors *SIGMA_T_B_1 = new TGraphErrors(16,x,sigma_T_1_B,err_x,err_sigma_T_1_B);SIGMA_T_B_1->SetName("SIGMA_T_B_1"); SIGMA_T_B_1->SetTitle("#delta t_{third} pos.1 - fit pos.1");
    
    TGraphErrors *SIGMA_L_A_0 = new TGraphErrors(16,x,sigma_L_0_A,err_x,err_sigma_L_0_A);SIGMA_L_A_0->SetName("SIGMA_L_A_0"); SIGMA_L_A_0->SetTitle("#delta t_{low} pos.0 - fit pos.0 #oplus 1");
    TGraphErrors *SIGMA_L_A_1 = new TGraphErrors(16,x,sigma_L_1_A,err_x,err_sigma_L_1_A);SIGMA_L_A_1->SetName("SIGMA_L_A_1"); SIGMA_L_A_1->SetTitle("#delta t_{low} pos.1 - fit pos.0 #oplus 1");
    TGraphErrors *SIGMA_H_A_0 = new TGraphErrors(16,x,sigma_H_0_A,err_x,err_sigma_H_0_A);SIGMA_H_A_0->SetName("SIGMA_H_A_0"); SIGMA_H_A_0->SetTitle("#delta t_{high} pos.0 - fit pos.0 #oplus 1");
    TGraphErrors *SIGMA_H_A_1 = new TGraphErrors(16,x,sigma_H_1_A,err_x,err_sigma_H_1_A);SIGMA_H_A_1->SetName("SIGMA_H_A_1"); SIGMA_H_A_1->SetTitle("#delta t_{high} pos.1 - fit pos.0 #oplus 1");
    TGraphErrors *SIGMA_T_A_0 = new TGraphErrors(16,x,sigma_T_0_A,err_x,err_sigma_T_0_A);SIGMA_T_A_0->SetName("SIGMA_T_A_0"); SIGMA_T_A_0->SetTitle("#delta t_{third} pos.0 - fit pos.0 #oplus 1");
    TGraphErrors *SIGMA_T_A_1 = new TGraphErrors(16,x,sigma_T_1_A,err_x,err_sigma_T_1_A);SIGMA_T_A_1->SetName("SIGMA_T_A_1"); SIGMA_T_A_1->SetTitle("#delta t_{third} pos.1 - fit pos.0 #oplus 1");
    
    
    TGraphErrors *FRAC_L_B_0 = new TGraphErrors(16,x,frac_L_B_0,err_x,err_frac_L_B_0);FRAC_L_B_0->SetName("FRAC_L_B_0"); FRAC_L_B_0->SetTitle("f_{L} pos.0 - fit pos 0");
    TGraphErrors *FRAC_H_B_0 = new TGraphErrors(16,x,frac_H_B_0,err_x,err_frac_H_B_0);FRAC_H_B_0->SetName("FRAC_H_B_0"); FRAC_H_B_0->SetTitle("f_{H} pos.0 - fit pos 0");
    TGraphErrors *FRAC_L_B_1 = new TGraphErrors(16,x,frac_L_B_1,err_x,err_frac_L_B_1);FRAC_L_B_1->SetName("FRAC_L_B_1"); FRAC_L_B_1->SetTitle("f_{L} pos.1 - fit pos 1");
    TGraphErrors *FRAC_H_B_1 = new TGraphErrors(16,x,frac_H_B_1,err_x,err_frac_H_B_1);FRAC_H_B_1->SetName("FRAC_H_B_1"); FRAC_H_B_1->SetTitle("f_{H} pos.1 - fit pos 1");
    TGraphErrors *FRAC_L_A_0 = new TGraphErrors(16,x,frac_L_A_0,err_x,err_frac_L_A_0);FRAC_L_A_0->SetName("FRAC_L_A_0"); FRAC_L_A_0->SetTitle("f_{L} pos.0 - fit pos 0 #oplus 1");
    TGraphErrors *FRAC_H_A_0 = new TGraphErrors(16,x,frac_H_A_0,err_x,err_frac_H_A_0);FRAC_H_A_0->SetName("FRAC_H_A_0"); FRAC_H_A_0->SetTitle("f_{H} pos.0 - fit pos 0 #oplus 1");
    TGraphErrors *FRAC_L_A_1 = new TGraphErrors(16,x,frac_L_A_1,err_x,err_frac_L_A_1);FRAC_L_A_1->SetName("FRAC_L_A_1"); FRAC_L_A_1->SetTitle("f_{L} pos.1 - fit pos 0 #oplus 1");
    TGraphErrors *FRAC_H_A_1 = new TGraphErrors(16,x,frac_H_A_1,err_x,err_frac_H_A_1);FRAC_H_A_1->SetName("FRAC_H_A_1"); FRAC_H_A_1->SetTitle("f_{H} pos.1 - fit pos 0 #oplus 1");
    
    
    TGraphErrors *FRAC_0 = new TGraphErrors(16,x,frac_0,err_x,err_frac_0);FRAC_0->SetName("FRAC_0"); FRAC_0->SetTitle("F_{0}fit pos 0 #oplus 1");
    
    TGraphErrors *REF_sigma_L = new TGraphErrors(16,x,ref_sigma_L,err_x,err_x);REF_sigma_L->SetName("REF_sigma_L"); REF_sigma_L->SetTitle("#delta t_{low}^{ref} pos.0 - fit pos.0");
    TGraphErrors *REF_sigma_H = new TGraphErrors(16,x,ref_sigma_H,err_x,err_x);REF_sigma_H->SetName("REF_sigma_H"); REF_sigma_H->SetTitle("#delta t_{high}^{ref} pos.0 - fit pos.0");
    TGraphErrors *REF_sigma_T = new TGraphErrors(16,x,ref_sigma_T,err_x,err_x);REF_sigma_T->SetName("REF_sigma_T"); REF_sigma_T->SetTitle("#delta t_{third}^{ref} pos.0 - fit pos.0");

    TGraphErrors *SIGMA_FRAC_L_0 = new TGraphErrors(16,frac_L_B_0,sigma_L_0_B,err_frac_L_B_0,err_sigma_L_0_B);SIGMA_FRAC_L_0->SetName("SIGMA_FRAC_L_0"); SIGMA_FRAC_L_0->SetTitle("#delta t [ns] vs f_{L} - pos.0");
    TGraphErrors *SIGMA_FRAC_H_0 = new TGraphErrors(16,frac_L_B_0,sigma_H_0_B,err_frac_L_B_0,err_sigma_H_0_B);SIGMA_FRAC_H_0->SetName("SIGMA_FRAC_H_0"); SIGMA_FRAC_H_0->SetTitle("#delta t [ns] vs f_{L} - pos.0");
    TGraphErrors *SIGMA_FRAC_L_1 = new TGraphErrors(16,frac_L_B_1,sigma_L_1_B,err_frac_L_B_1,err_sigma_L_1_B);SIGMA_FRAC_L_1->SetName("SIGMA_FRAC_L_1"); SIGMA_FRAC_L_1->SetTitle("#delta t [ns] vs f_{L} - pos.1");
    TGraphErrors *SIGMA_FRAC_H_1 = new TGraphErrors(16,frac_L_B_1,sigma_H_1_B,err_frac_L_B_1,err_sigma_H_1_B);SIGMA_FRAC_H_1->SetName("SIGMA_FRAC_H_1"); SIGMA_FRAC_H_1->SetTitle("#delta t [ns] vs f_{L} - pos.1");
    TGraphErrors *SIGMA_FRAC_L_0_A = new TGraphErrors(16,frac_L_B_0,sigma_L_0_A,err_frac_L_B_0,err_sigma_L_0_A);SIGMA_FRAC_L_0_A->SetName("SIGMA_FRAC_L_0_A"); SIGMA_FRAC_L_0_A->SetTitle("#delta t [ns] vs f_{L} - pos.0");
    TGraphErrors *SIGMA_FRAC_H_0_A = new TGraphErrors(16,frac_L_B_0,sigma_H_0_A,err_frac_L_B_0,err_sigma_H_0_A);SIGMA_FRAC_H_0_A->SetName("SIGMA_FRAC_H_0_A"); SIGMA_FRAC_H_0_A->SetTitle("#delta t [ns] vs f_{L} - pos.0");
    TGraphErrors *SIGMA_FRAC_L_1_A = new TGraphErrors(16,frac_L_B_1,sigma_L_1_A,err_frac_L_B_1,err_sigma_L_1_A);SIGMA_FRAC_L_1_A->SetName("SIGMA_FRAC_L_1_A"); SIGMA_FRAC_L_1_A->SetTitle("#delta t [ns] vs f_{L} - pos.1");
    TGraphErrors *SIGMA_FRAC_H_1_A = new TGraphErrors(16,frac_L_B_1,sigma_H_1_A,err_frac_L_B_1,err_sigma_H_1_A);SIGMA_FRAC_H_1_A->SetName("SIGMA_FRAC_H_1_A"); SIGMA_FRAC_H_1_A->SetTitle("#delta t [ns] vs f_{L} - pos.1");
    SIGMA_FRAC_H_0->SetMarkerStyle(20);SIGMA_FRAC_H_0->SetMarkerColor(4);SIGMA_FRAC_H_0->SetMarkerSize(2.0);
    SIGMA_FRAC_L_0->SetMarkerStyle(20);SIGMA_FRAC_L_0->SetMarkerColor(4);SIGMA_FRAC_L_0->SetMarkerSize(2.0);
    SIGMA_FRAC_H_1->SetMarkerStyle(20);SIGMA_FRAC_H_1->SetMarkerColor(4);SIGMA_FRAC_H_1->SetMarkerSize(2.0);
    SIGMA_FRAC_L_1->SetMarkerStyle(20);SIGMA_FRAC_L_1->SetMarkerColor(4);SIGMA_FRAC_L_1->SetMarkerSize(2.0);
    SIGMA_FRAC_H_0_A->SetMarkerStyle(20);SIGMA_FRAC_H_0_A->SetMarkerColor(2);SIGMA_FRAC_H_0_A->SetMarkerSize(2.0);
    SIGMA_FRAC_L_0_A->SetMarkerStyle(20);SIGMA_FRAC_L_0_A->SetMarkerColor(2);SIGMA_FRAC_L_0_A->SetMarkerSize(2.0);
    SIGMA_FRAC_H_1_A->SetMarkerStyle(20);SIGMA_FRAC_H_1_A->SetMarkerColor(2);SIGMA_FRAC_H_1_A->SetMarkerSize(2.0);
    SIGMA_FRAC_L_1_A->SetMarkerStyle(20);SIGMA_FRAC_L_1_A->SetMarkerColor(2);SIGMA_FRAC_L_1_A->SetMarkerSize(2.0);


    FRAC_0->SetMarkerStyle(20); FRAC_0->SetMarkerColor(4);
    DELTA_means_ref->SetMarkerStyle(29);
    DELTA_means_ref->SetMarkerColor(kOrange-2);DELTA_means_ref->SetMarkerSize(2.0);
    DELTA_means_pos0_MC_B2->SetMarkerStyle(29);DELTA_means_pos0_MC_B2->SetMarkerColor(8);DELTA_means_pos0_MC_B2->SetMarkerSize(2.0);
    DELTA_means_pos1_MC_B2->SetMarkerStyle(29);DELTA_means_pos1_MC_B2->SetMarkerColor(8);DELTA_means_pos1_MC_B2->SetMarkerSize(2.0);
    DELTA_means_pos0_MC_PD->SetMarkerStyle(29);DELTA_means_pos0_MC_PD->SetMarkerColor(2);DELTA_means_pos0_MC_PD->SetMarkerSize(2.0);
    DELTA_means_pos0_3G->SetMarkerStyle(29);DELTA_means_pos0_3G->SetMarkerColor(kOrange-2);DELTA_means_pos0_3G->SetMarkerSize(2.0);
    DELTA_means_pos1_3G->SetMarkerStyle(29);DELTA_means_pos1_3G->SetMarkerColor(kOrange-2);DELTA_means_pos1_3G->SetMarkerSize(2.0);


    
    FRAC_SIG->SetMarkerColor(kOrange-2);
    FRAC_SIG->SetMarkerSize(2.0);
    FRAC_SP->SetMarkerColor(kOrange-2);
    FRAC_SP->SetMarkerSize(2.0);
    FRAC_SIG->SetMarkerColor(4);
    FRAC_SP->SetMarkerColor(kOrange-2);
    FRAC_SIG->SetMarkerStyle(20);
    FRAC_SP->SetMarkerStyle(20);
    AMP_SIG->SetMarkerStyle(20);
    AMP_SIG->SetMarkerSize(2.0);
    AMP_SIG->SetMarkerColor(8);

    DELTA_LSP_T0->SetMarkerColor(kOrange-2);
    DELTA_LSP_T0->SetMarkerStyle(20);
    DELTA_LSP_T0->SetMarkerSize(2);
    DELTA_HSP_LSP->SetMarkerColor(kOrange-2);
    DELTA_HSP_LSP->SetMarkerStyle(20);
    DELTA_HSP_LSP->SetMarkerSize(2);
    
    REF_sigma_L->SetMarkerStyle(29);
    REF_sigma_L->SetMarkerColor(kOrange-2);
    REF_sigma_H->SetMarkerStyle(29);
    REF_sigma_H->SetMarkerColor(4);
    REF_sigma_T->SetMarkerStyle(29);
    REF_sigma_T->SetMarkerColor(2);
    
    MEAN_L_B_0->SetMarkerStyle(20);
    MEAN_L_B_1->SetMarkerStyle(21);
    MEAN_L_A_0->SetMarkerStyle(20);
    MEAN_L_A_1->SetMarkerStyle(21);
    MEAN_L_B_0->SetMarkerSize(2);
    MEAN_L_B_1->SetMarkerSize(2);
    MEAN_L_A_0->SetMarkerSize(2);
    MEAN_L_A_1->SetMarkerSize(2);
    DELTA_H_B_0->SetMarkerStyle(22);DELTA_H_B_0->SetMarkerSize(2.0);
    DELTA_H_B_1->SetMarkerStyle(23);DELTA_H_B_1->SetMarkerSize(2.0);
    DELTA_H_A_0->SetMarkerStyle(22);
    DELTA_H_A_1->SetMarkerStyle(23);
    DELTA_T_B_0->SetMarkerStyle(22);
    DELTA_T_B_1->SetMarkerStyle(23);
    DELTA_T_A_0->SetMarkerStyle(22);
    DELTA_T_A_1->SetMarkerStyle(23);
    
    MEAN_L_B_0->SetMarkerColor(4);
    MEAN_L_B_1->SetMarkerColor(4);
    MEAN_L_A_0->SetMarkerColor(2);
    MEAN_L_A_1->SetMarkerColor(2);
    DELTA_H_B_0->SetMarkerColor(4);
    DELTA_H_B_1->SetMarkerColor(4);
    DELTA_H_A_0->SetMarkerColor(2);
    DELTA_H_A_1->SetMarkerColor(2);
    DELTA_T_B_0->SetMarkerColor(8);
    DELTA_T_B_1->SetMarkerColor(8);
    DELTA_T_A_0->SetMarkerColor(4);
    DELTA_T_A_1->SetMarkerColor(4);
    
    SIGMA_L_B_0->SetMarkerStyle(20);
    SIGMA_H_B_0->SetMarkerStyle(21);
    SIGMA_L_A_0->SetMarkerStyle(22);
    SIGMA_H_A_0->SetMarkerStyle(23);
    SIGMA_T_B_0->SetMarkerStyle(21);
    SIGMA_T_A_0->SetMarkerStyle(22);
    
    SIGMA_L_B_1->SetMarkerStyle(20);
    SIGMA_H_B_1->SetMarkerStyle(21);
    SIGMA_L_A_1->SetMarkerStyle(22);
    SIGMA_H_A_1->SetMarkerStyle(23);
    SIGMA_T_B_1->SetMarkerStyle(23);
    SIGMA_T_A_1->SetMarkerStyle(20);
    
    SIGMA_L_B_0->SetMarkerSize(2);//(20);
    SIGMA_H_B_0->SetMarkerSize(2);//(21);
    SIGMA_L_A_0->SetMarkerSize(2);//(22);
    SIGMA_H_A_0->SetMarkerSize(2);//(23);
    SIGMA_T_B_0->SetMarkerSize(2);//(21);
    SIGMA_T_A_0->SetMarkerSize(2);//(22);
    
    SIGMA_L_B_1->SetMarkerSize(2);//(20);
    SIGMA_H_B_1->SetMarkerSize(2);//(21);
    SIGMA_L_A_1->SetMarkerSize(2);//(22);
    SIGMA_H_A_1->SetMarkerSize(2);//(23);
    SIGMA_T_B_1->SetMarkerSize(2);//(23);
    SIGMA_T_A_1->SetMarkerSize(2);//(20);
    
    SIGMA_L_B_0->SetMarkerColor(1);
    SIGMA_H_B_0->SetMarkerColor(8);
    SIGMA_L_A_0->SetMarkerColor(4);
    SIGMA_H_A_0->SetMarkerColor(2);
    SIGMA_T_B_0->SetMarkerColor(8);
    SIGMA_T_A_0->SetMarkerColor(4);
    
    SIGMA_L_B_1->SetMarkerColor(2);
    SIGMA_H_B_1->SetMarkerColor(4);
    SIGMA_L_A_1->SetMarkerColor(8);
    SIGMA_H_A_1->SetMarkerColor(1);
    SIGMA_T_B_1->SetMarkerColor(kOrange-2);
    SIGMA_T_A_1->SetMarkerColor(2);
    
    FRAC_L_B_0->SetMarkerStyle(20);
    FRAC_H_B_0->SetMarkerStyle(21);
    FRAC_L_A_0->SetMarkerStyle(22);
    FRAC_H_A_0->SetMarkerStyle(23);
    FRAC_L_B_0->SetMarkerColor(1);
    FRAC_H_B_0->SetMarkerColor(2);
    FRAC_L_A_0->SetMarkerColor(8);
    FRAC_H_A_0->SetMarkerColor(4);
    FRAC_L_B_1->SetMarkerStyle(20);
    FRAC_H_B_1->SetMarkerStyle(21);
    FRAC_L_A_1->SetMarkerStyle(22);
    FRAC_H_A_1->SetMarkerStyle(23);
    FRAC_L_B_1->SetMarkerColor(1);
    FRAC_H_B_1->SetMarkerColor(2);
    FRAC_L_A_1->SetMarkerColor(8);
    FRAC_H_A_1->SetMarkerColor(4);
    
    
    FRAC_L_B_1->SetMarkerSize(2);//(20);
    FRAC_H_B_1->SetMarkerSize(2);//(21);
    FRAC_L_A_1->SetMarkerSize(2);//(22);
    FRAC_H_A_1->SetMarkerSize(2);//(23);
    FRAC_L_B_0->SetMarkerSize(2);//(20);
    FRAC_H_B_0->SetMarkerSize(2);//(21);
    FRAC_L_A_0->SetMarkerSize(2);//(22);
    FRAC_H_A_0->SetMarkerSize(2);//(23);
    
    TMultiGraph *mg_MEANS_L_0 = new TMultiGraph();mg_MEANS_L_0->SetName("mg_MEANS_L_0");
    mg_MEANS_L_0->SetTitle("T_{0} vs Channel - pos.0");
    mg_MEANS_L_0->Add(MEAN_L_B_0);
    //mg_MEANS_L_0->Add(MEAN_L_A_0);
    
    TCanvas *c_absolute_positions_01_b = new TCanvas("c_absolute_position_01_b","c_absolute_position_01_b",0,0,1124,700);
    //c_absolute_positions_01->Divide(2,1);
    //c_absolute_positions_01->cd(1);
    mg_MEANS_L_0->Draw("AP");
    mg_MEANS_L_0->GetXaxis()->SetTitle("Channel");
    mg_MEANS_L_0->GetYaxis()->SetTitle("T_{0} [ns]");
    mg_MEANS_L_0->GetYaxis()->SetRangeUser(22.75,23.2);
    gPad->BuildLegend();
    gPad->Update();
    gPad->SetGridy();
  
    
    TMultiGraph *mg_DELTA_B_pos0 = new TMultiGraph();mg_DELTA_B_pos0->SetName("mg_DELTA_B_pos0");
    mg_DELTA_B_pos0->SetTitle("#Delta t (#equiv t_{H}-t_{L}) vs Channel - pos.0");
    //mg_DELTA_B->Add(DELTA_means_ref);
    mg_DELTA_B_pos0->Add(DELTA_means_pos0_3G);
    mg_DELTA_B_pos0->Add(DELTA_H_B_0);
    //mg_DELTA_B_pos0->Add(DELTA_means_pos0_MC_B2);
    //mg_DELTA_B_pos0->Add(DELTA_means_pos0_MC_PD);
    mg_DELTA_B_pos0->SetMinimum(0.1);
    mg_DELTA_B_pos0->SetMaximum(0.4);
    TCanvas *c_pos0_deltas_data_MC = new TCanvas("c_pos0_deltas_data_MC","c_pos0_deltas_data_MC",0,0,1124,700);
    mg_DELTA_B_pos0->Draw("AP");
    gPad->BuildLegend();
    gPad->Update();
    gPad->SetGridy();
    mg_DELTA_B_pos0->GetXaxis()->SetTitle("Channel");
    mg_DELTA_B_pos0->GetYaxis()->SetTitle("#Delta t [ns]");
    mg_DELTA_B_pos0->GetYaxis()->SetNdivisions(20);
    gPad->Update();
    
    
    
    
    TMultiGraph *mg_SIGMA_B_L = new TMultiGraph();mg_SIGMA_B_L->SetName("mg_SIGMA_B_L");
    mg_SIGMA_B_L->SetTitle("Time resolution (#delta t) low time vs Channel");
    mg_SIGMA_B_L->Add(REF_sigma_L);
    mg_SIGMA_B_L->Add(SIGMA_L_B_0);
    //mg_SIGMA_B_L->Add(SIGMA_L_A_0);
    TMultiGraph *mg_SIGMA_B_H = new TMultiGraph();mg_SIGMA_B_H->SetName("mg_SIGMA_B_H");
    mg_SIGMA_B_H->SetTitle("Time resolution (#delta t) high time vs Channel");
    mg_SIGMA_B_H->Add(REF_sigma_H);
    mg_SIGMA_B_H->Add(SIGMA_H_B_0);
    //mg_SIGMA_B_H->Add(SIGMA_H_A_0);
    TMultiGraph *mg_SIGMA_B_T = new TMultiGraph();mg_SIGMA_B_T->SetName("mg_SIGMA_B_T");
    mg_SIGMA_B_T->SetTitle("Time resolution (#delta t) third time vs Channel");
    mg_SIGMA_B_T->Add(REF_sigma_T);
    //mg_SIGMA_B_T->Add(SIGMA_T_B_0);
    //mg_SIGMA_B_T->Add(SIGMA_T_A_0);
    mg_SIGMA_B_T->Add(SIGMA_T_B_0);
    //mg_SIGMA_B_T->Add(SIGMA_T_A_1);
    
    
    
    
    ///////NEW MG_GRAPHS/////
    TMultiGraph *mg_DELTA_0 = new TMultiGraph();mg_DELTA_0->SetName("mg_DELTA_0");
    mg_DELTA_0->SetTitle("#Delta t (#equiv t_{H}-t_{L}) vs Channel - pos.0 ");
    //mg_DELTA_0->Add(DELTA_means_ref);
    mg_DELTA_0->Add(DELTA_H_B_0);
    //mg_DELTA_0->Add(DELTA_H_A_0);
    //mg_DELTA_0->Add(DELTA_T_B_0);
    //mg_DELTA_0->Add(DELTA_T_A_0);

    
    
    TMultiGraph *mg_FRAC_0 = new TMultiGraph();mg_FRAC_0->SetName("mg_FRAC_0");
    mg_FRAC_0->SetTitle("Fraction of paths  vs Channel - pos.0 ");
    //mg_FRAC_0->Add(FRAC_H_B_0);
    //mg_FRAC_0->Add(FRAC_H_A_0);
    mg_FRAC_0->Add(FRAC_L_B_0);
    //mg_FRAC_0->Add(FRAC_L_A_0);

    
    TMultiGraph *mg_SIGMA_L_0 = new TMultiGraph();mg_SIGMA_L_0->SetName("mg_SIGMA_L_0");
    mg_SIGMA_L_0->SetTitle("Time resolution (#delta t) low time vs Channel - pos. 0");
    //mg_SIGMA_L_0->Add(REF_sigma_L);
    mg_SIGMA_L_0->Add(SIGMA_L_B_0);
    //mg_SIGMA_L_0->Add(SIGMA_L_A_0);
    TMultiGraph *mg_SIGMA_H_0 = new TMultiGraph();mg_SIGMA_H_0->SetName("mg_SIGMA_H_0");
    mg_SIGMA_H_0->SetTitle("Time resolution (#delta t) high time vs Channel - pos. 0");
    //mg_SIGMA_H_0->Add(REF_sigma_H);
    mg_SIGMA_H_0->Add(SIGMA_H_B_0);
    //mg_SIGMA_H_0->Add(SIGMA_H_A_0);
    TMultiGraph *mg_SIGMA_T_0 = new TMultiGraph();mg_SIGMA_T_0->SetName("mg_SIGMA_T_0");
    mg_SIGMA_T_0->SetTitle("Time resolution (#delta t) third time vs Channel - pos. 0");
    //mg_SIGMA_T_0->Add(REF_sigma_T);
    mg_SIGMA_T_0->Add(SIGMA_T_B_0);
    //mg_SIGMA_T_0->Add(SIGMA_T_A_0);
    
    
    TMultiGraph *mg_FRAC_SIGMA_L_0 = new TMultiGraph();mg_FRAC_SIGMA_L_0->SetName("mg_FRAC_SIGMA_L_0");
    mg_FRAC_SIGMA_L_0->SetTitle("#delta t [ns] vs f_{L} - pos.0 - low time peak");
    mg_FRAC_SIGMA_L_0->Add(SIGMA_FRAC_L_0);
    mg_FRAC_SIGMA_L_0->Add(SIGMA_FRAC_L_0_A);
    TMultiGraph *mg_FRAC_SIGMA_H_0 = new TMultiGraph();mg_FRAC_SIGMA_H_0->SetName("mg_FRAC_SIGMA_H_0");
    mg_FRAC_SIGMA_H_0->SetTitle("#delta t [ns] vs f_{L} - pos.0 - high time peak");
    mg_FRAC_SIGMA_H_0->Add(SIGMA_FRAC_H_0);
    mg_FRAC_SIGMA_H_0->Add(SIGMA_FRAC_H_0_A);
 


    TCanvas *c_frac_sig = new TCanvas("c_frac_sigma","c_frac_sigma",0,0,1124,700);
    c_frac_sig->Divide(1,2);
    c_frac_sig->cd(1);
    mg_FRAC_SIGMA_L_0->Draw("AP");
    gPad->Update();
    mg_FRAC_SIGMA_L_0->GetXaxis()->SetTitle("f_{L}");
    mg_FRAC_SIGMA_L_0->GetYaxis()->SetTitle("#delta t [ns]");
    c_frac_sig->cd(2);
    mg_FRAC_SIGMA_H_0->Draw("AP");
    gPad->Update();
    mg_FRAC_SIGMA_H_0->GetXaxis()->SetTitle("f_{L}");
    mg_FRAC_SIGMA_H_0->GetYaxis()->SetTitle("#delta t [ns]");
  
    //////
    
    TCanvas *c_DELTA_FRAC = new TCanvas("mean peak separation","mean peak separation",0,0,1124,700);
    //c_DELTA_FRAC->Divide(1,2);
    c_DELTA_FRAC->cd();
    mg_DELTA_0->Draw("AP");
    gPad->Update();
    mg_DELTA_0->GetXaxis()->SetTitle("Channel");
    mg_DELTA_0->GetYaxis()->SetTitle("#Delta t[ns]");
    gPad->Update();
    mg_DELTA_0->Draw("AP");
    gPad->BuildLegend();
    
    
    TCanvas *c_FRAC = new TCanvas("c_FRAC","c_FRAC",0,0,1124,700);
    //   c_FRAC->Divide(1,2);
    c_FRAC->cd();
    mg_FRAC_0->Draw("AP");
    gPad->Update();
    mg_FRAC_0->GetXaxis()->SetTitle("Channel");
    mg_FRAC_0->GetYaxis()->SetTitle("#Delta t[ns]");
    mg_FRAC_0->GetYaxis()->SetRangeUser(0,0.9);
    gPad->SetGridy();
    gPad->Update();
    mg_FRAC_0->Draw("AP");
    gPad->BuildLegend();
  
    /*
      TCanvas *c1 = new TCanvas("c1","c1");
      c1->Divide(2,2);
      c1->cd(1);
      mg_DELTA_B->Draw("AP");
      gPad->Update();
      mg_DELTA_B->GetXaxis()->SetTitle("Channel");
      mg_DELTA_B->GetYaxis()->SetTitle("#Delta t[ns]");
      gPad->Update();
      mg_DELTA_B->Draw("AP");
      gPad->BuildLegend();
      c1->cd(2);
      FRAC_0->Draw("AP");
      mg_SIGMA_B_L->Draw("AP");
      gPad->Update();
      mg_SIGMA_B_L->GetXaxis()->SetTitle("Channel");
      mg_SIGMA_B_L->GetYaxis()->SetTitle("#delta t [ns]");
      gPad->Update();
      mg_SIGMA_B_L->Draw("AP");
      gPad->BuildLegend();
      c1->cd(3);
      mg_SIGMA_B_H->Draw("AP");
      gPad->Update();
      mg_SIGMA_B_H->GetXaxis()->SetTitle("Channel");
      mg_SIGMA_B_H->GetYaxis()->SetTitle("#delta t [ns]");
      gPad->Update();
      mg_SIGMA_B_H->Draw("AP");
      gPad->BuildLegend();
      
      c1->cd(4);
      mg_SIGMA_B_T->Draw("AP");
      gPad->Update();
      mg_SIGMA_B_T->GetXaxis()->SetTitle("Channel");
      mg_SIGMA_B_T->GetYaxis()->SetTitle("#delta t [ns]");
      gPad->Update();
      mg_SIGMA_B_T->Draw("AP");
      gPad->BuildLegend();
    */
    TCanvas* c_Res_pos0 = new TCanvas("c_Res_pos0","c_Res_pos0",0,0,1124,700);
    c_Res_pos0->Divide(2,1);
    c_Res_pos0->cd(1);
    mg_SIGMA_L_0->Draw("AP");
    line->Draw("same");
    gPad->Update();
    mg_SIGMA_L_0->GetXaxis()->SetTitle("Channel");
    mg_SIGMA_L_0->GetYaxis()->SetTitle("#delta t [ns]");
    gPad->Update();
    gPad->BuildLegend();
    c_Res_pos0->cd(2);
    mg_SIGMA_H_0->Draw("AP");
    line->Draw("same");
    gPad->Update();
    mg_SIGMA_H_0->GetXaxis()->SetTitle("Channel");
    mg_SIGMA_H_0->GetYaxis()->SetTitle("#delta t [ns]");
    gPad->Update();
    gPad->BuildLegend();
    /*
      c_Res_pos0->cd(3);
      mg_SIGMA_T_0->Draw("AP");
      line->Draw("same");
      gPad->BuildLegend();
      gPad->Update();
      mg_SIGMA_T_0->GetXaxis()->SetTitle("Channel");
      mg_SIGMA_T_0->GetYaxis()->SetTitle("#delta t [ns]");
      gPad->Update();
    */
  
    /*
      c_Res_pos1->cd(3);
      mg_SIGMA_T_1->Draw("AP");
      line->Draw("same");
      gPad->BuildLegend();
      gPad->Update();
      mg_SIGMA_T_1->GetXaxis()->SetTitle("Channel");
      mg_SIGMA_T_1->GetYaxis()->SetTitle("#delta t [ns]");
      gPad->Update();
    */
    /*
      int X[16];
      int Fit_status_0[16];
      int Fit_status_1[16];
      int Fit_status[16];
      for(int g=0;g<16;g++){
      Fit_status_0[g]=Fits_status[0+g];
      Fit_status_1[g]=Fits_status[1+g];
      Fit_status[g]  =Fits_status[2+g];
      X[g]=g;
      cout<<Fit_status_0[g]<<"  "<<Fit_status_1[g]<<"  "<<Fit_status[g]<<"  "<<endl;
      }
      
      TGraph *GR_Fit_status_0 = new TGraph(16,X,Fit_status_0);GR_Fit_status_0->SetName("GR_Fit_status_0"); GR_Fit_status_0->SetTitle("Fits status - fit pos 0");
      TGraph *GR_Fit_status_1 = new TGraph(16,X,Fit_status_1);GR_Fit_status_1->SetTitle("Fits status - fit pos 1");
      TGraph *GR_Fit_status = new TGraph(16,X,Fit_status);GR_Fit_status->SetTitle("Fits status - fit pos 0 #oplus 1");
      GR_Fit_status_0->SetMarkerStyle(20);
      GR_Fit_status_0->SetMarkerColor(1);
      GR_Fit_status_1->SetMarkerStyle(21);
      GR_Fit_status_1->SetMarkerColor(2);
      GR_Fit_status->SetMarkerStyle(22);
      GR_Fit_status->SetMarkerColor(3);
      TMultiGraph *mg_FITSTATUS = new TMultiGraph();mg_FITSTATUS->SetName("mg_FITSTATUS");
      mg_FITSTATUS->SetTitle("Fits status");
      mg_FITSTATUS->Add(GR_Fit_status_0);
      mg_FITSTATUS->Add(GR_Fit_status_1);
      mg_FITSTATUS->Add(GR_Fit_status);
    */
    
    c_pos0_AllChannels->Draw();    c_pos0_AllChannels_pulls->Draw();  c_pos0_AllChannels_corr->Draw(); 
    
    TH2D* h2_DTH_0 = new TH2D("h2_DTH_0","h2_DTH_0",4,0,4,4,0,4);
    h2_DTH_0->Reset();
    h2_DTH_0->SetBinContent(1,1,Delta_H_0_B[12]);
    h2_DTH_0->SetBinContent(1,2,Delta_H_0_B[13]);
    h2_DTH_0->SetBinContent(1,3,Delta_H_0_B[14]);
    h2_DTH_0->SetBinContent(1,4,Delta_H_0_B[15]);
    h2_DTH_0->SetBinContent(2,1,Delta_H_0_B[8]);
    h2_DTH_0->SetBinContent(2,2,Delta_H_0_B[9]);
    h2_DTH_0->SetBinContent(2,3,Delta_H_0_B[10]);
    h2_DTH_0->SetBinContent(2,4,Delta_H_0_B[11]);
    h2_DTH_0->SetBinContent(3,1,Delta_H_0_B[4]);
    h2_DTH_0->SetBinContent(3,2,Delta_H_0_B[5]);
    h2_DTH_0->SetBinContent(3,3,Delta_H_0_B[6]);
    h2_DTH_0->SetBinContent(3,4,Delta_H_0_B[7]);
    h2_DTH_0->SetBinContent(4,1,Delta_H_0_B[0]);
    h2_DTH_0->SetBinContent(4,2,Delta_H_0_B[1]);
    h2_DTH_0->SetBinContent(4,3,Delta_H_0_B[2]);
    h2_DTH_0->SetBinContent(4,4,Delta_H_0_B[3]);
    h2_DTH_0->SetTitle("#Delta T - pos. 0");
    
    
    TH2D* h2_FracL_0 = new TH2D("h2_FracL_0","h2_FracL_0",4,0,4,4,0,4);
    h2_FracL_0->Reset();
    h2_FracL_0->SetBinContent(1,1,frac_L_B_0[12]);
    h2_FracL_0->SetBinContent(1,2,frac_L_B_0[13]);
    h2_FracL_0->SetBinContent(1,3,frac_L_B_0[14]);
    h2_FracL_0->SetBinContent(1,4,frac_L_B_0[15]);
    h2_FracL_0->SetBinContent(2,1,frac_L_B_0[8]);
    h2_FracL_0->SetBinContent(2,2,frac_L_B_0[9]);
    h2_FracL_0->SetBinContent(2,3,frac_L_B_0[10]);
    h2_FracL_0->SetBinContent(2,4,frac_L_B_0[11]);
    h2_FracL_0->SetBinContent(3,1,frac_L_B_0[4]);
    h2_FracL_0->SetBinContent(3,2,frac_L_B_0[5]);
    h2_FracL_0->SetBinContent(3,3,frac_L_B_0[6]);
    h2_FracL_0->SetBinContent(3,4,frac_L_B_0[7]);
    h2_FracL_0->SetBinContent(4,1,frac_L_B_0[0]);
    h2_FracL_0->SetBinContent(4,2,frac_L_B_0[1]);
    h2_FracL_0->SetBinContent(4,3,frac_L_B_0[2]);
    h2_FracL_0->SetBinContent(4,4,frac_L_B_0[3]);
    h2_FracL_0->SetTitle("f_{L} - pos. 0");
    
    TH2D* h2_DTH_1 = new TH2D("h2_DTH_1","h2_DTH_1",4,0,4,4,0,4);
    h2_DTH_1->Reset();
    h2_DTH_1->SetBinContent(1,1,Delta_H_1_B[12]);
    h2_DTH_1->SetBinContent(1,2,Delta_H_1_B[13]);
    h2_DTH_1->SetBinContent(1,3,Delta_H_1_B[14]);
    h2_DTH_1->SetBinContent(1,4,Delta_H_1_B[15]);
    h2_DTH_1->SetBinContent(2,1,Delta_H_1_B[8]);
    h2_DTH_1->SetBinContent(2,2,Delta_H_1_B[9]);
    h2_DTH_1->SetBinContent(2,3,Delta_H_1_B[10]);
    h2_DTH_1->SetBinContent(2,4,Delta_H_1_B[11]);
    h2_DTH_1->SetBinContent(3,1,Delta_H_1_B[4]);
    h2_DTH_1->SetBinContent(3,2,Delta_H_1_B[5]);
    h2_DTH_1->SetBinContent(3,3,Delta_H_1_B[6]);
    h2_DTH_1->SetBinContent(3,4,Delta_H_1_B[7]);
    h2_DTH_1->SetBinContent(4,1,Delta_H_1_B[0]);
    h2_DTH_1->SetBinContent(4,2,Delta_H_1_B[1]);
    h2_DTH_1->SetBinContent(4,3,Delta_H_1_B[2]);
    h2_DTH_1->SetBinContent(4,4,Delta_H_1_B[3]);
    h2_DTH_1->SetTitle("#Delta T - pos. 1");
    
    
    TH2D* h2_FracL_1 = new TH2D("h2_FracL_1","h2_FracL_1",4,0,4,4,0,4);
    h2_FracL_1->Reset();
    h2_FracL_1->SetBinContent(1,1,frac_L_B_1[12]);
    h2_FracL_1->SetBinContent(1,2,frac_L_B_1[13]);
    h2_FracL_1->SetBinContent(1,3,frac_L_B_1[14]);
    h2_FracL_1->SetBinContent(1,4,frac_L_B_1[15]);
    h2_FracL_1->SetBinContent(2,1,frac_L_B_1[8]);
    h2_FracL_1->SetBinContent(2,2,frac_L_B_1[9]);
    h2_FracL_1->SetBinContent(2,3,frac_L_B_1[10]);
    h2_FracL_1->SetBinContent(2,4,frac_L_B_1[11]);
    h2_FracL_1->SetBinContent(3,1,frac_L_B_1[4]);
    h2_FracL_1->SetBinContent(3,2,frac_L_B_1[5]);
    h2_FracL_1->SetBinContent(3,3,frac_L_B_1[6]);
    h2_FracL_1->SetBinContent(3,4,frac_L_B_1[7]);
    h2_FracL_1->SetBinContent(4,1,frac_L_B_1[0]);
    h2_FracL_1->SetBinContent(4,2,frac_L_B_1[1]);
    h2_FracL_1->SetBinContent(4,3,frac_L_B_1[2]);
    h2_FracL_1->SetBinContent(4,4,frac_L_B_1[3]);
    h2_FracL_1->SetTitle("f_{L} - pos. 1");
    
    /*
      TCanvas *cNew = new TCanvas("cNew","cNew");
      cNew->Divide(2,2);
      cNew->cd(1);
      h2_DTH_0->Draw("colz:text");
      cNew->cd(3);
      h2_FracL_0->Draw("colz:text");
      cNew->cd(2);
      h2_DTH_1->Draw("colz:text");
      cNew->cd(4);
      h2_FracL_1->Draw("colz:text");
      
      /*
      TCanvas *cD = new TCanvas();
      cD->cd();
      mg_DELTA_0->Draw("AP");
      gPad->Update();
      mg_DELTA_0->GetXaxis()->SetTitle("Channel");
      mg_DELTA_0->GetYaxis()->SetTitle("#Delta t[ns]");
      gPad->Update();
      mg_DELTA_0->Draw("AP");
      gPad->BuildLegend();
    */
  }

  TString out_path="fit_results";
  TFile *f_data = new TFile(out_path+"/"+input_tune+"_FR.root","recreate");
  f_data->cd();
  
  MEAN_L_B_0->Write();
  DELTA_H_B_0->Write();
  SIGMA_L_B_0->Write();
  SIGMA_H_B_0->Write();
  FRAC_SIG->Write();
  AMP_SIG->Write();
  DELTA_LSP_T0->Write();
  DELTA_HSP_LSP->Write();
  
  mg_MEANS_L_0->Write();
  mg_FRAC_0->Write();
  mg_SIGMA_L_0->Write();
  mg_SIGMA_H_0->Write();
  
  c_absolute_positions_01_b->Write();
  c_FRAC->Write();
  c_Res_pos0->Write();
  c_pos0_AllChannels->Write();
  c_pos0_AllChannels_pulls->Write();
  c_pos0_AllChannels_corr->Write();
  c_pos0_AllChannels_amp->Write();
  
  f_data->cd();
  f_data->Close();
  
  return relative_weight_Frac_0;
}





/*
  void draw_raw_data(){
  
  TCanvas *c_0 = new TCanvas("pos.0","pos.0");
  c_0->Divide(4,4);
  TCanvas *c_1 = new TCanvas("pos.1","pos.1");
  c_1->Divide(4,4);
  TCanvas *c_01 = new TCanvas("pos.0 #bigcup 1","pos. 0 #bigcup 1");
  c_01->Divide(4,4);
  
  int index_channel_pixel;
  for(int i=0;i<16;i++){
  if(i==0)  {index_channel_pixel=1 ;}
  if(i==1)  {index_channel_pixel=5 ;}
  if(i==2)  {index_channel_pixel=9 ;}
  if(i==3)  {index_channel_pixel=13 ;}
  if(i==4)  {index_channel_pixel=2 ;}
  if(i==5)  {index_channel_pixel=6 ;}
  if(i==6)  {index_channel_pixel=10 ;}
  if(i==7)  {index_channel_pixel=14 ;}
  if(i==8)  {index_channel_pixel=3 ;}
  if(i==9)  {index_channel_pixel=7 ;}
  if(i==10) {index_channel_pixel=11 ;}
  if(i==11) {index_channel_pixel=15 ;}
  if(i==12) {index_channel_pixel=4 ;}
  if(i==13) {index_channel_pixel=8 ;}
  if(i==14) {index_channel_pixel=12 ;}
  if(i==15) {index_channel_pixel=16 ;}
  TH1D* h_0 = f_input_histogram->Get(Form("fiber0-%d",i));
  TH1D* h_1 = f_input_histogram->Get(Form("fiber1-%d",i));
  TH1D* h_01 = h_0->Clone();
  h_01->Add(h_1);
  c_0->cd(index_channel_pixel); h_0->Draw(); h_0->GetXaxis()->SetRangeUser(8,9.8);h_0->GetXaxis()->SetLabelSize(0.075);h_0->SetTitle("Photon detection time [ns]");
  c_1->cd(index_channel_pixel); h_1->Draw(); h_1->GetXaxis()->SetRangeUser(8,9.8);h_1->GetXaxis()->SetLabelSize(0.075);h_1->SetTitle("Photon detection time [ns]");//h_1->GetXaxis()->SetTitle("t [ns]");
  c_01->cd(index_channel_pixel); h_01->SetLineColor(1); h_01->Draw(); h_0->SetLineColor(2); h_1->SetLineColor(4);  h_0->Draw("same"); h_1->Draw("same");h_01->GetXaxis()->SetRangeUser(8,9.8);h_01->GetXaxis()->SetLabelSize(0.075);;h_01->SetTitle("Photon detection time [ns]");
  }
  
  c_0->cd(); c_0->Draw();
  c_1->cd(); c_1->Draw();
  c_01->cd(); c_01->Draw();
  }
  
  
  
  
*/













