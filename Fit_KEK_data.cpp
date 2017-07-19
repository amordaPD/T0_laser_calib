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

void make_KEK_data_histos(int my_pixelID){
  Float_t upper_bound_hist=2;
  Int_t n_bins = 200;
  TFile *file_input = new TFile("run004881_TBC4855-4858_slot01_digits.root");
  TTree *t_input = (TTree*)file_input->Get("laser");
  TH1D *h_temp = new TH1D("h_temp","h_temp",500,-50,0);
  TCut cut = Form("-20<time&&time<0&&quality==1&&pixel==%i",my_pixelID);
  t_input->Project("h_temp","time",cut);
  //h_temp->Draw();
  Float_t max_bin = h_temp->GetMaximumBin();
  TAxis *xaxis = h_temp->GetXaxis(); 
  Double_t max_pos = xaxis->GetBinCenter(max_bin);
  //cout<< max_pos <<endl;
  TH1D *h_time = new TH1D("h_time","Time [ns]",n_bins,-1,upper_bound_hist);
  t_input->Project("h_time",Form("time-%f",max_pos),cut);
  /*
  Float_t time_sc;
  Int_t pixelID=-99;
  Int_t quality=-9;
  t_input->SetBranchAddress("time",&time_sc);
  t_input->SetBranchAddress("pixel",&pixelID);
  t_input->SetBranchAddress("quality",&quality);

  Int_t n_entries = t_input->GetEntries();
  for(int i=0; i<n_entries; i++){t_input->GetEntry(i); if(-20<time_sc&&time_sc<0&&quality==1&&pixelID==my_pixelID) {h_time->Fill(time_sc-max_pos);}} 
    
  //h_time->Draw("E");
  */

  TFile *file_input_MC = new TFile("ana_laser_s01_0reso_500k.root");
  TTree *tree_MC = (TTree*)file_input_MC->Get("tree_laser");
  TCut cut_MC = Form("propTime<1&&pixelID==%i",my_pixelID);
  
  TH1D *h_MC_tot_temp = new TH1D ("h_MC_tot_temp","h_MC_tot_temp",100,0.3,1);
  tree_MC->Project("h_MC_tot_temp","propTime",cut_MC);
  
  Float_t max_bin_MC = h_MC_tot_temp->GetMaximumBin();
  TAxis *xaxis_MC = h_MC_tot_temp->GetXaxis(); 
  Double_t max_pos_MC = xaxis_MC->GetBinCenter(max_bin_MC);
  TH1D *h_MC_tot = new TH1D ("h_MC_tot","h_MC_tot",n_bins,-1,upper_bound_hist);
  tree_MC->Project("h_MC_tot",Form("propTime-%f",max_pos_MC),cut_MC);
  TH1D *h_MC_f1 = new TH1D ("h_MC_f1","h_MC_f1",n_bins,-1,upper_bound_hist);
  TH1D *h_MC_f2 = new TH1D ("h_MC_f2","h_MC_f2",n_bins,-1,upper_bound_hist);
  TH1D *h_MC_f3 = new TH1D ("h_MC_f3","h_MC_f3",n_bins,-1,upper_bound_hist);
  TH1D *h_MC_f4 = new TH1D ("h_MC_f4","h_MC_f4",n_bins,-1,upper_bound_hist);
  TH1D *h_MC_f5 = new TH1D ("h_MC_f5","h_MC_f5",n_bins,-1,upper_bound_hist);
  TH1D *h_MC_f6 = new TH1D ("h_MC_f6","h_MC_f6",n_bins,-1,upper_bound_hist);  
  TH1D *h_MC_f7 = new TH1D ("h_MC_f7","h_MC_f7",n_bins,-1,upper_bound_hist);
  TH1D *h_MC_f8 = new TH1D ("h_MC_f8","h_MC_f8",n_bins,-1,upper_bound_hist);  
  TH1D *h_MC_f9 = new TH1D ("h_MC_f9","h_MC_f9",n_bins,-1,upper_bound_hist);
  tree_MC->Project("h_MC_f1",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==1&&pixelID==%i",my_pixelID));
  tree_MC->Project("h_MC_f2",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==2&&pixelID==%i",my_pixelID));
  tree_MC->Project("h_MC_f3",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==3&&pixelID==%i",my_pixelID));
  tree_MC->Project("h_MC_f4",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==4&&pixelID==%i",my_pixelID));
  tree_MC->Project("h_MC_f5",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==5&&pixelID==%i",my_pixelID));
  tree_MC->Project("h_MC_f6",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==6&&pixelID==%i",my_pixelID));
  tree_MC->Project("h_MC_f7",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==7&&pixelID==%i",my_pixelID));
  tree_MC->Project("h_MC_f8",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==8&&pixelID==%i",my_pixelID));
  tree_MC->Project("h_MC_f9",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==9&&pixelID==%i",my_pixelID));
  h_MC_tot->SetLineColor(1);
  h_MC_f1->SetLineColor(11);
  h_MC_f2->SetLineColor(2);
  h_MC_f3->SetLineColor(3);
  h_MC_f4->SetLineColor(4);
  h_MC_f5->SetLineColor(13);
  h_MC_f6->SetLineColor(6);
  h_MC_f7->SetLineColor(7);
  h_MC_f8->SetLineColor(8);
  h_MC_f9->SetLineColor(9);
  h_MC_f1->Scale(h_time->GetMaximum()/h_MC_tot->GetMaximum());
  h_MC_f2->Scale(h_time->GetMaximum()/h_MC_tot->GetMaximum());
  h_MC_f3->Scale(h_time->GetMaximum()/h_MC_tot->GetMaximum());
  h_MC_f4->Scale(h_time->GetMaximum()/h_MC_tot->GetMaximum());
  h_MC_f5->Scale(h_time->GetMaximum()/h_MC_tot->GetMaximum());
  h_MC_f6->Scale(h_time->GetMaximum()/h_MC_tot->GetMaximum());
  h_MC_f7->Scale(h_time->GetMaximum()/h_MC_tot->GetMaximum());
  h_MC_f8->Scale(h_time->GetMaximum()/h_MC_tot->GetMaximum());
  h_MC_f9->Scale(h_time->GetMaximum()/h_MC_tot->GetMaximum());
  h_MC_tot->Scale(h_time->GetMaximum()/h_MC_tot->GetMaximum());

  TCanvas *c = new TCanvas("c","c");
  h_time->Draw("E");
  h_MC_tot->Draw("same");
  h_MC_f1->Draw("same");
  h_MC_f2->Draw("same");
  h_MC_f3->Draw("same");
  h_MC_f4->Draw("same");
  h_MC_f5->Draw("same");
  h_MC_f6->Draw("same");
  h_MC_f7->Draw("same");
  h_MC_f8->Draw("same");
  h_MC_f9->Draw("same");
  
  
  
  TFile *f_data = new TFile("KEK_data_histos.root","recreate");
  f_data->cd();
  h_time->Write();
  h_MC_tot->Write();
  h_MC_f1->Write();
  h_MC_f2->Write();
  h_MC_f3->Write();
  h_MC_f4->Write();
  h_MC_f5->Write();
  h_MC_f6->Write();
  h_MC_f7->Write();
  h_MC_f8->Write();
  h_MC_f9->Write();
  c->Write();
  f_data->cd();
  f_data->Close();
  delete f_data;
  delete h_time;
  delete h_MC_tot;
  delete h_MC_f1;
  delete h_MC_f2;
  delete h_MC_f3;
  delete h_MC_f4;
  delete h_MC_f5;
  delete h_MC_f6;
  delete h_MC_f7;
  delete h_MC_f8;
  delete h_MC_f9;
  delete c;
}


RooPlot* Fit_KEK_data(bool CB_model, int column_number, int row_number){
  /////////////////////////////////////////////////////////////////
  ////// HERE Starts the fit part /////////////////////////////////
  /////////////////////////////////////////////////////////////////

  TFile *f_input = new TFile(Form("KEK_data_histos_col_%i.root",column_number));

  Int_t pixelID=column_number+64*(row_number-1);
  
  TH1D* h_time = f_input->Get(Form("histos_%i/h_time",pixelID));
 
  TH1D* h_MC_tot = f_input->Get(Form("histos_%i/h_MC_tot",pixelID));

  TCanvas *can0 = f_input->Get(Form("histos_%i/c",pixelID));
  can0->Draw();


  
  TString _draw_results="draw";

  bool do_prefit=true;
  bool use_NLL=true; //to set the use of fitTo method of RooAbsPdf or the explicit construction of the nll  ///true recomended
  int MN_output_print_level=-1;
  int MN_output_print_level_prefit;
  bool print_prefit_info=false;
  const char * Type_minim="Minuit";
  const char * Algo_minim="minimize";//
  const char * Type_minim_pf="Minuit";//"Minuit2";//
  const char * Algo_minim_pf="minimize";//"scan";//
  bool binned_fit = true; //recomended true if you don't want to wait 7 minutes for the fit output
  bool compute_FWHM = true;
  bool direct_parametrization =true;
  bool fit_for_first_peak=false;
  int  fiber_position_calibration_peak=0;
  bool fit_real_FiberCombs_data=true;
  int bkg_Chebychev_polynomial_degree=1;//set to n to have a n+1 degree Chebychev Polynomial!!!!!!!!!
  bool add_background_component=true;
  int amplitude_cut = -40;
  
  bool suppress_negligible_first_peak=false;
  bool do_simultaneous_fit=false;
  bool add_third_signal_pos0=false;
  
  bool simulate_CB_tail=false;
  
  
  if(!print_prefit_info){MN_output_print_level_prefit=-1;}else{MN_output_print_level_prefit=MN_output_print_level;}
  bool draw_results;
  if(_draw_results=="draw"){draw_results=true;}else if(_draw_results=="blind"){draw_results=false;}else{draw_results=false;}


  
  double my_low_x=-1;
  double my_up_x=2;
  TAxis *xaxis_MC = h_MC_tot->GetXaxis();
  xaxis_MC->SetRange(0,(TMath::Abs(my_low_x)/(my_up_x-my_low_x))*h_MC_tot->GetNbinsX()-1);
  Float_t max_bin_MC = h_MC_tot->GetMaximumBin();
  Double_t max_first_pos_MC = xaxis_MC->GetBinCenter(max_bin_MC);
  cout<<"ciao "<<max_first_pos_MC<<endl;
  
  RooRealVar x("Time","Time [ns]",my_low_x,my_up_x) ;
  RooDataHist ds_0("ds_0","ds_0",RooArgSet(x),Import(*h_time)) ;
  RooDataHist ds_0_H("ds_0_H","ds_0_H",RooArgSet(x),Import(*h_time)) ;
  

  double low_x_0;
  double up_x_0;
  double starting_mean_H_0;
  double low_mean_H_0;
  double up_mean_H_0;
  double starting_mean_L_0;
  double low_mean_L_0;
  double up_mean_L_0;
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

  
  double starting_alpha_CB_L;
  double low_alpha_CB_L;
  double up_alpha_CB_L;
  double starting_n_CB_L;
  double low_n_CB_L;
  double up_n_CB_L;

  double starting_alpha_CB_H;
  double low_alpha_CB_H;
  double up_alpha_CB_H;
  double starting_n_CB_H;
  double low_n_CB_H;
  double up_n_CB_H;

  double starting_alpha_CB_T;
  double low_alpha_CB_T;
  double up_alpha_CB_T;
  double starting_n_CB_T;
  double low_n_CB_T;
  double up_n_CB_T;
  /////POS 0
   low_x_0=my_low_x; up_x_0=my_up_x;
  
  low_delta_H_0=0.180;
   up_delta_H_0=0.5;
  
  low_delta_T_0=0;
   up_delta_T_0=0.5;
  
  low_sigma_L_0=0.035;
   up_sigma_L_0=0.5200;
  
  low_sigma_H_0=0.045;
   up_sigma_H_0=0.5150;
  
  low_sigma_T_0=0.050;
   up_sigma_T_0=0.500;
  
  low_alpha_0=0.01;
   up_alpha_0=0.8;
  
  low_beta_0=0.5;
   up_beta_0=1.0;

  
   starting_mean_L_0=22.6;
   starting_delta_H_0=TMath::Abs(max_first_pos_MC);//0.3;
   starting_delta_T_0=0.040;
   starting_sigma_L_0=0.080;
   starting_sigma_H_0=0.080;
   starting_sigma_T_0=0.1000;
   
   starting_alpha_0=0.25;
   starting_beta_0=0.9;
   
   starting_mean_H_0=0;//starting_mean_L_0+starting_delta_H_0;
   low_mean_L_0=starting_mean_L_0-0.315;
   up_mean_L_0=starting_mean_L_0+0.315;
   low_mean_H_0=starting_mean_H_0-0.05;//315;
   up_mean_H_0=starting_mean_H_0+0.05;//0.315;

 
   starting_alpha_CB_L=0.0;//-0.35;
  starting_n_CB_L=6;

  low_alpha_CB_L=-5;
   up_alpha_CB_L=5;
   
  low_n_CB_L=0;
   up_n_CB_L=200;


  starting_alpha_CB_H=-0.35;
  starting_n_CB_H=6;

  low_alpha_CB_H=-5;
   up_alpha_CB_H=0.0;
   
  low_n_CB_H=0;
   up_n_CB_H=200;

   starting_alpha_CB_T=-0.35;
  starting_n_CB_T=6;

  low_alpha_CB_T=-5;
   up_alpha_CB_T=0.0;
   
  low_n_CB_T=0;
   up_n_CB_T=20;

   /* 
  double starting_alpha_CB;
  double low_alpha_CB;
  double up_alpha_CB;
  double starting_n_CB;
  double low_n_CB;
  double up_n_CB;

  starting_alpha_CB=-0.35;
  starting_n_CB=6;

  low_alpha_CB=-5;
   up_alpha_CB=0.0;
   
  low_n_CB=0;
   up_n_CB=10;
  RooRealVar alpha_CB("alpha_CB","alpha parameter of CB",   starting_alpha_CB,low_alpha_CB,up_alpha_CB);
  RooRealVar     n_CB("n_CB",    "exponential decay of CB ",starting_n_CB,low_n_CB,up_n_CB);
   */
  RooRealVar alpha_CB_L("alpha_CB_L","alpha parameter of CB",   starting_alpha_CB_L,low_alpha_CB_L,up_alpha_CB_L);
  RooRealVar     n_CB_L("n_CB_L",    "exponential decay of CB ",starting_n_CB_L,low_n_CB_L,up_n_CB_L);
  RooRealVar alpha_CB_H("alpha_CB_H","alpha parameter of CB",   starting_alpha_CB_H,low_alpha_CB_H,up_alpha_CB_H);
  RooRealVar     n_CB_H("n_CB_H",    "exponential decay of CB ",starting_n_CB_H,low_n_CB_H,up_n_CB_H);
  RooRealVar alpha_CB_T("alpha_CB_T","alpha parameter of CB",   starting_alpha_CB_T,low_alpha_CB_T,up_alpha_CB_T);
  RooRealVar     n_CB_T("n_CB_T",    "exponential decay of CB ",starting_n_CB_T,low_n_CB_T,up_n_CB_T);
  
  RooRealVar Delta_T_0("Delta_T_0","Delta T pos 0",starting_delta_T_0,low_delta_T_0,up_delta_T_0);
  RooRealVar Delta_H_0("Delta_H_0","Delta H pos 0",starting_delta_H_0,low_delta_H_0,up_delta_H_0);
  RooRealVar mean_H_0("mean_H_0","mean_{H}^{0}",starting_mean_H_0,low_mean_H_0,up_mean_H_0);
  RooFormulaVar mean_L_0("mean_L_0","mean_L_0","mean_H_0-Delta_H_0",RooArgList(mean_H_0,Delta_H_0));
  RooFormulaVar mean_T_0("mean_T_0","mean_T_0","mean_H_0+Delta_T_0",RooArgList(mean_H_0,Delta_T_0));

  RooRealVar sigma_H_0("sigma_H_0","width of H gaussian background",starting_sigma_H_0,low_sigma_H_0,up_sigma_H_0);
  RooRealVar sigma_L_0("sigma_L_0","width of L gaussian background",starting_sigma_L_0,low_sigma_L_0,up_sigma_L_0);
  RooRealVar sigma_T_0("sigma_T_0","width of T gaussian background",starting_sigma_T_0,low_sigma_T_0,up_sigma_T_0);
  
  if(CB_model){
    RooGaussian PDF_L_0("PDF_L_0","gaussian L_0",x,mean_L_0,sigma_L_0) ;//,alpha_CB_L,n_CB_L) ;
    //RooCBShape PDF_L_0("PDF_L_0","gaussian L_0",x,mean_L_0,sigma_L_0,alpha_CB_L,n_CB_L) ;
    RooCBShape PDF_H_0("PDF_H_0","gaussian H_0",x,mean_H_0,sigma_H_0,alpha_CB_H,n_CB_H) ;
  }else{
    RooGaussian PDF_L_0("PDF_L_0","gaussian L_0",x,mean_L_0,sigma_L_0) ;//,alpha_CB_L,n_CB_L) ;
    RooGaussian PDF_H_0("PDF_H_0","gaussian H_0",x,mean_H_0,sigma_H_0) ;//,alpha_CB_H,n_CB_H) ;
  }
  RooCBShape PDF_T_0("PDF_T_0","gaussian T_0",x,mean_T_0,sigma_T_0,alpha_CB_T,n_CB_T) ;  
  RooRealVar alpha_0("alpha_0","alpha_0",starting_alpha_0,low_alpha_0,up_alpha_0);
  RooRealVar beta_0("beta_0","beta_0",starting_beta_0,low_beta_0,up_beta_0);
  RooFormulaVar Frac_L_0("Frac_L_0","Frac_L_0","alpha_0",RooArgList(alpha_0));
  RooFormulaVar Frac_H_0("Frac_H_0","Frac_H_0","beta_0-alpha_0*beta_0",RooArgList(beta_0,alpha_0));
  RooFormulaVar Frac_T_0("Frac_T_0","Frac_T_0","1-alpha_0-beta_0+alpha_0*beta_0",RooArgList(alpha_0,beta_0));
  RooArgList  pdfList_sig_0(PDF_L_0,PDF_H_0); if(add_third_signal_pos0) {pdfList_sig_0.add(PDF_T_0);}//
  RooArgList  fracList_sig_0(alpha_0); if(add_third_signal_pos0) {fracList_sig_0.add(beta_0);}
  RooAddPdf   PDF_sig_0("PDF_sig_0","PDF_sig_0",pdfList_sig_0,fracList_sig_0,kTRUE);
  
  RooRealVar a0_0("a0_0", "", 0.0, -10, 10);
  RooRealVar a1_0("a1_0", "", 0.0, -20, 20);
  RooRealVar a2_0("a2_0", "", 0.0015, -20, 20);
  RooArgList  coeffList_sig_0(a0_0);
  if(bkg_Chebychev_polynomial_degree>=1){
    coeffList_sig_0.add(a1_0);    
    if(bkg_Chebychev_polynomial_degree>=2){
      coeffList_sig_0.add(a2_0);
    }
  }
  if(!add_background_component){
  a0_0.setConstant(kTRUE);
  a1_0.setConstant(kTRUE);
  a2_0.setConstant(kTRUE);
  }
  Delta_H_0.setConstant(kTRUE);
  //Delta_T_0.setConstant(kTRUE);
  RooChebychev PDF_B_0("PDF_B_0","PDF_B_0",x,coeffList_sig_0);
  RooRealVar  Frac_sig_0("Frac_sig_0","fraction of sig events", 0.9, 0.7,1.0);
  RooArgList  pdfList_0(PDF_sig_0,PDF_B_0);
  RooArgList  fracList_0(Frac_sig_0);
  RooAddPdf   model_0("model_0","model_0",pdfList_0,fracList_0,kTRUE);
  RooAddPdf   model_0_b("model_0_b","model_0_b",pdfList_0,fracList_0,kTRUE);

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
    RooFitResult* fit_results_0_b = model_0_b.fitTo(ds_0_H,Save(),Minimizer(Type_minim_pf,Algo_minim_pf),Strategy(2),SumW2Error(kFALSE),PrintLevel(MN_output_print_level_prefit),PrintEvalErrors(-1),Warnings(kFALSE),Verbose(kFALSE));//);
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
  
  
  
  RooFitResult* fit_results_0;
  if(use_NLL){
    RooAbsReal* nll_model_0 = model_0.createNLL(ds_0_H,Extended(kFALSE)) ;
    //RooMinimizer RMN_0 = RooMinimizer(*nll_model_0);
    RooMinimizer RMN_0 (*nll_model_0);
    RMN_0.setErrorLevel(-1);
    RMN_0.setVerbose(kFALSE);
    RMN_0.setPrintEvalErrors(-1);
    RMN_0.setPrintLevel(MN_output_print_level);
    RMN_0.setStrategy(2);
    RMN_0.minimize(Type_minim,Algo_minim);
    fit_results_0=RMN_0.fit("hmrt") ;
  }else{fit_results_0 = model_0.fitTo(ds_0,Save(),Strategy(2),SumW2Error(kFALSE),InitialHesse(true),PrintLevel(MN_output_print_level),PrintEvalErrors(-1),Warnings(kFALSE));//,InitialHesse(true));//,Hesse(kFALSE));//,Extended(kFALSE),Verbose(kFALSE));//Minimizer(Type_minim,Algo_minim),);
  }
  
  fit_results_0->Print("v");
  RooPlot* xframe2_0 = x.frame(Title("Fit")) ;
  ds_0.plotOn(xframe2_0);//,DataError(RooAbsData::SumW2)) ;
  TH2 *h_correlation_0 = fit_results_0->correlationHist(); h_correlation_0->SetTitle("correlation matrix ");h_correlation_0->GetYaxis()->SetLabelSize(0.1); h_correlation_0->GetYaxis()->SetLabelFont(70);h_correlation_0->GetXaxis()->SetLabelSize(0.075); h_correlation_0->GetXaxis()->SetLabelFont(70); h_correlation_0->SetMarkerSize(2);
  model_0.plotOn(xframe2_0,Range("Fit_Range_0"),Components("PDF_L_0"),LineStyle(kDashed),LineColor(1)) ;
  model_0.plotOn(xframe2_0,Range("Fit_Range_0"),Components("PDF_H_0"),LineStyle(kDashed),LineColor(1)) ;
  if(add_third_signal_pos0) model_0.plotOn(xframe2_0,Range("Fit_Range_0"),Components("PDF_T_0"),LineStyle(kDashed),LineColor(1)) ;
  model_0.plotOn(xframe2_0,Range("Fit_Range_0"),Components("PDF_sig_0"),LineStyle(kDashed),LineColor(5)) ;
  model_0.plotOn(xframe2_0,Range("Fit_Range_0"),Components("PDF_B_0"),LineStyle(kDashed),LineColor(2)) ;
  model_0.plotOn(xframe2_0,Range("Fit_Range_0")) ;

  
  // Construct a histogram with the residuals of the data w.r.t. the curve
  RooHist* hresid_0 = xframe2_0->residHist() ;
  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist* hpull_0 = xframe2_0->pullHist() ;
  // Create a new frame to draw the pull distribution and add the distribution to the frame
  RooPlot* xframe4_0 = x.frame(Title("Pull Distribution ")) ;
  xframe4_0->addPlotable(hpull_0,"P") ;
  
  xframe2_0->GetXaxis()->SetLabelSize(0.05);
  //xframe2_0->GetXaxis()->SetLabelFont(70);
  xframe2_0->GetXaxis()->SetTitleSize(0.05);
  //xframe2_0->GetXaxis()->SetTitleFont(70);
  xframe4_0->GetXaxis()->SetLabelSize(0.05);
  //xframe4_0->GetXaxis()->SetLabelFont(70);
  xframe4_0->GetXaxis()->SetTitleSize(0.05);
  //xframe4_0->GetXaxis()->SetTitleFont(70);
  
  
  
  RooPlot* xframe2_0_log = x.frame(Title("Fit")) ;
  ds_0.plotOn(xframe2_0_log);//,DataError(RooAbsData::SumW2)) ;
  model_0.plotOn(xframe2_0_log,Range("Fit_Range_0"),Components("PDF_L_0"),LineStyle(kDashed),LineColor(1)) ;
  model_0.plotOn(xframe2_0_log,Range("Fit_Range_0"),Components("PDF_H_0"),LineStyle(kDashed),LineColor(1)) ;
  if(add_third_signal_pos0) model_0.plotOn(xframe2_0_log,Range("Fit_Range_0"),Components("PDF_T_0"),LineStyle(kDashed),LineColor(1)) ;
  model_0.plotOn(xframe2_0_log,Range("Fit_Range_0"),Components("PDF_sig_0"),LineStyle(kDashed),LineColor(5)) ;
  model_0.plotOn(xframe2_0_log,Range("Fit_Range_0"),Components("PDF_B_0"),LineStyle(kDashed),LineColor(2)) ;
  model_0.plotOn(xframe2_0_log,Range("Fit_Range_0")) ;
  model_0.paramOn(xframe2_0_log);
  model_0.paramOn(xframe2_0);
  xframe2_0_log->GetXaxis()->SetLabelSize(0.05);
  //xframe2_0_log->GetXaxis()->SetLabelFont(70);
  xframe2_0_log->GetXaxis()->SetTitleSize(0.05);
  //xframe2_0_log->GetXaxis()->SetTitleFont(70);







  
  TF1 *line_0s = new TF1("line_0s","0",-100,100);line_0s->SetLineColor(8);
  TF1 *line_1ps = new TF1("line_1ps","1",-100,100);line_1ps->SetLineColor(4);
  TF1 *line_1ns = new TF1("line_1ns","-1",-100,100);line_1ns->SetLineColor(4);
  TF1 *line_2ps = new TF1("line_2ps","2",-100,100);line_2ps->SetLineColor(kOrange-2);
  TF1 *line_2ns = new TF1("line_2ns","-2",-100,100);line_2ns->SetLineColor(kOrange-2);
  TF1 *line_3ps = new TF1("line_3ps","3",-100,100);line_3ps->SetLineColor(2);
  TF1 *line_3ns = new TF1("line_3ns","-3",-100,100);line_3ns->SetLineColor(2);
  if(draw_results){  
    TCanvas* c_Fit = new TCanvas("Fit results","Fit results",0,0,1124,700) ;
    c_Fit->Divide(2,2) ;
    gStyle->SetOptFit(0111); 
    c_Fit->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe2_0->GetYaxis()->SetTitleOffset(1.6) ; xframe2_0->Draw() ;
    c_Fit->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe2_0_log->GetYaxis()->SetTitleOffset(1.6) ; gPad->SetLogy(1); xframe2_0_log->Draw() ;
    c_Fit->cd(3) ; gPad->SetLeftMargin(0.15) ; xframe4_0->GetYaxis()->SetTitleOffset(1.6) ; gPad->SetGridy(); xframe4_0->GetYaxis()->SetRangeUser(-5,5); xframe4_0->Draw() ;
    line_3ns->Draw("same");line_3ps->Draw("same");
    line_2ns->Draw("same");line_2ps->Draw("same");
    line_1ns->Draw("same");line_1ps->Draw("same");
    line_0s->Draw("same");line_0s->Draw("same");
    c_Fit->cd(4) ; h_correlation_0->Draw("colz:text");
    
  
    TCanvas* c_Fit_0 = new TCanvas("c_Fit_0","c_Fit_0");
    xframe2_0->Draw() ;

  
  }








  return xframe2_0;






  
}


void  make_KEK_data_histos_column(int my_column){

  
  TFile *file_input = new TFile("run004881_TBC4855-4858_slot01_digits.root");
  TTree *t_input = (TTree*)file_input->Get("laser");
  TFile *file_input_MC = new TFile("ana_laser_s01_0reso_500k.root");
  TTree *tree_MC = (TTree*)file_input_MC->Get("laser");
  TFile *file_input_MC_ring = new TFile("ana_laser_s01_0reso_ring_500k.root");
  TTree *tree_MC_ring = (TTree*)file_input_MC_ring->Get("laser");
  
  TFile *f_data = new TFile(Form("KEK_data_histos_col_%i.root",my_column),"recreate");
  
  int my_pixelID=-9;
  my_pixelID=my_column;
  
  for(int g=1; g<=8;g++){
    cout<<"doing "<<g<<"th row"<<endl;
    my_pixelID=my_column+64*(g-1);
    Float_t upper_bound_hist=2;
    Int_t n_bins = 200;
    TH1D *h_temp = new TH1D("h_temp","h_temp",500,-50,0);
    TCut cut = Form("-20<time&&time<0&&quality==1&&pixel==%i",my_pixelID);
    t_input->Project("h_temp","time",cut);
    //h_temp->Draw();
    Float_t max_bin = h_temp->GetMaximumBin();
    TAxis *xaxis = h_temp->GetXaxis(); 
    Double_t max_pos = xaxis->GetBinCenter(max_bin);
    //cout<< max_pos <<endl;
    TH1D *h_time = new TH1D("h_time","Time [ns]",n_bins,-1,upper_bound_hist);
    t_input->Project("h_time",Form("time-%f",max_pos),cut);
    /*
      Float_t time_sc;
      Int_t pixelID=-99;
      Int_t quality=-9;
      t_input->SetBranchAddress("time",&time_sc);
      t_input->SetBranchAddress("pixel",&pixelID);
      t_input->SetBranchAddress("quality",&quality);
      
      Int_t n_entries = t_input->GetEntries();
      for(int i=0; i<n_entries; i++){t_input->GetEntry(i); if(-20<time_sc&&time_sc<0&&quality==1&&pixelID==my_pixelID) {h_time->Fill(time_sc-max_pos);}} 
      
      //h_time->Draw("E");
      */
    
   
    TCut cut_MC = Form("propTime<1&&pixel==%i",my_pixelID);
    
    TH1D *h_MC_tot_temp = new TH1D ("h_MC_tot_temp","h_MC_tot_temp",100,0.3,1);
    tree_MC->Project("h_MC_tot_temp","propTime",cut_MC);
    
    Float_t max_bin_MC = h_MC_tot_temp->GetMaximumBin();
    TAxis *xaxis_MC = h_MC_tot_temp->GetXaxis(); 
    Double_t max_pos_MC = xaxis_MC->GetBinCenter(max_bin_MC);
    TH1D *h_MC_tot = new TH1D ("h_MC_tot","h_MC_tot",n_bins,-1,upper_bound_hist);
    tree_MC->Project("h_MC_tot",Form("propTime-%f",max_pos_MC),cut_MC);
    TH1D *h_MC_f1 = new TH1D ("h_MC_f1","h_MC_f1",n_bins,-1,upper_bound_hist);
    TH1D *h_MC_f2 = new TH1D ("h_MC_f2","h_MC_f2",n_bins,-1,upper_bound_hist);
    TH1D *h_MC_f3 = new TH1D ("h_MC_f3","h_MC_f3",n_bins,-1,upper_bound_hist);
    TH1D *h_MC_f4 = new TH1D ("h_MC_f4","h_MC_f4",n_bins,-1,upper_bound_hist);
    TH1D *h_MC_f5 = new TH1D ("h_MC_f5","h_MC_f5",n_bins,-1,upper_bound_hist);
    TH1D *h_MC_f6 = new TH1D ("h_MC_f6","h_MC_f6",n_bins,-1,upper_bound_hist);  
    TH1D *h_MC_f7 = new TH1D ("h_MC_f7","h_MC_f7",n_bins,-1,upper_bound_hist);
    TH1D *h_MC_f8 = new TH1D ("h_MC_f8","h_MC_f8",n_bins,-1,upper_bound_hist);  
    TH1D *h_MC_f9 = new TH1D ("h_MC_f9","h_MC_f9",n_bins,-1,upper_bound_hist);
    tree_MC->Project("h_MC_f1",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==1&&pixel==%i",my_pixelID));
    tree_MC->Project("h_MC_f2",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==2&&pixel==%i",my_pixelID));
    tree_MC->Project("h_MC_f3",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==3&&pixel==%i",my_pixelID));
    tree_MC->Project("h_MC_f4",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==4&&pixel==%i",my_pixelID));
    tree_MC->Project("h_MC_f5",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==5&&pixel==%i",my_pixelID));
    tree_MC->Project("h_MC_f6",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==6&&pixel==%i",my_pixelID));
    tree_MC->Project("h_MC_f7",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==7&&pixel==%i",my_pixelID));
    tree_MC->Project("h_MC_f8",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==8&&pixel==%i",my_pixelID));
    tree_MC->Project("h_MC_f9",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==9&&pixel==%i",my_pixelID));
    h_MC_tot->SetLineColor(1);
    h_MC_f1->SetLineColor(11);
    h_MC_f2->SetLineColor(2);
    h_MC_f3->SetLineColor(3);
    h_MC_f4->SetLineColor(4);
    h_MC_f5->SetLineColor(13);
    h_MC_f6->SetLineColor(6);
    h_MC_f7->SetLineColor(7);
    h_MC_f8->SetLineColor(8);
    h_MC_f9->SetLineColor(9);
    h_MC_f1->Scale(h_time->GetMaximum()/h_MC_tot->GetMaximum());
    h_MC_f2->Scale(h_time->GetMaximum()/h_MC_tot->GetMaximum());
    h_MC_f3->Scale(h_time->GetMaximum()/h_MC_tot->GetMaximum());
    h_MC_f4->Scale(h_time->GetMaximum()/h_MC_tot->GetMaximum());
    h_MC_f5->Scale(h_time->GetMaximum()/h_MC_tot->GetMaximum());
    h_MC_f6->Scale(h_time->GetMaximum()/h_MC_tot->GetMaximum());
    h_MC_f7->Scale(h_time->GetMaximum()/h_MC_tot->GetMaximum());
    h_MC_f8->Scale(h_time->GetMaximum()/h_MC_tot->GetMaximum());
    h_MC_f9->Scale(h_time->GetMaximum()/h_MC_tot->GetMaximum());
    h_MC_tot->Scale(h_time->GetMaximum()/h_MC_tot->GetMaximum());


    //////////////////////MC with ring model //////////////////
    TH1D *h_MC_ring_tot = new TH1D ("h_MC_ring_tot","h_MC_ring_tot",n_bins,-1,upper_bound_hist);
    tree_MC_ring->Project("h_MC_ring_tot",Form("propTime-%f",max_pos_MC),cut_MC);
    TH1D *h_MC_ring_f1 = new TH1D ("h_MC_ring_f1","h_MC_ring_f1",n_bins,-1,upper_bound_hist);
    TH1D *h_MC_ring_f2 = new TH1D ("h_MC_ring_f2","h_MC_ring_f2",n_bins,-1,upper_bound_hist);
    TH1D *h_MC_ring_f3 = new TH1D ("h_MC_ring_f3","h_MC_ring_f3",n_bins,-1,upper_bound_hist);
    TH1D *h_MC_ring_f4 = new TH1D ("h_MC_ring_f4","h_MC_ring_f4",n_bins,-1,upper_bound_hist);
    TH1D *h_MC_ring_f5 = new TH1D ("h_MC_ring_f5","h_MC_ring_f5",n_bins,-1,upper_bound_hist);
    TH1D *h_MC_ring_f6 = new TH1D ("h_MC_ring_f6","h_MC_ring_f6",n_bins,-1,upper_bound_hist);  
    TH1D *h_MC_ring_f7 = new TH1D ("h_MC_ring_f7","h_MC_ring_f7",n_bins,-1,upper_bound_hist);
    TH1D *h_MC_ring_f8 = new TH1D ("h_MC_ring_f8","h_MC_ring_f8",n_bins,-1,upper_bound_hist);  
    TH1D *h_MC_ring_f9 = new TH1D ("h_MC_ring_f9","h_MC_ring_f9",n_bins,-1,upper_bound_hist);
    tree_MC_ring->Project("h_MC_ring_f1",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==1&&pixel==%i",my_pixelID));
    tree_MC_ring->Project("h_MC_ring_f2",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==2&&pixel==%i",my_pixelID));
    tree_MC_ring->Project("h_MC_ring_f3",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==3&&pixel==%i",my_pixelID));
    tree_MC_ring->Project("h_MC_ring_f4",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==4&&pixel==%i",my_pixelID));
    tree_MC_ring->Project("h_MC_ring_f5",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==5&&pixel==%i",my_pixelID));
    tree_MC_ring->Project("h_MC_ring_f6",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==6&&pixel==%i",my_pixelID));
    tree_MC_ring->Project("h_MC_ring_f7",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==7&&pixel==%i",my_pixelID));
    tree_MC_ring->Project("h_MC_ring_f8",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==8&&pixel==%i",my_pixelID));
    tree_MC_ring->Project("h_MC_ring_f9",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==9&&pixel==%i",my_pixelID));
    h_MC_ring_tot->SetLineColor(1);
    h_MC_ring_f1->SetLineColor(11);
    h_MC_ring_f2->SetLineColor(2);
    h_MC_ring_f3->SetLineColor(3);
    h_MC_ring_f4->SetLineColor(4);
    h_MC_ring_f5->SetLineColor(13);
    h_MC_ring_f6->SetLineColor(6);
    h_MC_ring_f7->SetLineColor(7);
    h_MC_ring_f8->SetLineColor(8);
    h_MC_ring_f9->SetLineColor(9);
    h_MC_ring_tot->SetMarkerStyle(20);
    h_MC_ring_f1->SetMarkerStyle(20);
    h_MC_ring_f2->SetMarkerStyle(20);
    h_MC_ring_f3->SetMarkerStyle(20);
    h_MC_ring_f4->SetMarkerStyle(20);
    h_MC_ring_f5->SetMarkerStyle(20);
    h_MC_ring_f6->SetMarkerStyle(20);
    h_MC_ring_f7->SetMarkerStyle(20);
    h_MC_ring_f8->SetMarkerStyle(20);
    h_MC_ring_f9->SetMarkerStyle(20);
    h_MC_ring_f1->Scale(h_time->GetMaximum()/h_MC_ring_tot->GetMaximum());
    h_MC_ring_f2->Scale(h_time->GetMaximum()/h_MC_ring_tot->GetMaximum());
    h_MC_ring_f3->Scale(h_time->GetMaximum()/h_MC_ring_tot->GetMaximum());
    h_MC_ring_f4->Scale(h_time->GetMaximum()/h_MC_ring_tot->GetMaximum());
    h_MC_ring_f5->Scale(h_time->GetMaximum()/h_MC_ring_tot->GetMaximum());
    h_MC_ring_f6->Scale(h_time->GetMaximum()/h_MC_ring_tot->GetMaximum());
    h_MC_ring_f7->Scale(h_time->GetMaximum()/h_MC_ring_tot->GetMaximum());
    h_MC_ring_f8->Scale(h_time->GetMaximum()/h_MC_ring_tot->GetMaximum());
    h_MC_ring_f9->Scale(h_time->GetMaximum()/h_MC_ring_tot->GetMaximum());
    h_MC_ring_tot->Scale(h_time->GetMaximum()/h_MC_ring_tot->GetMaximum());
    
    TCanvas *c = new TCanvas("c","c");
    h_time->Draw("E");
    h_MC_tot->Draw("same");
    h_MC_f1->Draw("same");
    h_MC_f2->Draw("same");
    h_MC_f3->Draw("same");
    h_MC_f4->Draw("same");
    h_MC_f5->Draw("same");
    h_MC_f6->Draw("same");
    h_MC_f7->Draw("same");
    h_MC_f8->Draw("same");
    h_MC_f9->Draw("same");
    /*
    h_MC_ring_tot->Draw("same");
    h_MC_ring_f1->Draw("same");
    h_MC_ring_f2->Draw("same");
    h_MC_ring_f3->Draw("same");
    h_MC_ring_f4->Draw("same");
    h_MC_ring_f5->Draw("same");
    h_MC_ring_f6->Draw("same");
    h_MC_ring_f7->Draw("same");
    h_MC_ring_f8->Draw("same");
    h_MC_ring_f9->Draw("same");
    */
    
    
    f_data->cd();
    f_data->mkdir(Form("histos_%i",my_pixelID));
    f_data->cd(Form("histos_%i",my_pixelID));
    h_time->Write();
    h_MC_tot->Write();
    h_MC_f1->Write();
    h_MC_f2->Write();
    h_MC_f3->Write();
    h_MC_f4->Write();
    h_MC_f5->Write();
    h_MC_f6->Write();
    h_MC_f7->Write();
    h_MC_f8->Write();
    h_MC_f9->Write();
    h_MC_ring_tot->Write();
    h_MC_ring_f1->Write();
    h_MC_ring_f2->Write();
    h_MC_ring_f3->Write();
    h_MC_ring_f4->Write();
    h_MC_ring_f5->Write();
    h_MC_ring_f6->Write();
    h_MC_ring_f7->Write();
    h_MC_ring_f8->Write();
    h_MC_ring_f9->Write();
    c->Write();
    f_data->cd();
    delete h_temp;
    delete h_time;
    delete h_MC_tot;
    delete h_MC_tot_temp;
    delete h_MC_f1;
    delete h_MC_f2;
    delete h_MC_f3;
    delete h_MC_f4;
    delete h_MC_f5;
    delete h_MC_f6;
    delete h_MC_f7;
    delete h_MC_f8;
    delete h_MC_f9;
    delete h_MC_ring_tot;
    delete h_MC_ring_f1;
    delete h_MC_ring_f2;
    delete h_MC_ring_f3;
    delete h_MC_ring_f4;
    delete h_MC_ring_f5;
    delete h_MC_ring_f6;
    delete h_MC_ring_f7;
    delete h_MC_ring_f8;
    delete h_MC_ring_f9;
    delete c;
  }
    f_data->Close();
    delete f_data;
}


void loop_fit_column(int column_number){

  TCanvas *cc = new TCanvas("cc","cc");
  cc->Divide(1,8);


  for(int h=1;h<=8;h++){
    
    RooPlot* frame_temp=NULL;
    frame_temp = Fit_KEK_data(true,column_number,h);
    cc->cd(9-h);
    frame_temp->Draw() ;
    //delete frame_temp;
  }
  
}
