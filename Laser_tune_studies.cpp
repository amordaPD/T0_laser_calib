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
  int n_POIs;
  RooPlot* xframe2_amplitude_0;
  RooPlot* xframe2_fit_0;
  RooPlot* xframe2_pull_0;
  TH2 *h_correlation_0;
  
};



//TString input_filename ="runScan-T74-pos0_large_scale_out";


Fit_results Fit_head(string _draw_results="draw", TString _fiber, TString _tune  ){





  

  TFile *f_input_histogram = new TFile("raw_data/fbk/fbk-B31.8-"+_fiber+"-"+_tune+"-times.root"); 
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
  //bool fix_params = true;
  vector<float> POIs;
  POIs.clear();
  // S e t u p   m o d e l 
  // ---------------------
  
  double my_low_x=13.9;//2.4;
  double my_up_x=14.2;
  RooRealVar x("time","Time [ns]",my_low_x,my_up_x) ;
  RooRealVar amp("amplitude","Amplitude [ADC counts]",0,10000) ;

  ///////fixing starting values and boundaries for two positions


  

  TTree *tree = (TTree*)f_input_histogram->Get("times");
  RooDataSet ds_0("ds_0","ds_0", RooArgSet(x,amp),Import(*tree));//,Cut(Form("Amplitude>%d",amplitude_cut)));
  


  double starting_mean_L_0;
  double low_mean_L_0;
  double up_mean_L_0;
 
  double starting_sigma_L_0;
  double low_sigma_L_0;
  double up_sigma_L_0;




  
  double starting_alpha_CB;
  double low_alpha_CB;
  double up_alpha_CB;
  double starting_n_CB;
  double low_n_CB;
  double up_n_CB;
  
  
  /////POS 0
  //  low_x_0=my_low_x; up_x_0=my_up_x;
  
  
  low_sigma_L_0=0.01;
  up_sigma_L_0=0.500;
  
  starting_mean_L_0=13.725;
  starting_sigma_L_0=0.0410;
  
  low_mean_L_0=starting_mean_L_0-0.35;
  up_mean_L_0=starting_mean_L_0+0.35;
  
  
  starting_alpha_CB=-1.5;//-0.5;
  starting_n_CB=4;
  
  low_alpha_CB=-5;
  up_alpha_CB=0.0;
  
  low_n_CB=0;
  up_n_CB=20;
  
  
  RooRealVar alpha_CB("alpha_CB","alpha parameter of CB",   starting_alpha_CB,low_alpha_CB,up_alpha_CB);
  RooRealVar     n_CB("n_CB",    "exponential decay of CB ",starting_n_CB,low_n_CB,up_n_CB);
  RooRealVar mean_L_0("mean_L_0","mean of L gaussian background pos 0",starting_mean_L_0,low_mean_L_0,up_mean_L_0);
  RooRealVar sigma_L_0("sigma_L_0","width of L gaussian background",starting_sigma_L_0,low_sigma_L_0,up_sigma_L_0);
  RooCBShape PDF_L_0("PDF_L_0","gaussian L_0",x,mean_L_0,sigma_L_0,alpha_CB,n_CB) ;
  RooCBShape model_0("model_0","gaussian L_0",x,mean_L_0,sigma_L_0,alpha_CB,n_CB) ;
  RooCBShape model_0_b("model_0_b","gaussian L_0",x,mean_L_0,sigma_L_0,alpha_CB,n_CB) ;
  
  RooRealVar a0_0("a0_0", "", 0.0, -10, 10);
  RooRealVar a1_0("a1_0", "", 0.0, -20, 20);
  RooRealVar a2_0("a2_0", "", 0.0015, -20, 20);
  RooChebychev PDF_B_0("PDF_B_0","PDF_B_0",x,RooArgList(a0_0,a1_0));//,a2_0));//
  /*
  RooRealVar  Frac_sig_0("Frac_sig_0","fraction of sig events", 0.9, 0.7,1.0);
  RooAddPdf   model_0("model_0","model_0",RooArgList(PDF_L_0,PDF_B_0),RooArgList(Frac_sig_0),kTRUE);
  RooAddPdf   model_0_b("model_0_b","model_0_b",RooArgList(PDF_L_0,PDF_B_0),RooArgList(Frac_sig_0),kTRUE);
  
  
  */
  
 
  
  
  
  
  
  
  
  
  
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
    RooFitResult* fit_results_0_b = model_0_b.fitTo(ds_0,Save(),Minimizer(Type_minim_pf,Algo_minim_pf),Strategy(2),SumW2Error(kFALSE),PrintLevel(MN_output_print_level_prefit),PrintEvalErrors(-1),Warnings(kFALSE),Verbose(kFALSE));//);
    sigma_L_0.setConstant(kFALSE);
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
    RooAbsReal* nll_model_0 = model_0.createNLL(ds_0,Extended(kFALSE));
    RooMinimizer RMN_0 = RooMinimizer(*nll_model_0);
    RMN_0.setErrorLevel(-1);
    RMN_0.setVerbose(kFALSE);
    RMN_0.setPrintEvalErrors(-1);
    RMN_0.setPrintLevel(MN_output_print_level);
    RMN_0.setStrategy(2);
    RMN_0.minimize(Type_minim,Algo_minim);
    fit_results_0=RMN_0.fit("hmr") ;
  }else{
    fit_results_0 = model_0.fitTo(ds_0,Save(),Strategy(2),SumW2Error(kFALSE),InitialHesse(true),PrintLevel(MN_output_print_level),PrintEvalErrors(-1),Warnings(kFALSE));
  }
  
  fit_results_0->Print("v");
  RooPlot* xframe2_0 = x.frame(Title(_fiber+_tune)) ;
  ds_0.plotOn(xframe2_0);//,DataError(RooAbsData::SumW2)) ;
  TH2 *h_correlation_0 = fit_results_0->correlationHist(); h_correlation_0->SetTitle("correlation matrix");h_correlation_0->GetYaxis()->SetLabelSize(0.1); h_correlation_0->GetYaxis()->SetLabelFont(70); h_correlation_0->SetMarkerSize(2);
  model_0.plotOn(xframe2_0,Components("PDF_L_0"),LineStyle(kDashed),LineColor(1)) ;
  model_0.plotOn(xframe2_0,Components("PDF_B_0"),LineStyle(kDashed),LineColor(2)) ;
  model_0.plotOn(xframe2_0) ;
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
  RooPlot* xframe4_0 = x.frame(Title(_fiber+_tune)) ;
  xframe4_0->addPlotable(hpull_0,"P") ;
  
  xframe2_0->GetXaxis()->SetLabelSize(0.05);
  xframe2_0->GetXaxis()->SetLabelFont(70);
  xframe2_0->GetXaxis()->SetTitleSize(0.05);
  xframe2_0->GetXaxis()->SetTitleFont(70);
  xframe4_0->GetXaxis()->SetLabelSize(0.05);
  xframe4_0->GetXaxis()->SetLabelFont(70);
  xframe4_0->GetXaxis()->SetTitleSize(0.05);
  xframe4_0->GetXaxis()->SetTitleFont(70);
  
  
  
  RooPlot* xframe2_0_log = x.frame(Title(_fiber+_tune)) ;
  ds_0.plotOn(xframe2_0_log);//,DataError(RooAbsData::SumW2)) ;
  model_0.plotOn(xframe2_0_log,Components("PDF_L_0"),LineStyle(kDashed),LineColor(1)) ;
  model_0.plotOn(xframe2_0_log,Components("PDF_B_0"),LineStyle(kDashed),LineColor(2)) ;
  model_0.plotOn(xframe2_0_log) ;
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
  double T0_res(-9);
  double T0_res_err(-9);
  
  
  T0=mean_L_0.getVal();
  T0_err=mean_L_0.getError();
  T0_res=sigma_L_0.getVal();
  T0_res_err= sigma_L_0.getError();
  
  ///////SINGLE FIBER FIT VALUES////////
  /*0   */   POIs.push_back(T0); 
  /*1   */   POIs.push_back(T0_err);
  /*2   */   POIs.push_back(T0_res);
  /*3   */   POIs.push_back(T0_res_err);
  
  
  
  
    /*
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
  
    */
  
  
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
  
  
  
  
  for(int g=0; g<POIs.size();g++){
    my_fit_results.Results[g]=POIs[g];
  }

  my_fit_results.n_POIs=POIs.size();
  my_fit_results.xframe2_fit_0=xframe2_0;
  my_fit_results.xframe2_pull_0=xframe4_0;
  my_fit_results.h_correlation_0=h_correlation_0;
  
  
  return my_fit_results;
  
  //return POIs;
}



vector<float> loop_channels(int deep_fixed_params,TString input_tune, bool plot_summaries){ //rel_weight is the relative weight of the 1 to 0 pois, i.e. total dataset= pos.0 + rel_weight*pos.1 - in the real case rel_weight=1.0
  
  
  TCanvas *c_pos0_AllChannels = new TCanvas("c_pos0_AllChannels","c_pos0_AllChannels",0,0,1124,700);
  TCanvas *c_pos0_AllChannels_pulls = new TCanvas("c_pos0_AllChannels_pulls","c_pos0_AllChannels_pulls",0,0,1124,700);
  TCanvas *c_pos0_AllChannels_corr = new TCanvas("c_pos0_AllChannels_corr","c_pos0_AllChannels_corr",0,0,1124,700);
  TCanvas *c_pos0_AllChannels_amp = new TCanvas("c_pos0_AllChannels_amp","c_pos0_AllChannels_amp",0,0,1124,700);
  
  
  
  vector<float> relative_weight_Frac_0; relative_weight_Frac_0.clear();
  c_pos0_AllChannels->Divide(6,5);
  
  c_pos0_AllChannels_pulls->Divide(6,5);
  c_pos0_AllChannels_corr->Divide(6,5);
  c_pos0_AllChannels_amp->Divide(6,5);
  
  
  
  TString _tunes[30];
  _tunes[0]="T50";
  _tunes[1]="T52";
  _tunes[2]="T54";
  _tunes[3]="T56";
  _tunes[4]="T58";
  _tunes[5]="T60";
  _tunes[6]="T62";
  _tunes[7]="T64";
  _tunes[8]="T66";
  _tunes[9]="T68";
  _tunes[10]="T70";
  _tunes[11]="T71";
  _tunes[12]="T72";
  _tunes[13]="T73";
  _tunes[14]="T74";
  _tunes[15]="T75";
  _tunes[16]="T76";
  _tunes[17]="T77";
  _tunes[18]="T78";
  _tunes[19]="T79";
  _tunes[20]="T80";
  _tunes[21]="T81";
  _tunes[22]="T82";
  _tunes[23]="T84";
  _tunes[24]="T86";
  _tunes[25]="T88";
  _tunes[26]="T90";
  _tunes[27]="T92";
  _tunes[28]="T94";
  _tunes[29]="T96";
  
  
  
  
  
  
  
  vector<float> POIs;
  vector<float> FRCs;
  vector<float> AMPs;
  vector<float> SPSs;//Secondary peaks separations
  float      mean_L_0_B[30];
  float  err_mean_L_0_B[30];
  float     sigma_L_0_B[30];
  float err_sigma_L_0_B[30];
  
  
  ////////////////////////FIT RESULTS WITH 3 GAUSSIAN SIGNALS/////////////////
  
  
  float x[30];
  float err_x[30];


  x[0]=50;
  x[1]=52;
  x[2]=54;
  x[3]=56;
  x[4]=58;
  x[5]=60;
  x[6]=62;
  x[7]=64;
  x[8]=66;
  x[9]=68;
  x[10]=70;
  x[11]=71;
  x[12]=72;
  x[13]=73;
  x[14]=74;
  x[15]=75;
  x[16]=76;
  x[17]=77;
  x[18]=78;
  x[19]=79;
  x[20]=80;
  x[21]=81;
  x[22]=82;
  x[23]=84;
  x[24]=86;
  x[25]=88;
  x[26]=90;
  x[27]=92;
  x[28]=94;
  x[29]=96;

  err_x[0]=0;
  err_x[1]=0;
  err_x[2]=0;
  err_x[3]=0;
  err_x[4]=0;
  err_x[5]=0;
  err_x[6]=0;
  err_x[7]=0;
  err_x[8]=0;
  err_x[9]=0;
  err_x[10]=0;
  err_x[11]=0;
  err_x[12]=0;
  err_x[13]=0;
  err_x[14]=0;
  err_x[15]=0;
  err_x[16]=0;
  err_x[17]=0;
  err_x[18]=0;
  err_x[19]=0;
  err_x[20]=0;
  err_x[21]=0;
  err_x[22]=0;
  err_x[23]=0;
  err_x[24]=0;
  err_x[25]=0;
  err_x[26]=0;
  err_x[27]=0;
  err_x[28]=0;
  err_x[29]=0;
  
  for(Int_t i=0; i<30; i++){
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
  }
  
  
  if(plot_summaries){
    
    
    
    TF1 *line = new TF1("line","0.100",-1,16);
    line->SetTitle("100 ps");
    // TLine *line = new TLine(0,100,16,100);
    line->SetLineColor(2);
 
    
    TGraphErrors *MEAN_L_B_0 = new TGraphErrors(16,x,mean_L_0_B,err_x,err_mean_L_0_B);MEAN_L_B_0->SetName("MEAN_L_B_0"); MEAN_L_B_0->SetTitle("T_{0} - fit pos.0");
   
    
    TGraphErrors *SIGMA_L_B_0 = new TGraphErrors(16,x,sigma_L_0_B,err_x,err_sigma_L_0_B);SIGMA_L_B_0->SetName("SIGMA_L_B_0"); SIGMA_L_B_0->SetTitle("#delta t_{low} pos.0 - fit pos.0");

    MEAN_L_B_0->SetMarkerStyle(20);
    MEAN_L_B_0->SetMarkerSize(2);
    MEAN_L_B_0->SetMarkerColor(4);
    SIGMA_L_B_0->SetMarkerStyle(20);
    SIGMA_L_B_0->SetMarkerSize(2);
    SIGMA_L_B_0->SetMarkerColor(1);
    
    TMultiGraph *mg_MEANS_L_0 = new TMultiGraph();mg_MEANS_L_0->SetName("mg_MEANS_L_0");
    mg_MEANS_L_0->SetTitle("T_{0} vs Channel - pos.0");
    mg_MEANS_L_0->Add(MEAN_L_B_0);
  
    
    TCanvas *c_absolute_positions_01_b = new TCanvas("c_absolute_position_01_b","c_absolute_position_01_b",0,0,1124,700);
    mg_MEANS_L_0->Draw("AP");
    mg_MEANS_L_0->GetXaxis()->SetTitle("Channel");
    mg_MEANS_L_0->GetYaxis()->SetTitle("T_{0} [ns]");
    mg_MEANS_L_0->GetYaxis()->SetRangeUser(22.75,23.2);
    gPad->BuildLegend();
    gPad->Update();
    gPad->SetGridy();
  
    
    TMultiGraph *mg_SIGMA_B_L = new TMultiGraph();mg_SIGMA_B_L->SetName("mg_SIGMA_B_L");
    mg_SIGMA_B_L->SetTitle("Time resolution (#delta t) low time vs Channel");
    mg_SIGMA_B_L->Add(REF_sigma_L);
    mg_SIGMA_B_L->Add(SIGMA_L_B_0);
    //mg_SIGMA_B_L->Add(SIGMA_L_A_0);
 
    
    TCanvas* c_Res_pos0 = new TCanvas("c_Res_pos0","c_Res_pos0",0,0,1124,700);
    c_Res_pos0->cd();
    mg_SIGMA_L_0->Draw("AP");
    line->Draw("same");
    gPad->Update();
    mg_SIGMA_L_0->GetXaxis()->SetTitle("Channel");
    mg_SIGMA_L_0->GetYaxis()->SetTitle("#delta t [ns]");
    gPad->Update();
    gPad->BuildLegend();

    
    c_pos0_AllChannels->Draw();    c_pos0_AllChannels_pulls->Draw();  c_pos0_AllChannels_corr->Draw(); 
  
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













