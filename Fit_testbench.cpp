////////TIME CALIBRATION CODE///////****
///to run:
///gROOT->ProcessLine(".L Fit_testbench.cpp"); loop_channels(2,true); > ridaje.log


#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooMinimizer.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooMinuit.h"
#include "TAxis.h"
#include <vector>

using namespace RooFit ;
vector<int> Fits_status;
struct Fit_results{
  float Results[80];
  int n_POIs;
  RooPlot* xframe2_fit_0;
  RooPlot* xframe2_pull_0;
  TH2 *h_correlation_0;
  RooPlot* xframe2_fit_1;
  RooPlot* xframe2_pull_1;
  TH2 *h_correlation_1;
  RooPlot* xframe2_fit_01;
  RooPlot* xframe2_pull_01;
  TH2 *h_correlation_01;
  
};


//  ofstream cout("fits_report.txt");

TFile *f_input_histogram_noPrism = new TFile("./flat_ntuples/noprism-T50-V2490_out.root");
TFile *f_input_histogram_pos0    = new TFile("./flat_ntuples/run100317-2-T77-nuovalente-p0_out.root");
TFile *f_input_histogram_pos1    = new TFile("./flat_ntuples/run100317-3-T77-nuovalente-p1_out.root");
TFile *f_input_histogram_full_ds = new TFile("./flat_ntuples/run100317-4-T77-nuovalente-p0e1_out.root");

//vector<float> Fit_head(string _draw_results, int fix_params, int ch ){
Fit_results Fit_head(string _draw_results="draw", int fix_params=2, int ch =0 ){ 
  bool do_prefit=true;
  bool do_prefit_fullSpectrum = true;
  bool use_NLL=true; //to set the use of fitTo method of RooAbsPdf or the explicit construction of the nll  ///true recomended
  int MN_output_print_level=-1;
  int MN_output_print_level_prefit;
  bool print_prefit_info=false;
  const char * Type_minim="Minuit";
  const char * Algo_minim="minimize";//
  const char * Type_minim_pf="Minuit";//"Minuit2";//
  const char * Algo_minim_pf="minimize";//"scan";//
  bool do_PixByPix_CBparams_fit = true;
  bool do_simultaneous_fit=false;
  bool add_third_signal=false;
  bool simulate_CB_tail=false;
  bool fit_real_FiberCombs_data=true;
  bool binned_fit = true; //recomended true if you don't want to wait 7 minutes fot the fit output
  bool compute_FWHM = true;
  int bkg_Chebychev_polynomial_degree=0;//set to n to have  degree n+1!!!!!!!!!
  int amplitude_cut = -40;

  
  
  if(!print_prefit_info){MN_output_print_level_prefit=-1;}else{MN_output_print_level_prefit=MN_output_print_level;}
  bool draw_results;
  if(_draw_results=="draw"){draw_results=true;}else if(_draw_results=="blind"){draw_results=false;}else{draw_results=false;}
  Fit_results my_fit_results;
  //TString channel_0 = Form("fiber%d-10-%d",0,ch);
  TString hist_channel_0 = Form("fiber0-0-%d",ch);
  TString hist_channel_1 = Form("fiber0-0-%d",ch);
  //bool fix_params = true;
  vector<float> POIs;
  POIs.clear();


  
  // S e t u p   m o d e l 
  // ---------------------
 
  // double my_low_x=8;
  //double my_up_x=11;
  double my_low_x=22;
  double my_up_x=25;
  RooRealVar x("Time","Time [ns]",my_low_x,my_up_x) ;
  RooRealVar amp("Amplitude","Amplitude [ADC counts]",-200,0) ;
  RooRealVar CH("Channel","PMT Channel",0,16) ;
  RooRealVar weight("weight","Event weight",0,1) ;



  TTree *tree_0 = (TTree*)f_input_histogram_pos0->Get("tree_input");
  TTree *tree_1 = (TTree*)f_input_histogram_pos1->Get("tree_input");
  TTree *tree_ds = (TTree*)f_input_histogram_full_ds->Get("tree_input");
  //RooDataSet ds_0_amp("ds_0_amp","ds_0_amp", RooArgSet(x,amp,CH),Import(*tree_0),Cut(Form("Amplitude<%d &&Channel==%d",amplitude_cut,ch)));
  //RooDataSet ds_1_amp("ds_1_amp","ds_1_amp", RooArgSet(x,amp,CH),Import(*tree_1),Cut(Form("Amplitude<%d &&Channel==%d",amplitude_cut,ch)));
  
  if(binned_fit){
    
    TH1 *h_input_histogram_0 = (TH1*)f_input_histogram_pos0->Get(hist_channel_0);
    RooDataHist ds_0("ds_0","ds_0",RooArgSet(x),Import(*h_input_histogram_0)) ;
    TH1 *h_input_histogram_1 = (TH1*)f_input_histogram_pos1->Get(hist_channel_1);
    RooDataHist ds_1("ds_1","ds_1",RooArgSet(x),Import(*h_input_histogram_1)) ;
    cout<<Form("dataset 0 ch %d info :",ch)<<ds_0.Print("v")<<endl;
    cout<<Form("dataset 1 ch %d info :",ch)<<ds_1.Print("v")<<endl;
        
    TH1*h_input_histogram;
    if(fit_real_FiberCombs_data){
      h_input_histogram = (TH1*)f_input_histogram_full_ds->Get(Form("fiber0-0-%d",ch));
    }else{
      h_input_histogram = (TH1*)f_input_histogram_pos0->Get(hist_channel_0);
      h_input_histogram->Add(h_input_histogram_1,1);
    }
    RooDataHist DS("DS","DS",RooArgSet(x),Import(*h_input_histogram)) ;
    
  }else{
    RooDataSet ds_0("ds_0","ds_0", RooArgSet(x,amp,CH),Import(*tree_0),Cut(Form("Amplitude<%d &&Channel==%d",amplitude_cut,ch)));//,WeightVar("weight"));//
    RooDataSet ds_1("ds_1","ds_1", RooArgSet(x,amp,CH),Import(*tree_1),Cut(Form("Amplitude<%d &&Channel==%d",amplitude_cut,ch)));//,WeightVar("weight"));
    if(fit_real_FiberCombs_data){
      RooDataSet DS("DS","DS", RooArgSet(x,amp,CH),Import(*tree_0),Cut(Form("Amplitude<%d &&Channel==%d",amplitude_cut,ch)));//WeightVar("weight"));
      DS.append(ds_1);
    }else{
      RooDataSet DS("DS","DS", RooArgSet(x,amp,CH),Import(*tree_ds),Cut(Form("Amplitude<%d &&Channel==%d",amplitude_cut,ch)));//,WeightVar("weight"));
    }
  }

  



  
  ///////fixing starting values and boundaries for two positions






  double low_x_0;
  double up_x_0;
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



  double low_x_1;
  double up_x_1;
  double starting_mean_L_1;
  double      low_mean_L_1;
  double       up_mean_L_1;
  double starting_delta_H_1;
  double      low_delta_H_1;
  double       up_delta_H_1;
  double starting_delta_T_1;
  double      low_delta_T_1;
  double       up_delta_T_1;
  double starting_sigma_T_1;
  double      low_sigma_T_1;
  double       up_sigma_T_1;
  double starting_sigma_L_1;
  double      low_sigma_L_1;
  double       up_sigma_L_1;
  double starting_sigma_H_1;
  double      low_sigma_H_1;
  double       up_sigma_H_1; 
  double starting_Frac_H_1;
  double      low_Frac_L_1;
  double       up_Frac_L_1;
  double starting_Frac_L_1;
  double      low_Frac_H_1;
  double       up_Frac_H_1;
  double starting_Frac_T_1;
  double      low_Frac_T_1;
  double       up_Frac_T_1;
  double starting_alpha_1;
  double      low_alpha_1;
  double       up_alpha_1;
  double starting_beta_1;
  double      low_beta_1;
  double       up_beta_1;







  
  double starting_alpha_CB;
  double low_alpha_CB;
  double up_alpha_CB;
  double starting_n_CB;
  double low_n_CB;
  double up_n_CB;
  
  
  /////POS 0
   low_x_0=my_low_x; up_x_0=my_up_x;
  
  low_delta_H_0=0.180;
   up_delta_H_0=0.5;
  
  low_delta_T_0=-0.2;
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

  
   starting_mean_L_0=22.7;//8.5;
   if(ch==3||ch==7) starting_mean_L_0=22.6;//8.4;
   starting_delta_H_0=0.280;
   starting_delta_T_0=0.200;
   starting_sigma_L_0=0.080;
   starting_sigma_H_0=0.080;
   starting_sigma_T_0=0.1000;
   starting_alpha_0=0.5;
   starting_beta_0=0.5;
  
  low_mean_L_0=starting_mean_L_0-0.15;
   up_mean_L_0=starting_mean_L_0+0.15;



   //////POS 1 ////////
   low_x_1=my_low_x; up_x_1=my_up_x;
  
  low_delta_H_1=0.180;
   up_delta_H_1=0.5;
  
  low_delta_T_1=0;
   up_delta_T_1=0.5;
  
  low_sigma_L_1=0.035;
   up_sigma_L_1=0.5200;
  
  low_sigma_H_1=0.045;
   up_sigma_H_1=0.5150;
  
  low_sigma_T_1=0.050;
   up_sigma_T_1=0.500;
  
  low_alpha_1=0.01;
   up_alpha_1=0.8;
  
  low_beta_1=0.5;
   up_beta_1=1.0;

  
   starting_mean_L_1=8.7;
   starting_delta_H_1=0.300;
  starting_delta_T_1=0.3;
  starting_sigma_L_1=0.080;
  starting_sigma_H_1=0.080;
  starting_sigma_T_1=0.120;
  starting_alpha_1=0.5;
   starting_beta_1=0.5;
 
  if(ch==4){starting_mean_L_1=8.5;}
  if(ch==5){starting_mean_L_1=8.5;}
  if(ch==7){starting_mean_L_1=8.5;}
  if(ch==10){starting_mean_L_1=8.5;}
  if(ch==11){starting_mean_L_1=8.5;}
  if(ch==12){starting_mean_L_1=8.5;}
  if(ch==13){starting_mean_L_1=8.5;}
  if(ch==14){starting_mean_L_1=8.5;}
  if(ch==15){starting_mean_L_1=8.5;starting_delta_H_1=0.2;}
    


  
  low_mean_L_1=starting_mean_L_1-0.15;
   up_mean_L_1=starting_mean_L_1+0.15;


  starting_alpha_CB=-0.5;
  starting_n_CB=4;

  low_alpha_CB=-5;
   up_alpha_CB=0.0;
   
  low_n_CB=0;
   up_n_CB=10;

  
  RooRealVar alpha_CB("alpha_CB","alpha parameter of CB",   starting_alpha_CB,low_alpha_CB,up_alpha_CB);
  RooRealVar     n_CB("n_CB",    "exponential decay of CB ",starting_n_CB,low_n_CB,up_n_CB);


  if(do_PixByPix_CBparams_fit){
    
    
    RooRealVar xR("Time","Time [ns]",26.0,27.8) ;//26.2,27.2) ;
    RooRealVar meanR("meanR","mean of L gaussian background pos 0",26.5,26,27);
    RooRealVar sigmaR("sigmaR","width of H gaussian background",0.060,0.0350,0.180);
    RooRealVar alphaR("alphaR","alphaR",-1.4,-5,0);
    RooRealVar nR("nR","nR",2,0,10);
    RooCBShape  modelR_sig("modelR_sig","gaussian",xR,meanR,sigmaR,alphaR,nR) ;
    
    
    
    RooRealVar a0R_bkg("a0R_bkg", "", 0.166, -100, 100);
    RooRealVar a1R_bkg("a1R_bkg", "", -0.027, -10, 10);
    RooChebychev modelR_bkg("modelR_bkg","PDFR_bkg",xR,RooArgList(a0R_bkg,a1R_bkg));//,a2_bkg));//
    
    
    RooRealVar  FracR("FracR","fraction of sig events", 0.8, 0.05,1.0);
    RooArgList  pdfListR(modelR_sig,modelR_bkg);
    RooArgList  fracListR(FracR);
    RooAddPdf   modelR("modelR","modelR",pdfListR,fracListR);
    
    TH1 *h_input_histogramR = (TH1*)f_input_histogram_noPrism->Get(Form("fiber0-0-%d",ch));
    RooDataHist dsR("dsR","dsR",RooArgSet(xR),Import(*h_input_histogramR)) ;
    RooFitResult* fit_resultsR = modelR.fitTo(dsR,Save(),Extended(kFALSE),SumW2Error(kFALSE),PrintLevel(-1),PrintEvalErrors(-1),Warnings(kFALSE),Verbose(kFALSE));
    
    starting_alpha_CB=alphaR.getVal();
    starting_n_CB=nR.getVal(); 


    
    if(draw_results){
      RooPlot* xframe2R = xR.frame() ;
      dsR.plotOn(xframe2R);//,DataError(RooAbsData::SumW2)) ;
      modelR.plotOn(xframe2R,Components("modelR_sig"),LineStyle(kDashed),LineColor(8)) ;
      modelR.plotOn(xframe2R,Components("modelR_bkg"),LineStyle(kDashed),LineColor(2)) ;
      modelR.plotOn(xframe2R,LineColor(4)) ;
      RooHist* hpullR = xframe2R->pullHist() ;
      // Create a new frame to draw the pull distribution and add the distribution to the frame
      RooPlot* xframe4R = xR.frame(Title("Pull Distribution")) ;
      xframe4R->addPlotable(hpullR,"P") ;
      
      
      TCanvas* c_FitR = new TCanvas("Fit results single signal","Fit results single signal",0,0,1124,700) ;
      c_FitR->Divide(1,2);
      gStyle->SetOptFit(1111); 
      modelR.paramOn(xframe2R);
      c_FitR->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe2R->GetYaxis()->SetTitleOffset(1.6) ; xframe2R->Draw() ;
      TPaveLabel *t1R = new TPaveLabel(0.7,0.6,0.9,0.68, Form("#chi^{2}/nDOF = %f", xframe2R->chiSquare()));
      t1R->Draw("same");
      
      c_FitR->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe2R->GetYaxis()->SetTitleOffset(1.6) ; gPad->SetGridy(); xframe4R->Draw() ;
      fit_resultsR->Print("v");
    }
    
  }




  
  
  RooRealVar x_0("Time","Time [ns]",low_x_0,up_x_0);//eee
  RooRealVar mean_L_0("mean_L_0","mean of L gaussian background pos 0",starting_mean_L_0,low_mean_L_0,up_mean_L_0);
  RooRealVar Delta_H_0("Delta_H_0","Delta H pos 0",starting_delta_H_0,low_delta_H_0,up_delta_H_0);
  RooRealVar Delta_T_0("Delta_T_0","Delta T pos 0",starting_delta_T_0,low_delta_T_0,up_delta_T_0);
  RooFormulaVar mean_H_0("mean_H_0","mean_H_0","mean_L_0+Delta_H_0",RooArgList(mean_L_0,Delta_H_0));
  RooFormulaVar mean_T_0("mean_T_0","mean_T_0","mean_H_0+Delta_T_0",RooArgList(mean_H_0,Delta_T_0));
  RooRealVar sigma_L_0("sigma_L_0","width of L gaussian background",starting_sigma_L_0,low_sigma_L_0,up_sigma_L_0);
  RooRealVar sigma_H_0("sigma_H_0","width of H gaussian background",starting_sigma_H_0,low_sigma_H_0,up_sigma_H_0);
  RooRealVar sigma_T_0("sigma_T_0","width of T gaussian background",starting_sigma_T_0,low_sigma_T_0,up_sigma_T_0);
  /*RooGaussian PDF_L_0("PDF_L_0","gaussian L_0",x,mean_L_0,sigma_L_0) ;
  RooGaussian PDF_H_0("PDF_H_0","gaussian H_0",x,mean_H_0,sigma_H_0) ;
  RooGaussian PDF_T_0("PDF_T_0","gaussian T_0",x,mean_T_0,sigma_T_0) ;*/
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
  RooAddPdf   PDF_sig_0("PDF_sig_0","PDF_sig_0",pdfList_sig_0,fracList_sig_0,kTRUE);


  float starting_delta_means_L(0.118);
  float      low_delta_means_L(0.000);
  float       up_delta_means_L(0.500);
  if(ch==1||ch==2){starting_delta_means_L=0.1;}
  RooRealVar delta_means_L("delta_means_L","delta_means_L",starting_delta_means_L,low_delta_means_L,up_delta_means_L);//0.118,0.0,0.5);


  
  RooRealVar x_1("Time","Time [ns]",low_x_1,up_x_1) ;//
  RooFormulaVar mean_L_1("mean_L_1","mean_L_1","mean_L_0+delta_means_L",RooArgList(mean_L_0,delta_means_L));
  //  RooRealVar mean_L_1("mean_L_1","mean of L gaussian background pos 0",starting_mean_L_1,low_mean_L_1,up_mean_L_1);
  RooRealVar Delta_H_1("Delta_H_1","Delta H pos 1",starting_delta_H_1,low_delta_H_1,up_delta_H_1);
  RooRealVar Delta_T_1("Delta_T_1","Delta T pos 1",starting_delta_T_1,low_delta_T_1,up_delta_T_1);
  RooFormulaVar mean_H_1("mean_H_1","mean_H_1","mean_L_1+Delta_H_1",RooArgList(mean_L_1,Delta_H_1));
  RooFormulaVar mean_T_1("mean_T_1","mean_T_1","mean_H_1+Delta_T_1",RooArgList(mean_H_1,Delta_T_1));
  RooRealVar sigma_L_1("sigma_L_1","width of L gaussian background",starting_sigma_L_1,low_sigma_L_1,up_sigma_L_1);
  RooRealVar sigma_H_1("sigma_H_1","width of H gaussian background",starting_sigma_H_1,low_sigma_H_1,up_sigma_H_1);
  RooRealVar sigma_T_1("sigma_T_1","width of T gaussian background",starting_sigma_T_1,low_sigma_T_1,up_sigma_T_1);
  /*RooGaussian PDF_L_1("PDF_L_1","gaussian L_1",x,mean_L_1,sigma_L_1) ;
  RooGaussian PDF_H_1("PDF_H_1","gaussian H_1",x,mean_H_1,sigma_H_1) ; 
  RooGaussian PDF_T_1("PDF_T_1","gaussian T_1",x,mean_T_1,sigma_T_1) ;*/
  RooCBShape PDF_L_1("PDF_L_1","gaussian L_1",x,mean_L_1,sigma_L_1,alpha_CB,n_CB) ;
  RooCBShape PDF_H_1("PDF_H_1","gaussian H_1",x,mean_H_1,sigma_H_1,alpha_CB,n_CB) ; 
  RooCBShape PDF_T_1("PDF_T_1","gaussian T_1",x,mean_T_1,sigma_T_1,alpha_CB,n_CB) ;
  RooRealVar alpha_1("alpha_1","alpha_1",starting_alpha_1,low_alpha_1,up_alpha_1);
  RooRealVar beta_1("beta_1","beta_1",starting_beta_1,low_beta_1,up_beta_1);
  RooFormulaVar Frac_L_1("Frac_L_1","Frac_L_1","alpha_1",RooArgList(alpha_1));
  RooFormulaVar Frac_H_1("Frac_H_1","Frac_H_1","beta_1-alpha_1*beta_1",RooArgList(beta_1,alpha_1));
  RooFormulaVar Frac_T_1("Frac_T_1","Frac_T_1","1-alpha_1-beta_1+alpha_1*beta_1",RooArgList(alpha_1,beta_1));
  RooArgList  pdfList_sig_1(PDF_L_1,PDF_H_1);if(add_third_signal) {pdfList_sig_1.add(PDF_T_1);}//
  RooArgList  fracList_sig_1(alpha_1);if(add_third_signal) {fracList_sig_1.add(beta_1);}
  RooAddPdf   PDF_sig_1("PDF_sig_1","PDF_sig_1",pdfList_sig_1,fracList_sig_1,kTRUE);
  //RooGenericPdf PDF_sig_1("PDF_sig_1","PDF_sig_1","alpha_1*PDF_L_1+(1-alpha_1)*PDF_H_1",RooArgList(x,alpha_1,PDF_L_1,PDF_H_1));
  //RooGenericPdf PDF_sig_1("PDF_sig_1","PDF_sig_1","@0*@1+(1-@0)*@2",RooArgList(alpha_1,PDF_L_1,PDF_H_1));
  
  RooRealVar a0_0("a0_0", "", 0.0, -10, 10);
  RooRealVar a1_0("a1_0", "", 0.0, -20, 20);
  RooRealVar a2_0("a2_0", "", 0.0015, -20, 20);
  RooArgList  coeffList_sig_0(a0_0);

  RooRealVar a0_1("a0_1", "", 0.0, -10, 10);
  RooRealVar a1_1("a1_1", "", 0.0, -20, 20);
  RooRealVar a2_1("a2_1", "", 0.0015, -20, 20);
  RooArgList  coeffList_sig_1(a0_1);

  RooRealVar a0_bkg("a0_bkg", "", 0.0, -10, 10);
  RooRealVar a1_bkg("a1_bkg", "", 0.0, -10, 10);
  RooRealVar a2_bkg("a2_bkg", "", 0.0015, -10, 10);
  RooArgList  coeffList_sig_bkg(a0_bkg);

  if(bkg_Chebychev_polynomial_degree>=1){
    coeffList_sig_0.add(a1_0);
    coeffList_sig_1.add(a1_1);
    coeffList_sig_bkg.add(a1_bkg);
    
    if(bkg_Chebychev_polynomial_degree>2){
      coeffList_sig_0.add(a2_0);
      coeffList_sig_1.add(a2_1);
      coeffList_sig_bkg.add(a2_bkg);
    }
  }
  RooChebychev PDF_B_0("PDF_B_0","PDF_B_0",x,coeffList_sig_0);//,a2_0));//
  RooChebychev PDF_B_1("PDF_B_1","PDF_B_1",x,coeffList_sig_1);//,a2_1));//
  RooChebychev PDF_bkg("PDF_bkg","PDF_bkg",x,coeffList_sig_bkg);//,a2_bkg));//
  RooRealVar  Frac_sig_0("Frac_sig_0","fraction of sig events", 0.9, 0.7,1.0);
  RooArgList  pdfList_0(PDF_sig_0,PDF_B_0);//
  RooArgList  fracList_0(Frac_sig_0);
  RooAddPdf   model_0("model_0","model_0",pdfList_0,fracList_0,kTRUE);

  
  RooRealVar  Frac_sig_1("Frac_sig_1","fraction of sig events", 0.9, 0.7,1.0);
  RooArgList  pdfList_1(PDF_sig_1,PDF_B_1);
  RooArgList  fracList_1(Frac_sig_1);
  RooAddPdf   model_1("model_1","model_1",pdfList_1,fracList_1,kTRUE);

  RooRealVar  Frac_0("Frac_0"      ,"fraction of pos 0 events", 0.3, 0.0,0.7);
  RooArgList  pdfList_sig(PDF_sig_0,PDF_sig_1);//
  RooArgList  fracList_sig(Frac_0);
  RooAddPdf   PDF_sig("PDF_sig","PDF_sig",pdfList_sig,fracList_sig,kTRUE);

  RooRealVar  Frac_sig("Frac_sig"      ,"fraction of signal events", 0.9, 0.8,1.0);
  RooArgList  pdfList(PDF_sig,PDF_bkg);//
  RooArgList  fracList(Frac_sig);
  RooAddPdf   model("model","model",pdfList,fracList);


  
  RooAddPdf   model_0_b("model_0_b","model_0_b",pdfList_0,fracList_0,kTRUE);
  RooAddPdf   model_1_b("model_1_b","model_1_b",pdfList_1,fracList_1,kTRUE);
  RooAddPdf   model_b("model_b","model_b",pdfList,fracList);








  
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
  



  
  /*
  if(do_simultaneous_fit){
    RooDataHist ds_0_s ("ds_0_s","ds_0_s",RooArgSet(x),Import(*h_input_histogram_0)) ;
    RooDataHist ds_1_s ("ds_1_s","ds_1_s",RooArgSet(x),Import(*h_input_histogram_1)) ;
    RooCategory sample("sample","sample") ;
    sample.defineType("fib_0") ;
    sample.defineType("fib_1") ;
    RooDataHist combData("combData","combined data",x,Index(sample),Import("fib_0",*h_input_histogram_0),Import("fib_1",*h_input_histogram_1)) ;
    RooSimultaneous sim_Pdf("sim_Pdf","simultaneous pdf",sample) ;
    
    // Associate model with the physics state and model_ctl with the control state
    sim_Pdf.addPdf(model_0,"fib_0") ;
    sim_Pdf.addPdf(model_1,"fib_1") ;
    
    
    
    // P e r f o r m   a   s i m u l t a n e o u s   f i t
    // ---------------------------------------------------
    
    // Perform simultaneous fit of model to data and model_ctl to data_ctl
    RooFitResult* fit_res_sim=sim_Pdf.fitTo(combData,Save()) ;
    TH2 *h_correlation_sim = fit_res_sim->correlationHist(); h_correlation_sim->SetTitle(Form("correlation matrix 0 ch %d",ch));h_correlation_sim->GetYaxis()->SetLabelSize(0.1); h_correlation_sim->GetYaxis()->SetLabelFont(70); h_correlation_sim->SetMarkerSize(2);
    
    
    RooPlot* frame1_simultaneous = x.frame(Title("Fiber 0")) ;
    
    // Plot all data tagged as physics sample
    combData.plotOn(frame1_simultaneous,Cut("sample==sample::fib_0")) ;
    
    // Plot "physics" slice of simultaneous pdf. 
    // NBL You _must_ project the sample index category with data using ProjWData 
    // as a RooSimultaneous makes no prediction on the shape in the index category 
    // and can thus not be integrated
    sim_Pdf.plotOn(frame1_simultaneous,Slice(sample,"fib_0"),ProjWData(sample,combData)) ;
    sim_Pdf.plotOn(frame1_simultaneous,Slice(sample,"fib_0"),Components("PDF_sig_0"),ProjWData(sample,combData),LineStyle(kDashed)) ;
    sim_Pdf.plotOn(frame1_simultaneous,Slice(sample,"fib_0"),Components("PDF_B_0"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(2));
    sim_Pdf.plotOn(frame1_simultaneous,Slice(sample,"fib_0"),Components("PDF_L_0"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(8));
    sim_Pdf.plotOn(frame1_simultaneous,Slice(sample,"fib_0"),Components("PDF_H_0"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(8));
    
    // The same plot for the control sample slice
    RooPlot* frame2_simultaneous = x.frame(Title("Fiber 1")) ;
    combData.plotOn(frame2_simultaneous,Cut("sample==sample::fib_1")) ;
    sim_Pdf.plotOn(frame2_simultaneous,Slice(sample,"fib_1"),ProjWData(sample,combData)) ;
    sim_Pdf.plotOn(frame2_simultaneous,Slice(sample,"fib_1"),Components("PDF_sig_1"),ProjWData(sample,combData),LineStyle(kDashed)) ;
    sim_Pdf.plotOn(frame2_simultaneous,Slice(sample,"fib_1"),Components("PDF_B_0"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(2));
    sim_Pdf.plotOn(frame2_simultaneous,Slice(sample,"fib_1"),Components("PDF_L_1"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(8));
    sim_Pdf.plotOn(frame2_simultaneous,Slice(sample,"fib_1"),Components("PDF_H_1"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(8));
    
    model_0.paramOn(frame1_simultaneous);
    model_1.paramOn(frame2_simultaneous);
    
    RooHist* hpull_0_sim = frame1_simultaneous->pullHist() ;
    // Create a new frame to draw the pull distribution and add the distribution to the frame
    RooPlot* xframe4_0_sim = x.frame(Title("Pull Distribution 0")) ;
    xframe4_0_sim->addPlotable(hpull_0_sim,"P") ;
    RooHist* hpull_1_sim = frame2_simultaneous->pullHist() ;
    // Create a new frame to draw the pull distribution and add the distribution to the frame
    RooPlot* xframe4_1_sim = x.frame(Title("Pull Distribution 1")) ;
    xframe4_1_sim->addPlotable(hpull_1_sim,"P") ;
    
    TCanvas* c_simultaneous = new TCanvas("c_simultaneouspdf","c_simultaneouspdf",800,400) ;
    c_simultaneous->Divide(2,3) ;
    c_simultaneous->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1_simultaneous->GetYaxis()->SetTitleOffset(1.4) ; frame1_simultaneous->Draw() ;
    c_simultaneous->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2_simultaneous->GetYaxis()->SetTitleOffset(1.4) ; frame2_simultaneous->Draw() ;
    c_simultaneous->cd(3) ; gPad->SetLeftMargin(0.15) ; h_correlation_sim->Draw("colz:text");
    c_simultaneous->cd(4) ; xframe4_0_sim->Draw();
    c_simultaneous->cd(5) ; xframe4_0_sim->Draw();
    
    
    } */
  
  // D o  F i t 
  // ------------------------
  
  
  
  
  
  
  
  double low_x_01;
  double up_x_01;
  
  low_x_01=TMath::Min(low_x_0,low_x_1);
  up_x_01=TMath::Max(up_x_0,up_x_1);
  x.setRange("Fit_Range",low_x_01,up_x_01);
  
  x.setRange("Fit_Range_0",low_x_0,up_x_0);
  x.setRange("Fit_Range_1",low_x_1,up_x_1);
  
  
  
  alpha_CB.setVal(starting_alpha_CB);
  n_CB.setVal(starting_n_CB);
  if(do_PixByPix_CBparams_fit){  n_CB.setConstant(kTRUE);}// alpha_CB.setConstant(kTRUE);}//



  
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
  
  
  
  RooFitResult* fit_results_0;
  if(use_NLL){
    RooAbsReal* nll_model_0 = model_0.createNLL(ds_0,Extended(kFALSE)) ;
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
  RooPlot* xframe2_0 = x.frame(Title(Form("pos. 0, channel %d",ch))) ;
  ds_0.plotOn(xframe2_0);//,DataError(RooAbsData::SumW2)) ;
  TH2 *h_correlation_0 = fit_results_0->correlationHist(); h_correlation_0->SetTitle(Form("correlation matrix 0 ch %d",ch));h_correlation_0->GetYaxis()->SetLabelSize(0.1); h_correlation_0->GetYaxis()->SetLabelFont(70); h_correlation_0->SetMarkerSize(2);
  model_0.plotOn(xframe2_0,Range("Fit_Range_0"),Components("PDF_L_0"),LineStyle(kDashed),LineColor(1)) ;
  model_0.plotOn(xframe2_0,Range("Fit_Range_0"),Components("PDF_H_0"),LineStyle(kDashed),LineColor(1)) ;
  if(add_third_signal) model_0.plotOn(xframe2_0,Range("Fit_Range_0"),Components("PDF_T_0"),LineStyle(kDashed),LineColor(1)) ;
  model_0.plotOn(xframe2_0,Range("Fit_Range_0"),Components("PDF_sig_0"),LineStyle(kDashed),LineColor(5)) ;
  model_0.plotOn(xframe2_0,Range("Fit_Range_0"),Components("PDF_B_0"),LineStyle(kDashed),LineColor(2)) ;
  model_0.plotOn(xframe2_0,Range("Fit_Range_0")) ;
  // Construct a histogram with the residuals of the data w.r.t. the curve
  RooHist* hresid_0 = xframe2_0->residHist() ;
  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist* hpull_0 = xframe2_0->pullHist() ;
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
  if(do_PixByPix_CBparams_fit){  n_CB.setConstant(kTRUE); }//alpha_CB.setConstant(kTRUE);}
  
  
  
  
  
  /*
    RooAbsReal* nll_model_1 = model_1.createNLL(ds_1,Save(),Range("Fit_Range_1"),SumW2Error(kFALSE),Extended(kFALSE)) ;
    RooMinuit RMN_1 = RooMinuit(*nll_model_1);
    RMN_1.setWarnLevel(0);
    RMN_1.setVerbose(kFALSE);
    mean_L_1.setConstant(kTRUE);
    Delta_H_1.setConstant(kTRUE);
    sigma_L_1.setConstant(kTRUE);
    sigma_H_1.setConstant(kTRUE);
    RMN_1.migrad() ;
    mean_L_1.setConstant(kFALSE);
    Delta_H_1.setConstant(kFALSE);
    sigma_L_1.setConstant(kFALSE);
    sigma_H_1.setConstant(kFALSE);
    RooFitResult* fit_results_1=RMN_1.fit("shr") ;
    
    
  */
  
  
  mean_L_0.setConstant(kTRUE);
  if(do_prefit){
    
    cout<<"___________________________________"<<endl;
    cout<<"                                   "<<endl;
    cout<<"-----------------------------------"<<endl;
    cout<<"***********************************"<<endl;
    cout<<"-----------------------------------"<<endl;
    cout<<"___________________________________"<<endl;
    cout<<"STATUS MIGRAD EDM 1 0"<<endl;
    
    //Delta_H_1.setConstant(kTRUE);
    sigma_L_1.setConstant(kTRUE);
    sigma_H_1.setConstant(kTRUE);
    RooFitResult* fit_results_1_b = model_1_b.fitTo(ds_1,Save(),Minimizer(Type_minim_pf,Algo_minim_pf),Strategy(2),SumW2Error(kFALSE),PrintLevel(MN_output_print_level_prefit),PrintEvalErrors(-1),Warnings(kFALSE),Verbose(kFALSE));//);
    sigma_L_1.setConstant(kFALSE);
    sigma_H_1.setConstant(kFALSE);
    //Delta_H_1.setConstant(kFALSE);
    if(print_prefit_info) fit_results_1_b->Print("v");
  }
  
  cout<<"___________________________________"<<endl;
  cout<<"                                   "<<endl;
  cout<<"-----------------------------------"<<endl;
  cout<<"***********************************"<<endl;
  cout<<"-----------------------------------"<<endl;
  cout<<"___________________________________"<<endl;
  cout<<"STATUS MIGRAD EDM 1"<<endl;
  
  
  RooFitResult* fit_results_1;
  if(use_NLL){
    RooAbsReal* nll_model_1 = model_1.createNLL(ds_1,Extended(kFALSE)) ;
    //RooMinimizer RMN_1 = RooMinimizer(*nll_model_1);
    RooMinimizer RMN_1(*nll_model_1);
    RMN_1.setErrorLevel(-1);
    RMN_1.setVerbose(kFALSE);
    RMN_1.setPrintEvalErrors(-1);
    RMN_1.setPrintLevel(MN_output_print_level);
    RMN_1.setStrategy(2);
    RMN_1.minimize(Type_minim,Algo_minim);
    fit_results_1=RMN_1.fit("hmrt") ;
  }else{
    fit_results_1 = model_1.fitTo(ds_1,Save(),InitialHesse(true),Strategy(2),SumW2Error(kFALSE),PrintLevel(MN_output_print_level),PrintEvalErrors(-1),Warnings(kFALSE));//,InitialHesse(true));//,Hesse(kFALSE));//,Extended(kFALSE),Verbose(kFALSE));//,Minimizer(Type_minim,Algo_minim));
  }
  fit_results_1->Print("v");
  RooPlot* xframe2_1 = x.frame(Title(Form("pos. 1, channel %d",ch))) ;
  ds_1.plotOn(xframe2_1);
  TH2 *h_correlation_1 = fit_results_1->correlationHist(); h_correlation_1->SetTitle(Form("correlation matrix 1 ch %d",ch));h_correlation_1->GetYaxis()->SetLabelSize(0.1); h_correlation_1->GetYaxis()->SetLabelFont(70); h_correlation_1->SetMarkerSize(2);
  model_1.plotOn(xframe2_1,Range("Fit_Range_1"),Components("PDF_L_1"),LineStyle(kDashed),LineColor(8)) ;
  model_1.plotOn(xframe2_1,Range("Fit_Range_1"),Components("PDF_H_1"),LineStyle(kDashed),LineColor(8)) ;
  if(add_third_signal) model_1.plotOn(xframe2_1,Range("Fit_Range_1"),Components("PDF_T_1"),LineStyle(kDashed),LineColor(8)) ;
  model_1.plotOn(xframe2_1,Range("Fit_Range_1"),Components("PDF_sig_1"),LineStyle(kDashed),LineColor(5)) ;
  model_1.plotOn(xframe2_1,Range("Fit_Range_1"),Components("PDF_B_1"),LineStyle(kDashed),LineColor(2)) ;
  model_1.plotOn(xframe2_1,Range("Fit_Range_1")) ;
  // Construct a histogram with the residuals of the data w.r.t. the curve
  RooHist* hresid_1 = xframe2_1->residHist() ;
  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist* hpull_1 = xframe2_1->pullHist() ;
  // Create a new frame to draw the pull distribution and add the distribution to the frame
  RooPlot* xframe4_1 = x.frame(Title(Form("Pull Distribution fiber 1 - ch %d",ch)));
  xframe4_1->addPlotable(hpull_1,"P") ;
  
  
  
  xframe2_1->GetXaxis()->SetLabelSize(0.05);
  xframe2_1->GetXaxis()->SetLabelFont(70);
  xframe2_1->GetXaxis()->SetTitleSize(0.05);
  xframe2_1->GetXaxis()->SetTitleFont(70);
  xframe4_1->GetXaxis()->SetLabelSize(0.05);
  xframe4_1->GetXaxis()->SetLabelFont(70);
  xframe4_1->GetXaxis()->SetTitleSize(0.05);
  xframe4_1->GetXaxis()->SetTitleFont(70);
  
  RooPlot* xframe2_1_log = x.frame(Title(Form("channel %d",ch))) ;
  ds_1.plotOn(xframe2_1_log);//,DataError(RooAbsData::SumW2)) ;
  model_1.plotOn(xframe2_1_log,Range("Fit_Range_1"),Components("PDF_L_1"),LineStyle(kDashed),LineColor(8)) ;
  model_1.plotOn(xframe2_1_log,Range("Fit_Range_1"),Components("PDF_H_1"),LineStyle(kDashed),LineColor(8)) ;
   if(add_third_signal)  model_1.plotOn(xframe2_1_log,Range("Fit_Range_1"),Components("PDF_T_1"),LineStyle(kDashed),LineColor(8)) ;
  model_1.plotOn(xframe2_1_log,Range("Fit_Range_1"),Components("PDF_sig_1"),LineStyle(kDashed),LineColor(5)) ;
  model_1.plotOn(xframe2_1_log,Range("Fit_Range_1"),Components("PDF_B_1"),LineStyle(kDashed),LineColor(2)) ;
  model_1.plotOn(xframe2_1_log,Range("Fit_Range_1")) ;
  Fits_status.push_back(fit_results_1->status());
  model_1.paramOn(xframe2_1_log);
  model_1.paramOn(xframe2_1);
  
  xframe2_1_log->GetXaxis()->SetLabelSize(0.05);
  xframe2_1_log->GetXaxis()->SetLabelFont(70);
  xframe2_1_log->GetXaxis()->SetTitleSize(0.05);
  xframe2_1_log->GetXaxis()->SetTitleFont(70);



  if(compute_FWHM){
    double Res;
    Double_t HW_CB;
    double absAlpha; 
    double n;
    if(absAlpha>=sqrt(2*log(2))) {
      HW_CB=sqrt(2*log(2));
    } else {
      absAlpha=TMath::Abs(alpha_CB.getVal());
      n= n_CB.getVal();
      Double_t A = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
      Double_t B= n/absAlpha - absAlpha;
      HW_CB = -B + TMath::Power(2*A,1/n); ////////the proper formula is multiplied by -1
      cout<<" eddaje "<<absAlpha<<"  "<<n<<"  "<<A<<"  "<<B<<"  "<<TMath::Power(2*A,1/n)<<"  "<<HW_CB<<endl;
    }
    Res=sigma_L_0.getVal();
    sigma_L_0.setVal(Res*0.5*(sqrt(2*log(2))+HW_CB));
    Res=sigma_L_1.getVal();
    sigma_L_1.setVal(Res*0.5*(sqrt(2*log(2))+HW_CB));
    Res=sigma_H_0.getVal();
    sigma_H_0.setVal(Res*0.5*(sqrt(2*log(2))+HW_CB));
    Res=sigma_H_1.getVal();
    sigma_H_1.setVal(Res*0.5*(sqrt(2*log(2))+HW_CB));
  }

  
  ///////SINGLE FIBER FIT VALUES////////
  /*0   */   POIs.push_back(mean_L_0.getVal()); 
  /*1   */   POIs.push_back(mean_L_0.getError());
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
  //  cout<<"eeeei  Delta H "<<Delta_H_0.getVal()<<" +-  "<<Delta_H_0.getError()<<endl;
  /*16   */  POIs.push_back(delta_means_L.getVal()); 
  /*17   */  POIs.push_back(delta_means_L.getError());//getError());
  /*18   */  POIs.push_back(Delta_H_1.getVal());
  /*19   */  POIs.push_back(Delta_H_1.getError());
  /*20   */  POIs.push_back(Delta_T_1.getVal());
  /*21   */  POIs.push_back(Delta_T_1.getError());
  /*22   */  POIs.push_back(sigma_L_1.getVal());
  /*23   */  POIs.push_back(sigma_L_1.getError());
  /*24   */  POIs.push_back(sigma_H_1.getVal());
  /*25   */  POIs.push_back(sigma_H_1.getError());
  /*26  */   POIs.push_back(sigma_T_1.getVal());
  /*27  */   POIs.push_back(sigma_T_1.getError());
  /*28  */   POIs.push_back(alpha_1.getVal());
  /*29  */   POIs.push_back(alpha_1.getError());
  /*30  */   POIs.push_back(beta_1.getVal());
  /*31  */   POIs.push_back(beta_1.getError());
  
  
  
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
  cout<<"eeeih Delta H 1   ch : "<<ch<<"   "<<Delta_H_1.getVal()<<" +- "<<Delta_H_1.getError()<<endl;
  cout<<"eeeih Delta H 0 1 ch : "<<ch<<"   "<<delta_means_L.getVal()<<" +- "<<delta_means_L.getError()<<endl;
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
  
  cout<<" #### Fit quality checks, channel "<<ch<<" , fiber position 1 ####"<<endl;
  cout<<"------- General quality parameters ---------"<<endl;
  cout<<"Chi2/ndof  : "<<  xframe2_1->chiSquare()<<endl;
  cout<<"EDM : "<<fit_results_1->edm()<<endl;
  cout<<"Fit status : "<<fit_results_1->status()<<endl;
  cout<<"-------------  Hit boundaries --------------"<<endl;
  if(((POIs[18+2]-low_delta_H_1)<0.0001) || ((-POIs[18+2]+up_delta_H_1)<0.001))   cout<<"delta H_1 hits boundary"<<endl;
  //if(((POIs[18+4]-low_delta_T_1)<0.0001) || ((-POIs[18+4]+up_delta_T_1)<0.001))   cout<<"delta T_1 hits boundary"<<endl;
  if(((POIs[18+6]-low_sigma_L_1)<0.0001) || ((-POIs[18+6]+up_sigma_L_1)<0.001))   cout<<"sigma L_1 hits boundary"<<endl;
  if(((POIs[18+8]-low_sigma_H_1)<0.0001) || ((-POIs[18+8]+up_sigma_H_1)<0.001))   cout<<"sigma H_1 hits boundary"<<endl;
  //if(((POIs[18+10]-low_sigma_T_1)<0.0001)|| ((-POIs[18+10]+up_sigma_T_1)<0.001))  cout<<"sigma T_1 hits boundary"<<endl;
  if(((POIs[18+12]-low_alpha_1)<0.0001)  || ((-POIs[18+12]+up_alpha_1)<0.001))    cout<<"alpha_1  hits boundary"<<endl;
  //if(((POIs[18+14]-low_beta_1)<0.0001)   || ((-POIs[18+14]+up_beta_1)<0.001))     cout<<"beta_1 hits boundary"<<endl;
  cout<<"--------------- Large/small errors --------------"<<endl;
  for(int h=0;h<8;h++){
    if(h!=2&&h!=5&&h!=7){
      if(TMath::Abs(POIs[16+2*h+1]/POIs[16+2*h])>0.50)    cout<<"large error at "<<PAR_NAMES[h]<<"_1"<<endl;
      if(TMath::Abs(POIs[16+2*h+1]/POIs[16+2*h])<0.00001) cout<<"small error at "<<PAR_NAMES[h]<<"_1"<<endl;
    }
  }
  cout<<"--------------- Correlation matrix --------------"<<endl;
  for(int ff=1;ff<h_correlation_1->GetXaxis()->GetNbins()+1;ff++){
    for(int gg=1;gg<h_correlation_1->GetXaxis()->GetNbins()-ff+1;gg++){
      if(TMath::Abs(h_correlation_1->GetBinContent(ff,gg))<0.001||TMath::Abs(h_correlation_1->GetBinContent(ff,gg))>0.99){
	cout<<"low/high "<<h_correlation_1->GetXaxis()->GetBinLabel(ff)<<"-"<<h_correlation_1->GetYaxis()->GetBinLabel(gg)<<" correlation : "<<h_correlation_1->GetBinContent(ff,gg)<< endl;
      }
    }
  }
  


  if(simulate_CB_tail){
    // TH1 *h_input_histogram_1_cut = (TH1*)f_input_histogram->Get(channel_1);
    //RooDataHist *ds_1_toy_core=
    // RooDataHist ds_1_s ("ds_1_s","ds_1_s",RooArgSet(x),Import(*h_input_histogram_1)) ;
    //x.setBins(110) ;
    ///    taking the core of the distribution
    RooDataHist *ds_1_toy_core=(RooDataHist*)ds_1.reduce("Time<9.8||Time>10.5");
    //generating the tail of the distribution
    x.setRange("tail",9.5,11) ;
    RooAbsReal* igx_1 = model_1.getNormIntegral(RooArgSet(x));//createIntegral(x,Range("tail")) ;
    cout << "gx_Int[x] = " << igx_1->getVal() << "    "<<model_1.expectedEvents(0)<<"   "<<model_1.getValV(0)<<endl ;
    
    RooPlot* xframe2_1_toy = x.frame(Title(Form("pos. 1, channel toy tail %d",ch))) ;
    ds_1_toy_core->plotOn(xframe2_1_toy);
    TCanvas *ggg = new TCanvas("ggg","ggg");
    xframe2_1_toy->Draw();
  }



  
  mean_L_0.setConstant(kFALSE);
  if(fix_params==0){ ///////everything is free
    Delta_H_0.setVal(starting_delta_H_0);
    Delta_H_1.setVal(starting_delta_H_1);
    mean_L_0.setVal(starting_mean_L_0);
    delta_means_L.setVal(starting_delta_means_L);
    //mean_L_1.setVal(starting_mean_L_1);
    alpha_0.setVal(starting_alpha_0);
    alpha_1.setVal(starting_alpha_1);
  }else if(fix_params==1){///////fix only separation between peaks
    Delta_H_0.setConstant(kTRUE);//
    Delta_H_1.setConstant(kTRUE);//
    mean_L_0.setVal(starting_mean_L_0);
    delta_means_L.setVal(starting_delta_means_L);
    //mean_L_1.setVal(starting_mean_L_1);
    alpha_0.setVal(starting_alpha_0);
    alpha_1.setVal(starting_alpha_1);
  }else if(fix_params==2){///////fix separation between peaks+separation between first peak positions ****RECOMENDED****
    Delta_H_0.setConstant(kTRUE);//
    Delta_H_1.setConstant(kTRUE);//
    delta_means_L.setConstant(kTRUE);
    mean_L_0.setVal(starting_mean_L_0);
    mean_L_0.setConstant(kFALSE);
    //mean_L_1.setConstant(kTRUE);
    alpha_0.setVal(starting_alpha_0);
    alpha_1.setVal(starting_alpha_1);
  } else if(fix_params==3){//fix absolute positions+separations+fractions
    Delta_H_0.setConstant(kTRUE);//
    Delta_H_1.setConstant(kTRUE);//
    mean_L_0.setConstant(kTRUE);
    delta_means_L.setConstant(kTRUE);
    //mean_L_1.setConstant(kTRUE);
    alpha_0.setConstant(kTRUE);
    alpha_1.setConstant(kTRUE);
  }
  
  alpha_CB.setVal(starting_alpha_CB);
  n_CB.setVal(starting_n_CB);
  if(do_PixByPix_CBparams_fit){  n_CB.setConstant(kTRUE);}// alpha_CB.setConstant(kTRUE);}

  sigma_L_0.setVal(starting_sigma_L_0);
  sigma_H_0.setVal(starting_sigma_H_0);
  sigma_L_1.setVal(starting_sigma_L_1);
  sigma_H_1.setVal(starting_sigma_H_1);
  
  
  if(do_prefit_fullSpectrum){
    
    cout<<"___________________________________"<<endl;
    cout<<"                                   "<<endl;
    cout<<"-----------------------------------"<<endl;
    cout<<"***********************************"<<endl;
    cout<<"-----------------------------------"<<endl;
    cout<<"___________________________________"<<endl;
    cout<<"STATUS MIGRAD EDM 10 0"<<endl;
    
    //alpha_CB.setConstant(kTRUE);
    //n_CB.setConstant(kTRUE);
    sigma_L_0.setConstant(kTRUE);
    sigma_H_0.setConstant(kTRUE);
    sigma_L_1.setConstant(kTRUE);
    sigma_H_1.setConstant(kTRUE);
    RooFitResult* fit_results_b = model_b.fitTo(DS,Save(),Minimizer(Type_minim_pf,Algo_minim_pf),SumW2Error(kFALSE),PrintLevel(MN_output_print_level_prefit),PrintEvalErrors(-1),Warnings(kFALSE),Verbose(kFALSE));//);
    sigma_L_1.setConstant(kFALSE);
    sigma_H_1.setConstant(kFALSE);
    sigma_L_0.setConstant(kFALSE);
    sigma_H_0.setConstant(kFALSE);
    // alpha_CB.setConstant(kFALSE);
    //n_CB.setConstant(kFALSE);
    
    
    
    if(print_prefit_info) fit_results_b->Print("v");
  }
  
  
  
  
  cout<<"___________________________________"<<endl;
  cout<<"                                   "<<endl;
  cout<<"-----------------------------------"<<endl;
  cout<<"***********************************"<<endl;
  cout<<"-----------------------------------"<<endl;
  cout<<"___________________________________"<<endl;
  cout<<"STATUS MIGRAD EDM 10"<<endl;
  
  
  
  
  RooFitResult* fit_results;
  if(use_NLL){  
    RooAbsReal* nll_model = model.createNLL(DS,Extended(kFALSE)) ;
    //RooMinimizer RMN = RooMinimizer(*nll_model);
    RooMinimizer RMN (*nll_model);
    RMN.setErrorLevel(-1);
    RMN.setVerbose(kFALSE);
    RMN.setPrintEvalErrors(-1);
    RMN.setPrintLevel(MN_output_print_level);
    RMN.setStrategy(2);
    RMN.minimize(Type_minim,Algo_minim);
    fit_results=RMN.fit("hmrt") ;
  }else{
    fit_results = model.fitTo(DS,Save(),InitialHesse(true),Strategy(2),SumW2Error(kFALSE),PrintLevel(MN_output_print_level),PrintEvalErrors(-1),Warnings(kFALSE));//,InitialHesse(true));//,Hesse(kFALSE));//Minimizer(Type_minim,Algo_minim)
  }
  /*
    alpha_0.setVal(starting_alpha_0);
    alpha_1.setVal(starting_alpha_1);
    mean_L_0.setVal(starting_mean_L_0);
    mean_L_1.setVal(starting_mean_L_1);
    Delta_H_0.setVal(starting_delta_H_0);
    Delta_H_1.setVal(starting_delta_H_1);
    
    
    
    RooAbsReal* nll_model = model.createNLL(DS,Save(),Range("Fit_Range"),SumW2Error(kFALSE),Extended(kFALSE)) ;
    RooMinuit RMN = RooMinuit(*nll_model);
    RMN.setWarnLevel(0);
    RMN.setVerbose(kFALSE);
    mean_L_0.setConstant(kTRUE);
    Delta_H_0.setConstant(kTRUE);
    sigma_L_0.setConstant(kTRUE);
    sigma_H_0.setConstant(kTRUE);
    mean_L_1.setConstant(kTRUE);
    Delta_H_1.setConstant(kTRUE);
    sigma_L_1.setConstant(kTRUE);
    sigma_H_1.setConstant(kTRUE);
    RMN.migrad() ;
    mean_L_0.setConstant(kFALSE);
    Delta_H_0.setConstant(kFALSE);
    sigma_L_0.setConstant(kFALSE);
    sigma_H_0.setConstant(kFALSE);
    mean_L_1.setConstant(kFALSE);
    Delta_H_1.setConstant(kFALSE);
    sigma_L_1.setConstant(kFALSE);
    sigma_H_1.setConstant(kFALSE);
    RooFitResult* fit_results=RMN.fit("shr") ;
    
    
    
  */
  
  fit_results->Print("v");
  
  Fits_status.push_back(fit_results->status());
  
  RooPlot* xframe2 = x.frame(Title(Form("pos.0 #oplus pos.1,  channel %d",ch))) ;
  DS.plotOn(xframe2);
  TH2 *h_correlation = fit_results->correlationHist();  h_correlation->SetTitle("correlation matrix 0 #oplus 1");h_correlation->SetTitle(Form("correlation matrix 0 #oplus 1 ch %d",ch));h_correlation->GetYaxis()->SetLabelSize(0.1); h_correlation->GetYaxis()->SetLabelFont(70); h_correlation->SetMarkerSize(2);
  model.plotOn(xframe2,Range("Fit_Range"),Components("PDF_L_0"),LineStyle(kDashed),LineColor(1)) ;
  model.plotOn(xframe2,Range("Fit_Range"),Components("PDF_H_0"),LineStyle(kDashed),LineColor(1)) ;
  //  model.plotOn(xframe2,Range("Fit_Range"),Components("PDF_T_0"),LineStyle(kDashed),LineColor(1)) ;
  model.plotOn(xframe2,Range("Fit_Range"),Components("PDF_L_1"),LineStyle(kDashed),LineColor(8)) ;
  model.plotOn(xframe2,Range("Fit_Range"),Components("PDF_H_1"),LineStyle(kDashed),LineColor(8)) ;
  //  model.plotOn(xframe2,Range("Fit_Range"),Components("PDF_T_1"),LineStyle(kDashed),LineColor(8)) ;
  model.plotOn(xframe2,Range("Fit_Range"),Components("PDF_bkg"),LineStyle(kDashed),LineColor(2)) ;
  model.plotOn(xframe2,Range("Fit_Range"),Components("PDF_sig"),LineStyle(kDashed),LineColor(5)) ;
  model.plotOn(xframe2,Range("Fit_Range"));
  // Construct a histogram with the residuals of the data w.r.t. the curve
  RooHist* hresid = xframe2->residHist() ;
  // Create a new frame to draw the pull distribution and add the distribution to the frame
  RooPlot* xframe3 = x.frame(Title("Pull Distribution")) ;
  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist* hpull = xframe2->pullHist() ;
  xframe3->addPlotable(hresid,"P") ;
  // Create a new frame to draw the pull distribution and add the distribution to the frame
  RooPlot* xframe4 = x.frame(Title(Form("Pull Distribution fiber 0 #oplus 1 - ch %d",ch)));
  xframe4->addPlotable(hpull,"P") ;
  xframe2->GetXaxis()->SetLabelSize(0.05);
  xframe2->GetXaxis()->SetLabelFont(70);
  xframe2->GetXaxis()->SetTitleSize(0.05);
  xframe2->GetXaxis()->SetTitleFont(70);
  xframe4->GetXaxis()->SetLabelSize(0.05);
  xframe4->GetXaxis()->SetLabelFont(70);
  xframe4->GetXaxis()->SetTitleSize(0.05);
  xframe4->GetXaxis()->SetTitleFont(70);
  
  
  RooPlot* xframe2_log = x.frame(Title(Form("channel %d",ch))) ;
  DS.plotOn(xframe2_log);//,DataError(RooAbsData::SumW2)) ;
  model.plotOn(xframe2_log,Components("PDF_L_0"),LineStyle(kDashed),LineColor(1)) ;
  model.plotOn(xframe2_log,Components("PDF_H_0"),LineStyle(kDashed),LineColor(1)) ;
  //  model.plotOn(xframe2_log,Components("PDF_T_0"),LineStyle(kDashed),LineColor(1)) ;
  model.plotOn(xframe2_log,Components("PDF_L_1"),LineStyle(kDashed),LineColor(8)) ;
  model.plotOn(xframe2_log,Components("PDF_H_1"),LineStyle(kDashed),LineColor(8)) ;
  //  model.plotOn(xframe2_log,Components("PDF_T_1"),LineStyle(kDashed),LineColor(8)) ;
  model.plotOn(xframe2_log,Components("PDF_bkg"),LineStyle(kDashed),LineColor(2)) ;
  model.plotOn(xframe2_log,Components("PDF_sig"),LineStyle(kDashed),LineColor(5)) ;
  model.plotOn(xframe2_log) ;
  
  
  model.paramOn(xframe2_log);
  model.paramOn(xframe2,Layout(0.8,0.5,1));
  xframe2_log->GetXaxis()->SetLabelSize(0.05);
  xframe2_log->GetXaxis()->SetLabelFont(70);
  xframe2_log->GetXaxis()->SetTitleSize(0.05);
  xframe2_log->GetXaxis()->SetTitleFont(70);


  
  
  TF1 *line_0s = new TF1("line_0s","0",-1,16);line_0s->SetLineColor(8);
  TF1 *line_1ps = new TF1("line_1ps","1",-1,16);line_1ps->SetLineColor(4);
  TF1 *line_1ns = new TF1("line_1ns","-1",-1,16);line_1ns->SetLineColor(4);
  TF1 *line_2ps = new TF1("line_2ps","2",-1,16);line_2ps->SetLineColor(kOrange-2);
  TF1 *line_2ns = new TF1("line_2ns","-2",-1,16);line_2ns->SetLineColor(kOrange-2);
  TF1 *line_3ps = new TF1("line_3ps","3",-1,16);line_3ps->SetLineColor(2);
  TF1 *line_3ns = new TF1("line_3ns","-3",-1,16);line_3ns->SetLineColor(2);
  if(draw_results){  
    TCanvas* c_Fit = new TCanvas("Fit results","Fit results",0,0,1124,700) ;
    c_Fit->Divide(3,4) ;
    gStyle->SetOptFit(0111); 
    c_Fit->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe2_0->GetYaxis()->SetTitleOffset(1.6) ; xframe2_0->Draw() ;
    c_Fit->cd(4) ; gPad->SetLeftMargin(0.15) ; xframe2_0_log->GetYaxis()->SetTitleOffset(1.6) ; gPad->SetLogy(1); xframe2_0_log->Draw() ;
    c_Fit->cd(7) ; gPad->SetLeftMargin(0.15) ; xframe4_0->GetYaxis()->SetTitleOffset(1.6) ; gPad->SetGridy(); xframe4_0->GetYaxis()->SetRangeUser(-5,5); xframe4_0->Draw() ;
    line_3ns->Draw("same");line_3ps->Draw("same");
    line_2ns->Draw("same");line_2ps->Draw("same");
    line_1ns->Draw("same");line_1ps->Draw("same");
    line_0s->Draw("same");line_0s->Draw("same");
    c_Fit->cd(10) ; gPad->SetLeftMargin(0.15) ; h_correlation_0->Draw("colz:text");
    
    c_Fit->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe2_1->GetYaxis()->SetTitleOffset(1.6) ; xframe2_1->Draw() ;
    c_Fit->cd(5) ; gPad->SetLeftMargin(0.15) ; xframe2_1_log->GetYaxis()->SetTitleOffset(1.6) ; gPad->SetLogy(1); xframe2_1_log->Draw() ;
    //c_Fit->cd(5) ; gPad->SetLeftMargin(0.15) ; xframe2_1->GetYaxis()->SetTitleOffset(1.6) ; gPad->SetLogy(1); xframe2_1->Draw() ;
    c_Fit->cd(8) ; gPad->SetLeftMargin(0.15) ; xframe4_1->GetYaxis()->SetTitleOffset(1.6) ; gPad->SetGridy(); xframe4_1->GetYaxis()->SetRangeUser(-5,5); xframe4_1->Draw() ;
    line_3ns->Draw("same");line_3ps->Draw("same");
    line_2ns->Draw("same");line_2ps->Draw("same");
    line_1ns->Draw("same");line_1ps->Draw("same");
    line_0s->Draw("same");line_0s->Draw("same");
    c_Fit->cd(11) ; gPad->SetLeftMargin(0.15) ; h_correlation_1->Draw("colz:text");
    
    TCanvas* c_Fit_zoom = new TCanvas("c_Fit_zoom","Fit results zoom",0,0,1124,700) ;
    c_Fit_zoom->Divide(2,2) ;
    gStyle->SetOptFit(0111); 
    c_Fit_zoom->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe2_0->GetYaxis()->SetTitleOffset(1.6) ; xframe2_0->Draw() ;
    c_Fit_zoom->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe2_1->GetYaxis()->SetTitleOffset(1.6) ; xframe2_1->Draw() ;
    TCanvas* c_Fit_0 = new TCanvas("c_Fit_0","c_Fit_0");
     xframe2_0->Draw() ;
    TCanvas* c_Fit_1 = new TCanvas("c_Fit_1","c_Fit_1");
     xframe2_1->Draw() ;

    c_Fit->cd();
    c_Fit->cd(3) ; gPad->SetLeftMargin(0.15) ; xframe2->GetYaxis()->SetTitleOffset(1.6) ; xframe2->Draw() ;
    c_Fit->cd(6) ; gPad->SetLeftMargin(0.15) ; xframe2_log->GetYaxis()->SetTitleOffset(1.6) ; gPad->SetLogy(1); xframe2_log->Draw() ;
    c_Fit->cd(9) ; gPad->SetLeftMargin(0.15) ; xframe4->GetYaxis()->SetTitleOffset(1.6) ; gPad->SetGridy(); xframe4->GetYaxis()->SetRangeUser(-5,5);xframe4->Draw() ;
    line_3ns->Draw("same");line_3ps->Draw("same");
    line_2ns->Draw("same");line_2ps->Draw("same");
    line_1ns->Draw("same");line_1ps->Draw("same");
    line_0s->Draw("same");line_0s->Draw("same");
    c_Fit->cd(12) ; gPad->SetLeftMargin(0.15) ; h_correlation->Draw("colz:text");
    c_Fit_zoom->cd();
    c_Fit_zoom->cd(3) ; gPad->SetLeftMargin(0.15) ; xframe2->GetYaxis()->SetTitleOffset(1.6) ; xframe2->Draw() ;
    TCanvas* c_Fit_01 = new TCanvas("c_Fit_01","c_Fit_01");
     xframe2->Draw() ;
  }


  

  if(compute_FWHM){
    double Res;
    Double_t HW_CB;
    double absAlpha; 
    double n;
    if(absAlpha>=sqrt(2*log(2))) {
      HW_CB=sqrt(2*log(2));
    } else {
      absAlpha=TMath::Abs(alpha_CB.getVal());
      n= n_CB.getVal();
      Double_t A = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
      Double_t B= n/absAlpha - absAlpha;
      HW_CB = -B + TMath::Power(2*A,1/n);////////the proper formula is multiplied by -1
    }
    Res=sigma_L_0.getVal();
    sigma_L_0.setVal(Res*0.5*(sqrt(2*log(2))+HW_CB));
    Res=sigma_L_1.getVal();
    sigma_L_1.setVal(Res*0.5*(sqrt(2*log(2))+HW_CB));
    Res=sigma_H_0.getVal();
    sigma_H_0.setVal(Res*0.5*(sqrt(2*log(2))+HW_CB));
    Res=sigma_H_1.getVal();
    sigma_H_1.setVal(Res*0.5*(sqrt(2*log(2))+HW_CB));
  }
  
  ///////BOTH FIBER FIT VALUES////////
  /*32   */   POIs.push_back(mean_L_0.getVal()); 
  /*33   */   POIs.push_back(mean_L_0.getError());
  /*34   */   POIs.push_back(Delta_H_0.getVal());
  /*35   */   POIs.push_back(Delta_H_0.getError());
  /*36   */   POIs.push_back(Delta_T_0.getVal());
  /*37   */   POIs.push_back(Delta_T_0.getError());
  /*38   */   POIs.push_back(sigma_L_0.getVal());
  /*39   */   POIs.push_back(sigma_L_0.getError());
  /*40   */   POIs.push_back(sigma_H_0.getVal());
  /*41   */   POIs.push_back(sigma_H_0.getError());
  /*42   */   POIs.push_back(sigma_T_0.getVal());
  /*43   */   POIs.push_back(sigma_T_0.getError());
  /*44   */   POIs.push_back(alpha_0.getVal());
  /*45   */   POIs.push_back(alpha_0.getError());
  /*46   */   POIs.push_back(beta_0.getVal());
  /*47   */   POIs.push_back(beta_0.getError());
  
  /*48   */  POIs.push_back(delta_means_L.getVal()); //mean_L_1.getVal()); 
  /*49   */  POIs.push_back(delta_means_L.getError());//mean_L_1.getVal()); //getError());
  /*50   */  POIs.push_back(Delta_H_1.getVal());
  /*51   */  POIs.push_back(Delta_H_1.getError());
  /*52   */  POIs.push_back(Delta_T_1.getVal());
  /*53   */  POIs.push_back(Delta_T_1.getError());
  /*54   */  POIs.push_back(sigma_L_1.getVal());
  /*55   */  POIs.push_back(sigma_L_1.getError());
  /*56   */  POIs.push_back(sigma_H_1.getVal());
  /*57   */  POIs.push_back(sigma_H_1.getError());
  /*58  */   POIs.push_back(sigma_T_1.getVal());
  /*59  */   POIs.push_back(sigma_T_1.getError());
  /*60  */   POIs.push_back(alpha_1.getVal());
  /*61  */   POIs.push_back(alpha_1.getError());
  /*62  */   POIs.push_back(beta_1.getVal());
  /*63  */   POIs.push_back(beta_1.getError());
  
  /*64 */   POIs.push_back(Frac_0.getVal());
  /*65 */   POIs.push_back(Frac_0.getError());
  /*66 */   POIs.push_back(Frac_sig.getVal());
  /*67 */   POIs.push_back(Frac_sig.getError());
  
  cout<<"   "<<endl;
  cout<<"***************************************************************************"<<endl;
  cout<<"***************************************************************************"<<endl;
  cout<<"<<<<< FULL SPECTRA FIT >>>>>"<<endl;
  cout<<" #### Fit quality checks, channel "<<ch<<" , fiber position 0+1 ####"<<endl;
  
  cout<<"   "<<endl;
  cout<<"***************************************************************************"<<endl;
  
  cout<<"------- General quality parameters ---------"<<endl;
  cout<<"Chi2/ndof  : "<<  xframe2->chiSquare()<<endl;
  cout<<"EDM : "<<fit_results->edm()<<endl;
  cout<<"Fit status : "<<fit_results->status()<<endl;
  cout<<"--------------------------------------------"<<endl;
  cout<<" #### Fit quality checks, channel "<<ch<<" , fiber position spectrum 0 ####"<<endl;
  cout<<"-------------  Hit boundaries --------------"<<endl;
  if(((POIs[34+2]-low_delta_H_0)<0.0001) || ((-POIs[34+2]+up_delta_H_0)<0.001))   cout<<"delta H_0 hits boundary"<<endl;
  //if(((POIs[34+4]-low_delta_T_0)<0.0001) || ((-POIs[34+4]+up_delta_T_0)<0.001))   cout<<"delta T_0 hits boundary"<<endl;
  if(((POIs[34+6]-low_sigma_L_0)<0.0001) || ((-POIs[34+6]+up_sigma_L_0)<0.001))   cout<<"sigma L_0 hits boundary"<<endl;
  if(((POIs[34+8]-low_sigma_H_0)<0.0001) || ((-POIs[34+8]+up_sigma_H_0)<0.001))   cout<<"sigma H_0 hits boundary"<<endl;
  //if(((POIs[34+10]-low_sigma_T_0)<0.0001)|| ((-POIs[34+10]+up_sigma_T_0)<0.001))  cout<<"sigma T_0 hits boundary"<<endl;
  if(((POIs[34+12]-low_alpha_0)<0.0001)  || ((-POIs[34+12]+up_alpha_0)<0.001))    cout<<"alpha_0  hits boundary"<<endl;
  //if(((POIs[34+14]-low_beta_0)<0.0001)   || ((-POIs[34+14]+up_beta_0)<0.001))     cout<<"beta_0 hits boundary"<<endl;
  cout<<"--------------- Large/small errors --------------"<<endl;
  for(int h=0;h<8;h++){
    if(h!=2&&h!=5&&h!=7){
      if(TMath::Abs(POIs[32+2*h+1]/POIs[32+2*h])>0.50)    cout<<"large error at "<<PAR_NAMES[h]<<"_0"<<endl;
      if(TMath::Abs(POIs[32+2*h+1]/POIs[32+2*h])<0.00001) cout<<"small error at "<<PAR_NAMES[h]<<"_0"<<endl;
    }
  }
  cout<<"   "<<endl;
  cout<<"***************************************************************************"<<endl;
  cout<<" #### Fit quality checks, channel "<<ch<<" , fiber position spectrum 1 ####"<<endl;
  
  cout<<"-------------  Hit boundaries --------------"<<endl;
  if(((POIs[50+2]-low_delta_H_1)<0.0001) || ((-POIs[50+2]+up_delta_H_1)<0.001))   cout<<"delta H_1 hits boundary"<<endl;
  //  if(((POIs[50+4]-low_delta_T_1)<0.0001) || ((-POIs[50+4]+up_delta_T_1)<0.001))   cout<<"delta T_1 hits boundary"<<endl;
  if(((POIs[50+6]-low_sigma_L_1)<0.0001) || ((-POIs[50+6]+up_sigma_L_1)<0.001))   cout<<"sigma L_1 hits boundary"<<endl;
  if(((POIs[50+8]-low_sigma_H_1)<0.0001) || ((-POIs[50+8]+up_sigma_H_1)<0.001))   cout<<"sigma H_1 hits boundary"<<endl;
  //  if(((POIs[50+10]-low_sigma_T_1)<0.0001)|| ((-POIs[50+10]+up_sigma_T_1)<0.001))  cout<<"sigma T_1 hits boundary"<<endl;
  if(((POIs[50+12]-low_alpha_1)<0.0001)  || ((-POIs[50+12]+up_alpha_1)<0.001))    cout<<"alpha_1  hits boundary"<<endl;
  //if(((POIs[50+14]-low_beta_1)<0.0001)   || ((-POIs[50+14]+up_beta_1)<0.001))     cout<<"beta_1 hits boundary"<<endl;
  cout<<"--------------- Large/small errors --------------"<<endl;
  for(int h=0;h<8;h++){
    if(h!=2&&h!=5&&h!=7){
      if(TMath::Abs(POIs[48+2*h+1]/POIs[48+2*h])>0.50)    cout<<"large error at "<<PAR_NAMES[h]<<"_1"<<endl;
      if(TMath::Abs(POIs[48+2*h+1]/POIs[48+2*h])<0.00001) cout<<"small error at "<<PAR_NAMES[h]<<"_1"<<endl;
    }
  }
  
  
  cout<<"--------------- Correlation matrix --------------"<<endl;
  
  
  for(int ff=1;ff<h_correlation->GetXaxis()->GetNbins()+1;ff++){
    for(int gg=1;gg<h_correlation->GetXaxis()->GetNbins()-ff+1;gg++){
      if(TMath::Abs(h_correlation->GetBinContent(ff,gg))<0.001||TMath::Abs(h_correlation->GetBinContent(ff,gg))>0.99){
	cout<<"low/high "<<h_correlation->GetXaxis()->GetBinLabel(ff)<<"-"<<h_correlation->GetYaxis()->GetBinLabel(gg)<<" correlation : "<<h_correlation->GetBinContent(ff,gg)<< endl;
      }
    }
  }
  cout<<"<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
  cout<<"<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
  cout<<"<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
  
  
  
  
  //my_fit_results.Results=POIs;
  for(int g=0; g<POIs.size();g++){
    my_fit_results.Results[g]=POIs[g];
  }
  my_fit_results.n_POIs=POIs.size();
  my_fit_results.xframe2_fit_0=xframe2_0;
  my_fit_results.xframe2_pull_0=xframe4_0;
  my_fit_results.h_correlation_0=h_correlation_0;
  my_fit_results.xframe2_fit_1=xframe2_1;
  my_fit_results.xframe2_pull_1=xframe4_1;
  my_fit_results.h_correlation_1=h_correlation_1;
  my_fit_results.xframe2_fit_01=xframe2;
  my_fit_results.xframe2_pull_01=xframe4;
  my_fit_results.h_correlation_01=h_correlation;
  
  //cout<<my_fit_results.Results[0]<<endl;
  
  
  //  cout<<"test"<<endl;
  //  delete  h_input_histogram_0; //cout<<"test 1"<<endl;
  //  delete  h_input_histogram_1; //cout<<"test 2"<<endl;
  
  
  return my_fit_results;//POIs;
  
  //return POIs;
}



vector<float> loop_channels(int deep_fixed_params,bool plot_summaries){ //rel_weight is the relative weight of the 1 to 0 pois, i.e. total dataset= pos.0 + rel_weight*pos.1 - in the real case rel_weight=1.0
  
  
  TCanvas *c_pos0_AllChannels = new TCanvas("c_pos0_AllChannels","c_pos0_AllChannels",0,0,1124,700);
  TCanvas *c_pos0_AllChannels_pulls = new TCanvas("c_pos0_AllChannels_pulls","c_pos0_AllChannels_pulls",0,0,1124,700);
  TCanvas *c_pos0_AllChannels_corr = new TCanvas("c_pos0_AllChannels_corr","c_pos0_AllChannels_corr",0,0,1124,700);
  TCanvas *c_pos1_AllChannels = new TCanvas("c_pos1_AllChannels","c_pos1_AllChannels",0,0,1124,700);
  TCanvas *c_pos1_AllChannels_pulls = new TCanvas("c_pos1_AllChannels_pulls","c_pos1_AllChannels_pulls",0,0,1124,700);
  TCanvas *c_pos1_AllChannels_corr = new TCanvas("c_pos1_AllChannels_corr","c_pos1_AllChannels_corr",0,0,1124,700);
  TCanvas *c_pos01_AllChannels = new TCanvas("c_pos01_AllChannels","c_pos01_AllChannels",0,0,1124,700);
  TCanvas *c_pos01_AllChannels_pulls = new TCanvas("c_pos01_AllChannels_pulls","c_pos01_AllChannels_pulls",0,0,1124,700);
  TCanvas *c_pos01_AllChannels_corr = new TCanvas("c_pos01_AllChannels_corr","c_pos01_AllChannels_corr",0,0,1124,700);
  
  
  
  vector<float> relative_weight_Frac_0; relative_weight_Frac_0.clear();
  c_pos0_AllChannels->Divide(4,4);
  c_pos1_AllChannels->Divide(4,4);
  c_pos01_AllChannels->Divide(4,4);
  c_pos0_AllChannels_pulls->Divide(4,4);
  c_pos1_AllChannels_pulls->Divide(4,4);
  c_pos01_AllChannels_pulls->Divide(4,4);
  c_pos0_AllChannels_corr->Divide(4,4);
  c_pos1_AllChannels_corr->Divide(4,4);
  c_pos01_AllChannels_corr->Divide(4,4);
  
  
  
  
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
  
  //vector<Fit_results> RESULTS;
  //RESULTS.clear();
  Fits_status.clear();
  for(Int_t i=0; i<16; i++){
    x[i]=i; err_x[i]=0;
    POIs.clear();
    cout<<"Channel : "<<i<<endl;
    Fit_results my_fit_Results;
    
    my_fit_Results=Fit_head("blind",deep_fixed_params,i);
    //POIs=my_fit_Results.Results;//Fit_head(i,fix_params,rel_weight,false);
    for(int hg=0;hg<my_fit_Results.n_POIs;hg++){
      POIs.push_back(my_fit_Results.Results[hg]);
    }//Fit_head(i,fix_params,rel_weight,false);
    //RESULTS.push_back(my_fit_Results);
    
    
    TF1 *line_0s = new TF1("line_0s","0",-1,16);line_0s->SetLineColor(8);
    TF1 *line_1ps = new TF1("line_1ps","1",-1,16);line_1ps->SetLineColor(4);
    TF1 *line_1ns = new TF1("line_1ns","-1",-1,16);line_1ns->SetLineColor(4);
    TF1 *line_2ps = new TF1("line_2ps","2",-1,16);line_2ps->SetLineColor(kOrange-2);
    TF1 *line_2ns = new TF1("line_2ns","-2",-1,16);line_2ns->SetLineColor(kOrange-2);
    TF1 *line_3ps = new TF1("line_3ps","3",-1,16);line_3ps->SetLineColor(2);
    TF1 *line_3ns = new TF1("line_3ns","-3",-1,16);line_3ns->SetLineColor(2);
    c_pos0_AllChannels->cd(index_channel_pixel[i]); my_fit_Results.xframe2_fit_0->Draw(); 
    c_pos1_AllChannels->cd(index_channel_pixel[i]);my_fit_Results.xframe2_fit_1->Draw();
    c_pos0_AllChannels_pulls->cd(index_channel_pixel[i]); gPad->SetGridy(); my_fit_Results.xframe2_pull_0->GetYaxis()->SetRangeUser(-5,5); my_fit_Results.xframe2_pull_0->Draw();
    line_3ns->Draw("same");line_3ps->Draw("same");
    line_2ns->Draw("same");line_2ps->Draw("same");
    line_1ns->Draw("same");line_1ps->Draw("same");
    line_0s->Draw("same");line_0s->Draw("same");
    c_pos1_AllChannels_pulls->cd(index_channel_pixel[i]); gPad->SetGridy(); my_fit_Results.xframe2_pull_1->GetYaxis()->SetRangeUser(-5,5);my_fit_Results.xframe2_pull_1->Draw();
    line_3ns->Draw("same");line_3ps->Draw("same");
    line_2ns->Draw("same");line_2ps->Draw("same");
    line_1ns->Draw("same");line_1ps->Draw("same");
    line_0s->Draw("same");line_0s->Draw("same");
    c_pos0_AllChannels_corr->cd(index_channel_pixel[i]);my_fit_Results.h_correlation_0->Draw("colz:text");
    c_pos1_AllChannels_corr->cd(index_channel_pixel[i]);my_fit_Results.h_correlation_1->Draw("colz:text");
    c_pos01_AllChannels->cd(index_channel_pixel[i]);my_fit_Results.xframe2_fit_01->Draw(); 
    c_pos01_AllChannels_pulls->cd(index_channel_pixel[i]); gPad->SetGridy(); my_fit_Results.xframe2_pull_01->GetYaxis()->SetRangeUser(-5,5); my_fit_Results.xframe2_pull_01->Draw(); 
    line_3ns->Draw("same");line_3ps->Draw("same");
    line_2ns->Draw("same");line_2ps->Draw("same");
    line_1ns->Draw("same");line_1ps->Draw("same");
    line_0s->Draw("same");line_0s->Draw("same");
    c_pos01_AllChannels_corr->cd(index_channel_pixel[i]);my_fit_Results.h_correlation_01->Draw("colz:text"); 
    
    
    
    
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
    
    mean_L_1_B[i]      =POIs[16];
    err_mean_L_1_B[i]  =POIs[17];
    sigma_L_1_B[i]     =POIs[22];
    err_sigma_L_1_B[i] =POIs[23];
    Delta_H_1_B[i]     =POIs[18];
    err_Delta_H_1_B[i] =POIs[19];
    sigma_H_1_B[i]     =POIs[24];
    err_sigma_H_1_B[i] =POIs[25];
    Delta_T_1_B[i]     =POIs[20];
    err_Delta_T_1_B[i] =POIs[21];
    sigma_T_1_B[i]     =POIs[26];
    err_sigma_T_1_B[i] =POIs[27];
    
    
    mean_L_0_A[i]      =POIs[32];
    err_mean_L_0_A[i]  =POIs[33];
    sigma_L_0_A[i]     =POIs[38];
    err_sigma_L_0_A[i] =POIs[39];
    Delta_H_0_A[i]     =POIs[34];
    err_Delta_H_0_A[i] =POIs[35];
    sigma_H_0_A[i]     =POIs[40];
    err_sigma_H_0_A[i] =POIs[41];
    Delta_T_0_A[i]     =POIs[36];
    err_Delta_T_0_A[i] =POIs[37];
    sigma_T_0_A[i]     =POIs[42];
    err_sigma_T_0_A[i] =POIs[43];
    
    mean_L_1_A[i]      =POIs[48];
    err_mean_L_1_A[i]  =POIs[49];
    sigma_L_1_A[i]     =POIs[54];
    err_sigma_L_1_A[i] =POIs[55];
    Delta_H_1_A[i]     =POIs[50];
    err_Delta_H_1_A[i] =POIs[51];
    sigma_H_1_A[i]     =POIs[56];
    err_sigma_H_1_A[i] =POIs[57];
    Delta_T_1_A[i]     =POIs[52];
    err_Delta_T_1_A[i] =POIs[53];
    sigma_T_1_A[i]     =POIs[58];
    err_sigma_T_1_A[i] =POIs[59];
    
    frac_L_B_0[i]         =POIs[12];
    err_frac_L_B_0[i]     =POIs[13];
    frac_L_B_1[i]         =POIs[28];
    err_frac_L_B_1[i]     =POIs[29];
    
    frac_L_A_0[i]         =POIs[44];
    err_frac_L_A_0[i]     =POIs[45];
    frac_L_A_1[i]         =POIs[60];
    err_frac_L_A_1[i]     =POIs[61];
    
    
    frac_H_B_0[i]         =POIs[14];
    err_frac_H_B_0[i]     =POIs[15];
    frac_H_B_1[i]         =POIs[30];
    err_frac_H_B_1[i]     =POIs[31];
    
    frac_H_A_0[i]         =POIs[46];
    err_frac_H_A_0[i]     =POIs[47];
    frac_H_A_1[i]         =POIs[62];
    err_frac_H_A_1[i]     =POIs[63];
    
    frac_0[i]             =POIs[64]; 
    err_frac_0[i]         =POIs[65];
    
    relative_weight_Frac_0.push_back(frac_0[i]);
    relative_weight_Frac_0.push_back(err_frac_0[i]);
    
    delta_ref_sigma_H[i]=ref_sigma_H[i]-sigma_H_0_B[i];
    delta_ref_sigma_L[i]=ref_sigma_L[i]-sigma_L_0_B[i];
  }
  
  
  if(plot_summaries){
    
    
    
    
    // if (relative_weight_Frac_0.size()!=32){cout<<"Warning, problem with the size of relative_weight_Frac_0"<<endl; return relative_weight_Frac_0;}
    TF1 *line = new TF1("line","0.100",-1,16);
    line->SetTitle("100 ps");
    // TLine *line = new TLine(0,100,16,100);
    line->SetLineColor(2);
    TGraphErrors *DELTA_means_ref = new TGraphErrors(16,x,ref_delta_mean,err_x,err_x);DELTA_means_ref->SetTitle("#Delta t_{ref} pos.0 - fit pos.0");
    TGraphErrors *DELTA_means_pos0_MC_B2 = new TGraphErrors(16,x,ref_pos0_MC_B2,err_x,err_x);DELTA_means_pos0_MC_B2->SetTitle("#Delta t_{ref} pos.0 - B2 MC");
    TGraphErrors *DELTA_means_pos1_MC_B2 = new TGraphErrors(16,x,ref_pos1_MC_B2,err_x,err_x);DELTA_means_pos1_MC_B2->SetTitle("#Delta t_{ref} pos.1 - B2 MC");
    TGraphErrors *DELTA_means_pos0_MC_PD = new TGraphErrors(16,x,ref_pos0_MC_PD,err_x,err_x);DELTA_means_pos0_MC_PD->SetTitle("#Delta t_{ref} pos.0 - PD MC");
    
    TGraphErrors *DELTA_means_pos0_3G = new TGraphErrors(16,x,ref_pos0_3G,err_x,err_ref_pos0_3G);DELTA_means_pos0_3G->SetTitle("#Delta t_{ref} pos.0 - 3 Gaussian model");
    TGraphErrors *DELTA_means_pos1_3G = new TGraphErrors(16,x,ref_pos1_3G,err_x,err_ref_pos1_3G);DELTA_means_pos1_3G->SetTitle("#Delta t_{ref} pos.1 - 3 Gaussian model");
    
    
    
    
    
    TGraphErrors *MEAN_L_B_0 = new TGraphErrors(16,x,mean_L_0_B,err_x,err_mean_L_0_B);MEAN_L_B_0->SetTitle("T_{L}^{pos.0} - fit pos.0");
    TGraphErrors *MEAN_L_B_1 = new TGraphErrors(16,x,mean_L_1_B,err_x,err_mean_L_1_B);MEAN_L_B_1->SetTitle("#Delta T_{L} #equiv T_{L}^{pos.1}-T_{L}^{pos.0} T- fit pos.1");
    TGraphErrors *MEAN_L_A_0 = new TGraphErrors(16,x,mean_L_0_A,err_x,err_mean_L_0_A);MEAN_L_A_0->SetTitle("T_{L}^{pos.0} - fit pos.0 #oplus 1");
    TGraphErrors *MEAN_L_A_1 = new TGraphErrors(16,x,mean_L_1_A,err_x,err_mean_L_1_A);MEAN_L_A_1->SetTitle("#Delta T_{L} #equiv T_{L}^{pos.1}-T_{L}^{pos.0} - fit pos.0 #oplus 1");
    
    TGraphErrors *DELTA_H_B_0 = new TGraphErrors(16,x,Delta_H_0_B,err_x,err_Delta_H_0_B);DELTA_H_B_0->SetTitle("#Delta t pos.0 - fit pos.0");
    TGraphErrors *DELTA_H_B_1 = new TGraphErrors(16,x,Delta_H_1_B,err_x,err_Delta_H_1_B);DELTA_H_B_1->SetTitle("#Delta t pos.1 - fit pos.1");
    TGraphErrors *DELTA_H_A_0 = new TGraphErrors(16,x,Delta_H_0_A,err_x,err_Delta_H_0_A);DELTA_H_A_0->SetTitle("#Delta t_{LOW-HIGH} pos.0 - fit pos.0 #oplus 1");
    TGraphErrors *DELTA_H_A_1 = new TGraphErrors(16,x,Delta_H_1_A,err_x,err_Delta_H_1_A);DELTA_H_A_1->SetTitle("#Delta t_{LOW-HIGH} pos.1 - fit pos.0 #oplus 1");
    TGraphErrors *DELTA_T_B_0 = new TGraphErrors(16,x,Delta_T_0_B,err_x,err_Delta_T_0_B);DELTA_T_B_0->SetTitle("#Delta t_{HIGH-THIRD} pos.0 - fit pos.0");
    TGraphErrors *DELTA_T_B_1 = new TGraphErrors(16,x,Delta_T_1_B,err_x,err_Delta_T_1_B);DELTA_T_B_1->SetTitle("#Delta t_{HIGH-THIRD} pos.1 - fit pos.1");
    TGraphErrors *DELTA_T_A_0 = new TGraphErrors(16,x,Delta_T_0_A,err_x,err_Delta_T_0_A);DELTA_T_A_0->SetTitle("#Delta t_{HIGH-THIRD} pos.0 - fit pos.0 #oplus 1");
    TGraphErrors *DELTA_T_A_1 = new TGraphErrors(16,x,Delta_T_1_A,err_x,err_Delta_T_1_A);DELTA_T_A_1->SetTitle("#Delta t_{HIGH-THIRD} pos.1 - fit pos.0 #oplus 1");
    
    TGraphErrors *SIGMA_L_B_0 = new TGraphErrors(16,x,sigma_L_0_B,err_x,err_sigma_L_0_B);SIGMA_L_B_0->SetTitle("#delta t_{low} pos.0 - fit pos.0");
    TGraphErrors *SIGMA_L_B_1 = new TGraphErrors(16,x,sigma_L_1_B,err_x,err_sigma_L_1_B);SIGMA_L_B_1->SetTitle("#delta t_{low} pos.1 - fit pos.1");
    TGraphErrors *SIGMA_H_B_0 = new TGraphErrors(16,x,sigma_H_0_B,err_x,err_sigma_H_0_B);SIGMA_H_B_0->SetTitle("#delta t_{high} pos.0 - fit pos.0");
    TGraphErrors *SIGMA_H_B_1 = new TGraphErrors(16,x,sigma_H_1_B,err_x,err_sigma_H_1_B);SIGMA_H_B_1->SetTitle("#delta t_{high} pos.1 - fit pos.1");
    TGraphErrors *SIGMA_T_B_0 = new TGraphErrors(16,x,sigma_T_0_B,err_x,err_sigma_T_0_B);SIGMA_T_B_0->SetTitle("#delta t_{third} pos.0 - fit pos.0");
    TGraphErrors *SIGMA_T_B_1 = new TGraphErrors(16,x,sigma_T_1_B,err_x,err_sigma_T_1_B);SIGMA_T_B_1->SetTitle("#delta t_{third} pos.1 - fit pos.1");
    
    TGraphErrors *SIGMA_L_A_0 = new TGraphErrors(16,x,sigma_L_0_A,err_x,err_sigma_L_0_A);SIGMA_L_A_0->SetTitle("#delta t_{low} pos.0 - fit pos.0 #oplus 1");
    TGraphErrors *SIGMA_L_A_1 = new TGraphErrors(16,x,sigma_L_1_A,err_x,err_sigma_L_1_A);SIGMA_L_A_1->SetTitle("#delta t_{low} pos.1 - fit pos.0 #oplus 1");
    TGraphErrors *SIGMA_H_A_0 = new TGraphErrors(16,x,sigma_H_0_A,err_x,err_sigma_H_0_A);SIGMA_H_A_0->SetTitle("#delta t_{high} pos.0 - fit pos.0 #oplus 1");
    TGraphErrors *SIGMA_H_A_1 = new TGraphErrors(16,x,sigma_H_1_A,err_x,err_sigma_H_1_A);SIGMA_H_A_1->SetTitle("#delta t_{high} pos.1 - fit pos.0 #oplus 1");
    TGraphErrors *SIGMA_T_A_0 = new TGraphErrors(16,x,sigma_T_0_A,err_x,err_sigma_T_0_A);SIGMA_T_A_0->SetTitle("#delta t_{third} pos.0 - fit pos.0 #oplus 1");
    TGraphErrors *SIGMA_T_A_1 = new TGraphErrors(16,x,sigma_T_1_A,err_x,err_sigma_T_1_A);SIGMA_T_A_1->SetTitle("#delta t_{third} pos.1 - fit pos.0 #oplus 1");
    
    
    TGraphErrors *FRAC_L_B_0 = new TGraphErrors(16,x,frac_L_B_0,err_x,err_frac_L_B_0);FRAC_L_B_0->SetTitle("f_{L} pos.0 - fit pos 0");
    TGraphErrors *FRAC_H_B_0 = new TGraphErrors(16,x,frac_H_B_0,err_x,err_frac_H_B_0);FRAC_H_B_0->SetTitle("f_{H} pos.0 - fit pos 0");
    TGraphErrors *FRAC_L_B_1 = new TGraphErrors(16,x,frac_L_B_1,err_x,err_frac_L_B_1);FRAC_L_B_1->SetTitle("f_{L} pos.1 - fit pos 1");
    TGraphErrors *FRAC_H_B_1 = new TGraphErrors(16,x,frac_H_B_1,err_x,err_frac_H_B_1);FRAC_H_B_1->SetTitle("f_{H} pos.1 - fit pos 1");
    TGraphErrors *FRAC_L_A_0 = new TGraphErrors(16,x,frac_L_A_0,err_x,err_frac_L_A_0);FRAC_L_A_0->SetTitle("f_{L} pos.0 - fit pos 0 #oplus 1");
    TGraphErrors *FRAC_H_A_0 = new TGraphErrors(16,x,frac_H_A_0,err_x,err_frac_H_A_0);FRAC_H_A_0->SetTitle("f_{H} pos.0 - fit pos 0 #oplus 1");
    TGraphErrors *FRAC_L_A_1 = new TGraphErrors(16,x,frac_L_A_1,err_x,err_frac_L_A_1);FRAC_L_A_1->SetTitle("f_{L} pos.1 - fit pos 0 #oplus 1");
    TGraphErrors *FRAC_H_A_1 = new TGraphErrors(16,x,frac_H_A_1,err_x,err_frac_H_A_1);FRAC_H_A_1->SetTitle("f_{H} pos.1 - fit pos 0 #oplus 1");
    
    
    TGraphErrors *FRAC_0 = new TGraphErrors(16,x,frac_0,err_x,err_frac_0);FRAC_0->SetTitle("F_{0}fit pos 0 #oplus 1");
    
    TGraphErrors *REF_sigma_L = new TGraphErrors(16,x,ref_sigma_L,err_x,err_x);REF_sigma_L->SetTitle("#delta t_{low}^{ref} pos.0 - fit pos.0");
    TGraphErrors *REF_sigma_H = new TGraphErrors(16,x,ref_sigma_H,err_x,err_x);REF_sigma_H->SetTitle("#delta t_{high}^{ref} pos.0 - fit pos.0");
    TGraphErrors *REF_sigma_T = new TGraphErrors(16,x,ref_sigma_T,err_x,err_x);REF_sigma_T->SetTitle("#delta t_{third}^{ref} pos.0 - fit pos.0");

    TGraphErrors *SIGMA_FRAC_L_0 = new TGraphErrors(16,frac_L_B_0,sigma_L_0_B,err_frac_L_B_0,err_sigma_L_0_B);SIGMA_FRAC_L_0->SetTitle("#delta t [ns] vs f_{L} - pos.0");
    TGraphErrors *SIGMA_FRAC_H_0 = new TGraphErrors(16,frac_L_B_0,sigma_H_0_B,err_frac_L_B_0,err_sigma_H_0_B);SIGMA_FRAC_H_0->SetTitle("#delta t [ns] vs f_{L} - pos.0");
    TGraphErrors *SIGMA_FRAC_L_1 = new TGraphErrors(16,frac_L_B_1,sigma_L_1_B,err_frac_L_B_1,err_sigma_L_1_B);SIGMA_FRAC_L_1->SetTitle("#delta t [ns] vs f_{L} - pos.1");
    TGraphErrors *SIGMA_FRAC_H_1 = new TGraphErrors(16,frac_L_B_1,sigma_H_1_B,err_frac_L_B_1,err_sigma_H_1_B);SIGMA_FRAC_H_1->SetTitle("#delta t [ns] vs f_{L} - pos.1");
    TGraphErrors *SIGMA_FRAC_L_0_A = new TGraphErrors(16,frac_L_B_0,sigma_L_0_A,err_frac_L_B_0,err_sigma_L_0_A);SIGMA_FRAC_L_0_A->SetTitle("#delta t [ns] vs f_{L} - pos.0");
    TGraphErrors *SIGMA_FRAC_H_0_A = new TGraphErrors(16,frac_L_B_0,sigma_H_0_A,err_frac_L_B_0,err_sigma_H_0_A);SIGMA_FRAC_H_0_A->SetTitle("#delta t [ns] vs f_{L} - pos.0");
    TGraphErrors *SIGMA_FRAC_L_1_A = new TGraphErrors(16,frac_L_B_1,sigma_L_1_A,err_frac_L_B_1,err_sigma_L_1_A);SIGMA_FRAC_L_1_A->SetTitle("#delta t [ns] vs f_{L} - pos.1");
    TGraphErrors *SIGMA_FRAC_H_1_A = new TGraphErrors(16,frac_L_B_1,sigma_H_1_A,err_frac_L_B_1,err_sigma_H_1_A);SIGMA_FRAC_H_1_A->SetTitle("#delta t [ns] vs f_{L} - pos.1");
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
    
    TMultiGraph *mg_MEANS_L_0 = new TMultiGraph();
    mg_MEANS_L_0->SetTitle("T_{L}^{pos.0} vs Channel - pos.0");
    mg_MEANS_L_0->Add(MEAN_L_B_0);
    mg_MEANS_L_0->Add(MEAN_L_A_0);
    TMultiGraph *mg_DELTAMEANS_L_1 = new TMultiGraph();
    mg_DELTAMEANS_L_1->SetTitle("#Delta T_{L} #equiv T_{L}^{pos.1}-T_{L}^{pos.0} - vs Channel - pos.1");
    mg_DELTAMEANS_L_1->Add(MEAN_L_B_1);
    //mg_DELTAMEANS_L_1->Add(MEAN_L_A_1);
    
    TCanvas *c_absolute_positions_01_b = new TCanvas("c_absolute_position_01_b","c_absolute_position_01_b");
    //c_absolute_positions_01->Divide(2,1);
    //c_absolute_positions_01->cd(1);
    mg_MEANS_L_0->Draw("AP");
    mg_MEANS_L_0->GetXaxis()->SetTitle("Channel");
    mg_MEANS_L_0->GetYaxis()->SetTitle("T_{L}^{pos.0} [ns]");
    gPad->BuildLegend();
    gPad->Update();
    gPad->SetGridy();
    TCanvas *c_absolute_positions_01 = new TCanvas("c_absolute_position_01","c_absolute_position_01");
    //c_absolute_positions_01->cd(2);
    mg_DELTAMEANS_L_1->Draw("AP");
    mg_DELTAMEANS_L_1->GetXaxis()->SetTitle("Channel");
    mg_DELTAMEANS_L_1->GetYaxis()->SetTitle("#Delta T_{L} [ns]");
    gPad->BuildLegend();
    gPad->Update();
    gPad->SetGridy();
    
    TMultiGraph *mg_DELTA_B_pos0 = new TMultiGraph();
    mg_DELTA_B_pos0->SetTitle("#Delta t (#equiv t_{H}-t_{L}) vs Channel - pos.0");
    //mg_DELTA_B->Add(DELTA_means_ref);
    mg_DELTA_B_pos0->Add(DELTA_means_pos0_3G);
    mg_DELTA_B_pos0->Add(DELTA_H_B_0);
    //mg_DELTA_B_pos0->Add(DELTA_means_pos0_MC_B2);
    //mg_DELTA_B_pos0->Add(DELTA_means_pos0_MC_PD);
    mg_DELTA_B_pos0->SetMinimum(0.1);
    mg_DELTA_B_pos0->SetMaximum(0.4);
    TCanvas *c_pos0_deltas_data_MC = new TCanvas("c_pos0_deltas_data_MC","c_pos0_deltas_data_MC");
    mg_DELTA_B_pos0->Draw("AP");
    gPad->BuildLegend();
    gPad->Update();
    gPad->SetGridy();
    mg_DELTA_B_pos0->GetXaxis()->SetTitle("Channel");
    mg_DELTA_B_pos0->GetYaxis()->SetTitle("#Delta t [ns]");
    mg_DELTA_B_pos0->GetYaxis()->SetNdivisions(20);
    gPad->Update();
    
    
    TMultiGraph *mg_DELTA_B_pos1 = new TMultiGraph();
    mg_DELTA_B_pos1->SetTitle("#Delta t (#equiv t_{H}-t_{L}) vs Channel - pos.1");
    mg_DELTA_B_pos1->Add(DELTA_means_pos1_3G);
    mg_DELTA_B_pos1->Add(DELTA_H_B_1);
    //    mg_DELTA_B_pos1->Add(DELTA_means_pos1_MC_B2);
    mg_DELTA_B_pos1->SetMinimum(0.1);
    mg_DELTA_B_pos1->SetMaximum(0.4);
    TCanvas *c_pos1_deltas_data_MC = new TCanvas("c_pos1_deltas_data_MC","c_pos1_deltas_data_MC");
    mg_DELTA_B_pos1->Draw("AP");
    gPad->BuildLegend();
    gPad->Update();
    gPad->SetGridy();
    mg_DELTA_B_pos1->GetXaxis()->SetTitle("Channel");
    mg_DELTA_B_pos1->GetYaxis()->SetTitle("#Delta t [ns]");
    mg_DELTA_B_pos1->GetYaxis()->SetNdivisions(20);
    gPad->Update();
    
    
    
    TMultiGraph *mg_SIGMA_B_L = new TMultiGraph();
    mg_SIGMA_B_L->SetTitle("Time resolution (#delta t) low time vs Channel");
    mg_SIGMA_B_L->Add(REF_sigma_L);
    mg_SIGMA_B_L->Add(SIGMA_L_B_0);
    mg_SIGMA_B_L->Add(SIGMA_L_A_0);
    TMultiGraph *mg_SIGMA_B_H = new TMultiGraph();
    mg_SIGMA_B_H->SetTitle("Time resolution (#delta t) high time vs Channel");
    mg_SIGMA_B_H->Add(REF_sigma_H);
    mg_SIGMA_B_H->Add(SIGMA_H_B_0);
    mg_SIGMA_B_H->Add(SIGMA_H_A_0);
    TMultiGraph *mg_SIGMA_B_T = new TMultiGraph();
    mg_SIGMA_B_T->SetTitle("Time resolution (#delta t) third time vs Channel");
    mg_SIGMA_B_T->Add(REF_sigma_T);
    mg_SIGMA_B_T->Add(SIGMA_T_B_0);
    mg_SIGMA_B_T->Add(SIGMA_T_A_0);
    mg_SIGMA_B_T->Add(SIGMA_T_B_1);
    mg_SIGMA_B_T->Add(SIGMA_T_A_1);
    
    
    
    
    ///////NEW MG_GRAPHS/////
    TMultiGraph *mg_DELTA_0 = new TMultiGraph();
    mg_DELTA_0->SetTitle("#Delta t (#equiv t_{H}-t_{L}) vs Channel - pos.0 ");
    //mg_DELTA_0->Add(DELTA_means_ref);
    mg_DELTA_0->Add(DELTA_H_B_0);
    //mg_DELTA_0->Add(DELTA_H_A_0);
    //mg_DELTA_0->Add(DELTA_T_B_0);
    //mg_DELTA_0->Add(DELTA_T_A_0);
    
    TMultiGraph *mg_DELTA_1 = new TMultiGraph();
    mg_DELTA_1->SetTitle("#Delta t (#equiv t_{H}-t_{L}) vs Channel - pos.1 ");
    mg_DELTA_1->Add(DELTA_H_B_1);
    //mg_DELTA_1->Add(DELTA_H_A_1);
    //mg_DELTA_1->Add(DELTA_T_B_1);
    //mg_DELTA_1->Add(DELTA_T_A_1);
    
    
    TMultiGraph *mg_FRAC_0 = new TMultiGraph();
    mg_FRAC_0->SetTitle("Fraction of paths  vs Channel - pos.0 ");
    //mg_FRAC_0->Add(FRAC_H_B_0);
    //mg_FRAC_0->Add(FRAC_H_A_0);
    mg_FRAC_0->Add(FRAC_L_B_0);
    mg_FRAC_0->Add(FRAC_L_A_0);
    
    TMultiGraph *mg_FRAC_1 = new TMultiGraph();
    mg_FRAC_1->SetTitle("Fraction of paths  vs Channel - pos.1 ");
    //mg_FRAC_1->Add(FRAC_H_B_1);
    //mg_FRAC_1->Add(FRAC_H_A_1);
    mg_FRAC_1->Add(FRAC_L_B_1);
    mg_FRAC_1->Add(FRAC_L_A_1);
    
    TMultiGraph *mg_SIGMA_L_0 = new TMultiGraph();
    mg_SIGMA_L_0->SetTitle("Time resolution (#delta t) low time vs Channel - pos. 0");
    //mg_SIGMA_L_0->Add(REF_sigma_L);
    mg_SIGMA_L_0->Add(SIGMA_L_B_0);
    mg_SIGMA_L_0->Add(SIGMA_L_A_0);
    TMultiGraph *mg_SIGMA_H_0 = new TMultiGraph();
    mg_SIGMA_H_0->SetTitle("Time resolution (#delta t) high time vs Channel - pos. 0");
    //mg_SIGMA_H_0->Add(REF_sigma_H);
    mg_SIGMA_H_0->Add(SIGMA_H_B_0);
    mg_SIGMA_H_0->Add(SIGMA_H_A_0);
    TMultiGraph *mg_SIGMA_T_0 = new TMultiGraph();
    mg_SIGMA_T_0->SetTitle("Time resolution (#delta t) third time vs Channel - pos. 0");
    //mg_SIGMA_T_0->Add(REF_sigma_T);
    mg_SIGMA_T_0->Add(SIGMA_T_B_0);
    mg_SIGMA_T_0->Add(SIGMA_T_A_0);
    
    TMultiGraph *mg_SIGMA_L_1 = new TMultiGraph();
    mg_SIGMA_L_1->SetTitle("Time resolution (#delta t) low time vs Channel - pos. 1");
    mg_SIGMA_L_1->Add(SIGMA_L_B_1);
    mg_SIGMA_L_1->Add(SIGMA_L_A_1);
    TMultiGraph *mg_SIGMA_H_1 = new TMultiGraph();
    mg_SIGMA_H_1->SetTitle("Time resolution (#delta t) high time vs Channel - pos. 1");
    mg_SIGMA_H_1->Add(SIGMA_H_B_1);
    mg_SIGMA_H_1->Add(SIGMA_H_A_1);
    TMultiGraph *mg_SIGMA_T_1 = new TMultiGraph();
    mg_SIGMA_T_1->SetTitle("Time resolution (#delta t) third time vs Channel - pos. 1");
    mg_SIGMA_T_1->Add(SIGMA_T_B_1);
    mg_SIGMA_T_1->Add(SIGMA_T_A_1);
    
    TMultiGraph *mg_FRAC_SIGMA_L_0 = new TMultiGraph();
    mg_FRAC_SIGMA_L_0->SetTitle("#delta t [ns] vs f_{L} - pos.0 - low time peak");
    mg_FRAC_SIGMA_L_0->Add(SIGMA_FRAC_L_0);
    mg_FRAC_SIGMA_L_0->Add(SIGMA_FRAC_L_0_A);
    TMultiGraph *mg_FRAC_SIGMA_H_0 = new TMultiGraph();
    mg_FRAC_SIGMA_H_0->SetTitle("#delta t [ns] vs f_{L} - pos.0 - high time peak");
    mg_FRAC_SIGMA_H_0->Add(SIGMA_FRAC_H_0);
    mg_FRAC_SIGMA_H_0->Add(SIGMA_FRAC_H_0_A);
    TMultiGraph *mg_FRAC_SIGMA_L_1 = new TMultiGraph();
    mg_FRAC_SIGMA_L_1->SetTitle("#delta t [ns] vs f_{L} - pos.1 - low time peak");
    mg_FRAC_SIGMA_L_1->Add(SIGMA_FRAC_L_1);
    mg_FRAC_SIGMA_L_1->Add(SIGMA_FRAC_L_1_A);
    TMultiGraph *mg_FRAC_SIGMA_H_1 = new TMultiGraph();
    mg_FRAC_SIGMA_H_1->SetTitle("#delta t [ns] vs f_{L} - pos.1 - high time peak");
    mg_FRAC_SIGMA_H_1->Add(SIGMA_FRAC_H_1);
    mg_FRAC_SIGMA_H_1->Add(SIGMA_FRAC_H_1_A);


    TCanvas *c_frac_sig = new TCanvas("c_frac_sigma","c_frac_sigma");
    c_frac_sig->Divide(2,2);
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
    c_frac_sig->cd(3);
    mg_FRAC_SIGMA_L_1->Draw("AP");
    gPad->Update();
    mg_FRAC_SIGMA_L_1->GetXaxis()->SetTitle("f_{L}");
    mg_FRAC_SIGMA_L_1->GetYaxis()->SetTitle("#delta t [ns]");
    c_frac_sig->cd(4);
    mg_FRAC_SIGMA_H_1->Draw("AP");
    gPad->Update();
    mg_FRAC_SIGMA_H_1->GetXaxis()->SetTitle("f_{L}");
    mg_FRAC_SIGMA_H_1->GetYaxis()->SetTitle("#delta t [ns]");
    //////
    
    TCanvas *c_DELTA_FRAC = new TCanvas("mean peak separation","mean peak separation",0,0,1124,700);
    c_DELTA_FRAC->Divide(1,2);
    c_DELTA_FRAC->cd(1);
    mg_DELTA_0->Draw("AP");
    gPad->Update();
    mg_DELTA_0->GetXaxis()->SetTitle("Channel");
    mg_DELTA_0->GetYaxis()->SetTitle("#Delta t[ns]");
    gPad->Update();
    mg_DELTA_0->Draw("AP");
    gPad->BuildLegend();
    c_DELTA_FRAC->cd(2);
    mg_DELTA_1->Draw("AP");
    gPad->Update();
    mg_DELTA_1->GetXaxis()->SetTitle("Channel");
    mg_DELTA_1->GetYaxis()->SetTitle("#Delta t[ns]");
    gPad->Update();
    mg_DELTA_1->Draw("AP");
    gPad->BuildLegend();
    
    
    TCanvas *c_FRAC = new TCanvas("c_FRAC","c_FRAC",0,0,1124,700);
    c_FRAC->Divide(1,2);
    c_FRAC->cd(1);
    mg_FRAC_0->Draw("AP");
    gPad->Update();
    mg_FRAC_0->GetXaxis()->SetTitle("Channel");
    mg_FRAC_0->GetYaxis()->SetTitle("#Delta t[ns]");
    gPad->Update();
    mg_FRAC_0->Draw("AP");
    gPad->BuildLegend();
    c_FRAC->cd(2);
    mg_FRAC_1->Draw("AP");
    gPad->Update();
    mg_FRAC_1->GetXaxis()->SetTitle("Channel");
    mg_FRAC_1->GetYaxis()->SetTitle("#Delta t[ns]");
    gPad->Update();
    mg_FRAC_1->Draw("AP");
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
    
    TCanvas* c_Res_pos1 = new TCanvas("c_Res_pos1","c_Res_pos1",0,0,1124,700);
    c_Res_pos1->Divide(2,1);
    c_Res_pos1->cd(1);
    mg_SIGMA_L_1->Draw("AP");
    line->Draw("same");
    gPad->BuildLegend();
    gPad->Update();
    mg_SIGMA_L_1->GetXaxis()->SetTitle("Channel");
    mg_SIGMA_L_1->GetYaxis()->SetTitle("#delta t [ns]");
    gPad->Update();
    c_Res_pos1->cd(2);
    mg_SIGMA_H_1->Draw("AP");
    line->Draw("same");
    gPad->BuildLegend();
    gPad->Update();
    mg_SIGMA_H_1->GetXaxis()->SetTitle("Channel");
    mg_SIGMA_H_1->GetYaxis()->SetTitle("#delta t [ns]");
    gPad->Update();
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
      
      TGraph *GR_Fit_status_0 = new TGraph(16,X,Fit_status_0);GR_Fit_status_0->SetTitle("Fits status - fit pos 0");
      TGraph *GR_Fit_status_1 = new TGraph(16,X,Fit_status_1);GR_Fit_status_1->SetTitle("Fits status - fit pos 1");
      TGraph *GR_Fit_status = new TGraph(16,X,Fit_status);GR_Fit_status->SetTitle("Fits status - fit pos 0 #oplus 1");
      GR_Fit_status_0->SetMarkerStyle(20);
      GR_Fit_status_0->SetMarkerColor(1);
      GR_Fit_status_1->SetMarkerStyle(21);
      GR_Fit_status_1->SetMarkerColor(2);
      GR_Fit_status->SetMarkerStyle(22);
      GR_Fit_status->SetMarkerColor(3);
      TMultiGraph *mg_FITSTATUS = new TMultiGraph();
      mg_FITSTATUS->SetTitle("Fits status");
      mg_FITSTATUS->Add(GR_Fit_status_0);
      mg_FITSTATUS->Add(GR_Fit_status_1);
      mg_FITSTATUS->Add(GR_Fit_status);
    */
    
    c_pos0_AllChannels->Draw();    c_pos0_AllChannels_pulls->Draw();  c_pos0_AllChannels_corr->Draw(); 
    c_pos1_AllChannels->Draw();    c_pos1_AllChannels_pulls->Draw();  c_pos1_AllChannels_corr->Draw();
    c_pos01_AllChannels->Draw();   c_pos01_AllChannels_pulls->Draw();   c_pos01_AllChannels_corr->Draw();    
    
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













