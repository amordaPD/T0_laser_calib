////////TIME CALIBRATION CODE///////****
////////////////////////////////////
/////To RUN: gROOT->ProcessLine(".L Fit_single_fiber.cpp"); loop_channels(2,true); > ee.log

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
  RooPlot* xframe2_amplitude_0;
  RooPlot* xframe2_fit_0;
  RooPlot* xframe2_pull_0;
  TH2 *h_correlation_0;
  
};


//  ofstream cout("fits_report.txt");

//TFile *f_input_histogram = new TFile("flat_ntuples/run100317-1-T77-vecchialente_out.root");
//TFile *f_input_histogram = new TFile("flat_ntuples/run100317-2-T77-nuovalente-p0_out.root");
TFile *f_input_histogram = new TFile("flat_ntuples/run100317-3-T77-nuovalente-p1_out.root");
//TFile *f_input_histogram = new TFile("flat_ntuples/grease-1-T77_out.root");

//vector<float> Fit_head(string _draw_results, int fix_params, int ch ){
Fit_results Fit_head(string _draw_results="draw", int fix_params=2, int ch =0 ){
  
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
  bool no_grease=false;

  
  if(!print_prefit_info){MN_output_print_level_prefit=-1;}else{MN_output_print_level_prefit=MN_output_print_level;}
  bool draw_results;
  if(_draw_results=="draw"){draw_results=true;}else if(_draw_results=="blind"){draw_results=false;}else{draw_results=false;}
  Fit_results my_fit_results;
  TString channel_0 = Form("fiber0-0-%d",ch);
  //bool fix_params = true;
  vector<float> POIs;
  POIs.clear();
  // S e t u p   m o d e l 
  // ---------------------
 
  double my_low_x=21.5;
  double my_up_x=24.5;
  RooRealVar x("Time","Time [ns]",my_low_x,my_up_x) ;
  RooRealVar amp("Amplitude","Amplitude [ADC counts]",-200,0) ;
  RooRealVar CH("Channel","PMT Channel",0,16) ;
  ///////fixing starting values and boundaries for two positions


  
  TTree *tree = (TTree*)f_input_histogram->Get("tree_input");
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
   up_sigma_H_0=0.180;
  
  low_sigma_T_0=0.050;
   up_sigma_T_0=0.500;
  
  low_alpha_0=0.01;
   up_alpha_0=0.8;
  
  low_beta_0=0.5;
   up_beta_0=1.0;

   if(no_grease){
     starting_mean_L_0=22.7;
     if(ch==12||ch==14) starting_mean_L_0=22.5;
   }else{
     starting_mean_L_0=22.7;
     if(ch==14) starting_mean_L_0=22.5;

   }
   //if(ch==2) starting_mean_L_0=23.0;
   //   if(ch==3||ch==7) starting_mean_L_0=8.4;
   starting_delta_H_0=0.280;
   starting_delta_T_0=0.500;
   starting_sigma_L_0=0.0810;
   starting_sigma_H_0=0.0810;
   starting_sigma_T_0=0.1000;
   starting_alpha_0=0.5;
   starting_beta_0=0.5;
   
   //  low_mean_L_0=starting_mean_L_0-0.15;
   //   up_mean_L_0=starting_mean_L_0+0.15;
  low_mean_L_0=starting_mean_L_0-0.35;
   up_mean_L_0=starting_mean_L_0+0.35;
  starting_mean_H_0=starting_mean_L_0+starting_delta_H_0;
   low_mean_H_0=low_mean_L_0+starting_delta_H_0;
   up_mean_H_0=up_mean_L_0+starting_delta_H_0;





  starting_alpha_CB=-0.5;
  starting_n_CB=4;

  low_alpha_CB=-5;
   up_alpha_CB=0.0;
   
  low_n_CB=0;
   up_n_CB=10;

  
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
  RooRealVar Delta_H_0("Delta_H_0","Delta H pos 0",starting_delta_H_0,low_delta_H_0,up_delta_H_0);
  RooRealVar Delta_T_0("Delta_T_0","Delta T pos 0",starting_delta_T_0,low_delta_T_0,up_delta_T_0);
  RooFormulaVar mean_L_0("mean_L_0","mean_l_0","mean_H_0-Delta_H_0",RooArgList(mean_H_0,Delta_H_0));
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


  
  RooRealVar a0_0("a0_0", "", 0.0, -10, 10);
  RooRealVar a1_0("a1_0", "", 0.0, -20, 20);
  RooRealVar a2_0("a2_0", "", 0.0015, -20, 20);
  RooChebychev PDF_B_0("PDF_B_0","PDF_B_0",x,RooArgList(a0_0,a1_0));//,a2_0));//

  RooRealVar  Frac_sig_0("Frac_sig_0","fraction of sig events", 0.9, 0.7,1.0);
  RooArgList  pdfList_0(PDF_sig_0,PDF_B_0);//
  RooArgList  fracList_0(Frac_sig_0);
  RooAddPdf   model_0("model_0","model_0",pdfList_0,fracList_0,kTRUE);

  


  
  RooAddPdf   model_0_b("model_0_b","model_0_b",pdfList_0,fracList_0,kTRUE);







  
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
  
  
  
  RooFitResult* fit_results_0;
  if(use_NLL){
    RooAbsReal* nll_model_0 = model_0.createNLL(ds_0,Extended(kFALSE)) ;
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
  my_fit_results.n_POIs=POIs.size();
  my_fit_results.xframe2_amplitude_0=xframe2_0_amp;
  my_fit_results.xframe2_fit_0=xframe2_0;
  my_fit_results.xframe2_pull_0=xframe4_0;
  my_fit_results.h_correlation_0=h_correlation_0;
  
  

  //  delete  h_input_histogram_0;
  
  
  return my_fit_results;
  
  //return POIs;
}



vector<float> loop_channels(int deep_fixed_params,bool plot_summaries){ //rel_weight is the relative weight of the 1 to 0 pois, i.e. total dataset= pos.0 + rel_weight*pos.1 - in the real case rel_weight=1.0
  
  
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
    //mg_MEANS_L_0->Add(MEAN_L_A_0);
    
    TCanvas *c_absolute_positions_01_b = new TCanvas("c_absolute_position_01_b","c_absolute_position_01_b",0,0,1124,700);
    //c_absolute_positions_01->Divide(2,1);
    //c_absolute_positions_01->cd(1);
    mg_MEANS_L_0->Draw("AP");
    mg_MEANS_L_0->GetXaxis()->SetTitle("Channel");
    mg_MEANS_L_0->GetYaxis()->SetTitle("T_{L}^{pos.0} [ns]");
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
    TCanvas *c_pos0_deltas_data_MC = new TCanvas("c_pos0_deltas_data_MC","c_pos0_deltas_data_MC",0,0,1124,700);
    mg_DELTA_B_pos0->Draw("AP");
    gPad->BuildLegend();
    gPad->Update();
    gPad->SetGridy();
    mg_DELTA_B_pos0->GetXaxis()->SetTitle("Channel");
    mg_DELTA_B_pos0->GetYaxis()->SetTitle("#Delta t [ns]");
    mg_DELTA_B_pos0->GetYaxis()->SetNdivisions(20);
    gPad->Update();
    
    
    
    
    TMultiGraph *mg_SIGMA_B_L = new TMultiGraph();
    mg_SIGMA_B_L->SetTitle("Time resolution (#delta t) low time vs Channel");
    mg_SIGMA_B_L->Add(REF_sigma_L);
    mg_SIGMA_B_L->Add(SIGMA_L_B_0);
    //mg_SIGMA_B_L->Add(SIGMA_L_A_0);
    TMultiGraph *mg_SIGMA_B_H = new TMultiGraph();
    mg_SIGMA_B_H->SetTitle("Time resolution (#delta t) high time vs Channel");
    mg_SIGMA_B_H->Add(REF_sigma_H);
    mg_SIGMA_B_H->Add(SIGMA_H_B_0);
    //mg_SIGMA_B_H->Add(SIGMA_H_A_0);
    TMultiGraph *mg_SIGMA_B_T = new TMultiGraph();
    mg_SIGMA_B_T->SetTitle("Time resolution (#delta t) third time vs Channel");
    mg_SIGMA_B_T->Add(REF_sigma_T);
    //mg_SIGMA_B_T->Add(SIGMA_T_B_0);
    //mg_SIGMA_B_T->Add(SIGMA_T_A_0);
    mg_SIGMA_B_T->Add(SIGMA_T_B_0);
    //mg_SIGMA_B_T->Add(SIGMA_T_A_1);
    
    
    
    
    ///////NEW MG_GRAPHS/////
    TMultiGraph *mg_DELTA_0 = new TMultiGraph();
    mg_DELTA_0->SetTitle("#Delta t (#equiv t_{H}-t_{L}) vs Channel - pos.0 ");
    //mg_DELTA_0->Add(DELTA_means_ref);
    mg_DELTA_0->Add(DELTA_H_B_0);
    //mg_DELTA_0->Add(DELTA_H_A_0);
    //mg_DELTA_0->Add(DELTA_T_B_0);
    //mg_DELTA_0->Add(DELTA_T_A_0);

    
    
    TMultiGraph *mg_FRAC_0 = new TMultiGraph();
    mg_FRAC_0->SetTitle("Fraction of paths  vs Channel - pos.0 ");
    //mg_FRAC_0->Add(FRAC_H_B_0);
    //mg_FRAC_0->Add(FRAC_H_A_0);
    mg_FRAC_0->Add(FRAC_L_B_0);
    //mg_FRAC_0->Add(FRAC_L_A_0);

    
    TMultiGraph *mg_SIGMA_L_0 = new TMultiGraph();
    mg_SIGMA_L_0->SetTitle("Time resolution (#delta t) low time vs Channel - pos. 0");
    //mg_SIGMA_L_0->Add(REF_sigma_L);
    mg_SIGMA_L_0->Add(SIGMA_L_B_0);
    //mg_SIGMA_L_0->Add(SIGMA_L_A_0);
    TMultiGraph *mg_SIGMA_H_0 = new TMultiGraph();
    mg_SIGMA_H_0->SetTitle("Time resolution (#delta t) high time vs Channel - pos. 0");
    //mg_SIGMA_H_0->Add(REF_sigma_H);
    mg_SIGMA_H_0->Add(SIGMA_H_B_0);
    //mg_SIGMA_H_0->Add(SIGMA_H_A_0);
    TMultiGraph *mg_SIGMA_T_0 = new TMultiGraph();
    mg_SIGMA_T_0->SetTitle("Time resolution (#delta t) third time vs Channel - pos. 0");
    //mg_SIGMA_T_0->Add(REF_sigma_T);
    mg_SIGMA_T_0->Add(SIGMA_T_B_0);
    //mg_SIGMA_T_0->Add(SIGMA_T_A_0);
    
    
    TMultiGraph *mg_FRAC_SIGMA_L_0 = new TMultiGraph();
    mg_FRAC_SIGMA_L_0->SetTitle("#delta t [ns] vs f_{L} - pos.0 - low time peak");
    mg_FRAC_SIGMA_L_0->Add(SIGMA_FRAC_L_0);
    mg_FRAC_SIGMA_L_0->Add(SIGMA_FRAC_L_0_A);
    TMultiGraph *mg_FRAC_SIGMA_H_0 = new TMultiGraph();
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













