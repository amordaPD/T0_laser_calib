#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
using namespace RooFit ;

struct Fit_results{
  float mean[2];
  float delta[2];
  float Av_N_ph[2];
  float CT_prob[2];
  float s0[2];
  float s1[2];
};


Fit_results Fit_Laser_Stability(int fit_model=3, float T_central, float T_range){
  // S e t u p   m o d e l 
  // ---------------------
  Fit_results Results;
  //RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  TString string_out="";
  if(fit_model==1){string_out="Poisson_NoOCT";}
  else if (fit_model==2){string_out="Poisson_OCT";}
  else if (fit_model==3){string_out="Free_fractions";}
  
  gROOT->ProcessLine(".x ./myRooPdfs/RooExpGauss.cxx+") ;
  gROOT->ProcessLine(".x ./myRooPdfs/RooAsymGauss.cxx+") ;
  
  int x_max=500;//1000;
  RooRealVar x("amps","amps",0,x_max) ;//
  RooRealVar mean("mean","zero_value",25,20,30) ; //position of first peak
  RooRealVar Delta("Delta","Delta",40,30,120); //separation between two consecutive peaks - proportional to the gain of SiPM
  RooRealVar Av_N_ph("Av_N_ph","Average number of photons",0.2,0.,20); //Average number of photons
  RooRealVar q("q","q_{CT}", 0.1,0.0,1.0); //CT probability
  

  
  RooRealVar N_Obs_0("N_Obs_0","Observed number of photons",0);
  RooRealVar N_Obs_1("N_Obs_1","Observed number of photons",1);
  RooRealVar N_Obs_2("N_Obs_2","Observed number of photons",2);
  RooRealVar N_Obs_3("N_Obs_3","Observed number of photons",3);
  RooRealVar N_Obs_4("N_Obs_4","Observed number of photons",4);
  RooRealVar N_Obs_5("N_Obs_5","Observed number of photons",5);
  RooRealVar N_Obs_6("N_Obs_6","Observed number of photons",6);
  RooRealVar N_Obs_7("N_Obs_7","Observed number of photons",7);
  RooRealVar N_Obs_8("N_Obs_8","Observed number of photons",8);
  RooRealVar N_Obs_9("N_Obs_9","Observed number of photons",9);
  RooRealVar N_Obs_10("N_Obs_10","Observed number of photons",10);
  RooRealVar N_Obs_11("N_Obs_11","Observed number of photons",11);
  RooRealVar N_Obs_12("N_Obs_12","Observed number of photons",12);
  RooRealVar N_Obs_13("N_Obs_13","Observed number of photons",13);
  RooRealVar N_Obs_14("N_Obs_14","Observed number of photons",14);
  RooRealVar N_Obs_15("N_Obs_15","Observed number of photons",15);
  
  /*
    RooRealVar mean_0("mean_0","#mu_{0}",0,-5,5) ;
    RooRealVar mean_1("mean_1","#mu_{1}",100,95,130) ;
    RooRealVar mean_2("mean_2","#mu_{2}",200,195,230) ;
    RooRealVar mean_3("mean_3","#mu_{3}",300,295,320) ;
    RooRealVar mean_4("mean_4","#mu_{4}",400,395,450) ;
    RooRealVar mean_5("mean_5","#mu_{5}",500,495,520) ;
    RooRealVar mean_6("mean_6","#mu_{6}",600,590,630) ;
    RooRealVar mean_7("mean_7","#mu_{7}",700,650,730) ;
    RooRealVar mean_8("mean_8","#mu_{8}",800,795,850) ;
    RooRealVar mean_9("mean_9","#mu_{9}",900,895,920) ;
    RooRealVar mean_10("mean_10","#mu_{10}",1000,990,1030) ;
    RooRealVar mean_11("mean_11","#mu_{11}",1100,1090,1130) ;
  */
  /*
  RooRealVar sigma_0("sigma_0","#sigma_{0}",10,1,20) ;
  RooRealVar sigma_1("sigma_1","#sigma_{1}",10,1,20) ;
  RooRealVar sigma_2("sigma_2","#sigma_{2}",10,1,20) ;
  RooRealVar sigma_3("sigma_3","#sigma_{3}",10,1,20) ;
  RooRealVar sigma_4("sigma_4","#sigma_{4}",10,2,20) ;
  RooRealVar sigma_5("sigma_5","#sigma_{5}",10,2,20) ;
  RooRealVar sigma_6("sigma_6","#sigma_{6}",10,2,20) ;
  RooRealVar sigma_7("sigma_7","#sigma_{7}",10,2,20) ;
  RooRealVar sigma_8("sigma_8","#sigma_{8}",10,2,20) ;
  RooRealVar sigma_9("sigma_9","#sigma_{9}",10,2,20) ;
  RooRealVar sigma_10("sigma_10","#sigma_{10}",10,2,20) ;
  RooRealVar sigma_11("sigma_11","#sigma_{11}",10,2,20) ;
  RooRealVar sigma_12("sigma_12","#sigma_{12}",10,2,20) ;
  RooRealVar sigma_13("sigma_13","#sigma_{13}",10,2,20) ;
  RooRealVar sigma_14("sigma_14","#sigma_{14}",10,2,20) ;
  RooRealVar sigma_15("sigma_15","#sigma_{15}",10,2,20) ;
  */
  
  RooRealVar s0("s0","s0",5.,4.,6.0);//
  RooRealVar s1("s1","s1",0.1,0.05,1.5);
  RooRealVar sigma_0("sigma_0","#sigma_{0}",10,1,20) ;
  RooFormulaVar sigma_1("sigma_1","sigma_1","s0+0*s1",RooArgList(s0,s1));
  RooFormulaVar sigma_2("sigma_2","sigma_2","s0+1*s1",RooArgList(s0,s1));
  RooFormulaVar sigma_3("sigma_3","sigma_3","s0+2*s1",RooArgList(s0,s1));
  RooFormulaVar sigma_4("sigma_4","sigma_4","s0+3*s1",RooArgList(s0,s1));
  RooFormulaVar sigma_5("sigma_5","sigma_5","s0+4*s1",RooArgList(s0,s1));
  RooFormulaVar sigma_6("sigma_6","sigma_6","s0+5*s1",RooArgList(s0,s1));
  RooFormulaVar sigma_7("sigma_7","sigma_7","s0+6*s1",RooArgList(s0,s1));
  RooFormulaVar sigma_8("sigma_8","sigma_8","s0+7*s1",RooArgList(s0,s1));
  RooFormulaVar sigma_9("sigma_9","sigma_9","s0+8*s1",RooArgList(s0,s1));
  RooFormulaVar sigma_10("sigma_10","sigma_10","s0+9*s1",RooArgList(s0,s1));
  RooFormulaVar sigma_11("sigma_11","sigma_11","s0+10*s1",RooArgList(s0,s1));
  RooFormulaVar sigma_12("sigma_12","sigma_12","s0+11*s1",RooArgList(s0,s1));
  RooFormulaVar sigma_13("sigma_13","sigma_13","s0+12*s1",RooArgList(s0,s1));
  RooFormulaVar sigma_14("sigma_14","sigma_14","s0+13*s1",RooArgList(s0,s1));
  RooFormulaVar sigma_15("sigma_15","sigma_15","s0+14*s1",RooArgList(s0,s1));
  

  //RooRealVar    mean_0("mean_0","mean of gaussian background",0,-5,5) ;
  RooFormulaVar mean_1("mean_1","mean_1","mean+0*Delta",RooArgList(mean,Delta)) ;
  RooFormulaVar mean_2("mean_2","mean_2","mean+1*Delta",RooArgList(mean,Delta)) ;
  RooFormulaVar mean_3("mean_3","mean_3","mean+2*Delta",RooArgList(mean,Delta)) ;
  RooFormulaVar mean_4("mean_4","mean_4","mean+3*Delta",RooArgList(mean,Delta)) ;
  RooFormulaVar mean_5("mean_5","mean_5","mean+4*Delta",RooArgList(mean,Delta)) ;
  RooFormulaVar mean_6("mean_6","mean_6","mean+5*Delta",RooArgList(mean,Delta)) ;
  RooFormulaVar mean_7("mean_7","mean_7","mean+6*Delta",RooArgList(mean,Delta)) ;
  RooFormulaVar mean_8("mean_8","mean_8","mean+7*Delta",RooArgList(mean,Delta)) ;
  RooFormulaVar mean_9("mean_9","mean_9","mean+8*Delta",RooArgList(mean,Delta)) ;
  RooFormulaVar mean_10("mean_10","mean_10","mean+9*Delta",RooArgList(mean,Delta)) ;
  RooFormulaVar mean_11("mean_11","mean_11","mean+10*Delta",RooArgList(mean,Delta)) ;
  RooFormulaVar mean_12("mean_12","mean_12","mean+11*Delta",RooArgList(mean,Delta)) ;
  RooFormulaVar mean_13("mean_13","mean_13","mean+12*Delta",RooArgList(mean,Delta)) ;
  RooFormulaVar mean_14("mean_14","mean_14","mean+13*Delta",RooArgList(mean,Delta)) ;
  RooFormulaVar mean_15("mean_15","mean_15","mean+14*Delta",RooArgList(mean,Delta)) ;
  
  
  
  
  RooRealVar TP("TP","Transition Point of gaussExp",1,0.5,5);
  RooRealVar sigma_bkg("sigma_bkg","sigma_bkg",500,50,1000);
  RooRealVar mean_bkg("mean_bkg","mean_bkg",300,200,500);
  RooGaussian PDF_bkg("PDF_bkg","PDF background",x,mean_bkg,sigma_bkg) ;
  RooRealVar  Frac_sig("Frac_sig"  ,"fraction of signal events", 0.9 ,0.8,1);

  
  RooRealVar alpha_15("alpha_15","Observed number of photons",-0.5,-10,0);
  RooGaussian PDF_0("PDF_0","gaussian PDF background",x,mean,sigma_0) ;
  RooGaussian PDF_1("PDF_1","gaussian PDF background",x,mean_1,sigma_1) ;
  RooGaussian PDF_2("PDF_2","gaussian PDF background",x,mean_2,sigma_2) ;
  RooGaussian PDF_3("PDF_3","gaussian PDF background",x,mean_3,sigma_3) ;
  RooGaussian PDF_4("PDF_4","gaussian PDF background",x,mean_4,sigma_4) ;
  RooGaussian PDF_5("PDF_5","gaussian PDF background",x,mean_5,sigma_5) ;
  RooGaussian PDF_6("PDF_6","gaussian PDF background",x,mean_6,sigma_6) ;
  RooGaussian PDF_7("PDF_7","gaussian PDF background",x,mean_7,sigma_7) ;
  RooGaussian PDF_8("PDF_8","gaussian PDF background",x,mean_8,sigma_8) ;
  RooGaussian PDF_9("PDF_9","gaussian PDF background",x,mean_9,sigma_9) ;
  RooGaussian PDF_10("PDF_10","gaussian PDF background",x,mean_10,sigma_10) ;
  RooExpGauss PDF_11("PDF_11","gaussian PDF background",x,mean_11,sigma_11,alpha_15) ;
  RooGaussian PDF_12("PDF_12","gaussian PDF background",x,mean_12,sigma_12) ;
  RooGaussian PDF_13("PDF_13","gaussian PDF background",x,mean_13,sigma_13) ;
  RooGaussian PDF_14("PDF_14","gaussian PDF background",x,mean_14,sigma_14) ;
  RooExpGauss PDF_15("PDF_15","gaussian PDF background",x,mean_15,sigma_15,alpha_15) ;
  

  /*
  RooExpGauss PDF_0("PDF_0","gaussian PDF background",x,mean,sigma_0,TP) ;
  RooExpGauss PDF_1("PDF_1","gaussian PDF background",x,mean_1,sigma_1,TP) ;
  RooExpGauss PDF_2("PDF_2","gaussian PDF background",x,mean_2,sigma_2,TP) ;
  RooExpGauss PDF_3("PDF_3","gaussian PDF background",x,mean_3,sigma_3,TP) ;
  RooExpGauss PDF_4("PDF_4","gaussian PDF background",x,mean_4,sigma_4,TP) ;
  RooExpGauss PDF_5("PDF_5","gaussian PDF background",x,mean_5,sigma_5,TP) ;
  RooExpGauss PDF_6("PDF_6","gaussian PDF background",x,mean_6,sigma_6,TP) ;
  RooExpGauss PDF_7("PDF_7","gaussian PDF background",x,mean_7,sigma_7,TP) ;
  RooExpGauss PDF_8("PDF_8","gaussian PDF background",x,mean_8,sigma_8,TP) ;
  RooExpGauss PDF_9("PDF_9","gaussian PDF background",x,mean_9,sigma_9,TP) ;
  RooExpGauss PDF_10("PDF_10","gaussian PDF background",x,mean_10,sigma_10,TP) ;
  RooExpGauss PDF_11("PDF_11","gaussian PDF background",x,mean_11,sigma_11,TP) ;
  RooExpGauss PDF_12("PDF_12","gaussian PDF background",x,mean_12,sigma_12,TP) ;
  RooExpGauss PDF_13("PDF_13","gaussian PDF background",x,mean_13,sigma_13,TP) ;
  RooExpGauss PDF_14("PDF_14","gaussian PDF background",x,mean_14,sigma_14,TP) ;
  RooRealVar alpha_15("alpha_15","Observed number of photons",-0.5,-10,0);
  RooExpGauss PDF_15("PDF_15","gaussian PDF background",x,mean_15,sigma_15,alpha_15) ;
  */



  
  if(fit_model==1){
    
    ////////////////////////////////////////////////////////////////////
    ///////////////////////////// POISSON  /////////////////////////////
    ////////////////////////////////////////////////////////////////////
    RooPoisson  Frac_0( "Frac_0","fraction of 0 ph signal events"     ,N_Obs_0,Av_N_ph);
    RooPoisson  Frac_1( "Frac_1","fraction of 1 ph signal events"     ,N_Obs_1,Av_N_ph);
    RooPoisson  Frac_2( "Frac_2","fraction of 2 ph signal events"     ,N_Obs_2,Av_N_ph);
    RooPoisson  Frac_3( "Frac_3","fraction of 3 ph signal events"     ,N_Obs_3,Av_N_ph);
    RooPoisson  Frac_4( "Frac_4","fraction of 4 ph signal events"     ,N_Obs_4,Av_N_ph);
    RooPoisson  Frac_5( "Frac_5","fraction of 5 ph signal events"     ,N_Obs_5,Av_N_ph);
    RooPoisson  Frac_6( "Frac_6","fraction of 6 ph signal events"     ,N_Obs_6,Av_N_ph);
    RooPoisson  Frac_7( "Frac_7","fraction of 7 ph signal events"     ,N_Obs_7,Av_N_ph);
    RooPoisson  Frac_8( "Frac_8","fraction of 8 ph signal events"     ,N_Obs_8,Av_N_ph);
    RooPoisson  Frac_9( "Frac_9","fraction of 9 ph signal events"     ,N_Obs_9,Av_N_ph);
    RooPoisson  Frac_10( "Frac_10","fraction of 10 ph signal events"  ,N_Obs_10,Av_N_ph);
    RooRealVar  Frac_11( "Frac_11"  ,"Frac_{11}", 0.006,0.0001,0.5);
    //RooPoisson  Frac_11( "Frac_11","fraction of 11 ph signal events"  ,N_Obs_11,Av_N_ph);
    RooPoisson  Frac_12( "Frac_12","fraction of 12 ph signal events"  ,N_Obs_12,Av_N_ph);
    RooPoisson  Frac_13( "Frac_13","fraction of 13 ph signal events"  ,N_Obs_13,Av_N_ph);
    RooPoisson  Frac_14( "Frac_14","fraction of 13 ph signal events"  ,N_Obs_14,Av_N_ph);
    RooRealVar  Frac_15( "Frac_15"  ,"Frac_{15}", 0.006,0.0001,0.5);
 
    
  }else if(fit_model==2){
    
    
    ////////////////////////////////////////////////////////////////////
    //////////////////POISSON + OPTICAL CROSS TALKS ////////////////////
    ////////////////////////////////////////////////////////////////////
    RooPoisson  frac_0( "frac_0","fraction of 0 ph signal events"  ,N_Obs_0,Av_N_ph);
    RooPoisson  frac_1( "frac_1","fraction of 1 ph signal events"  ,N_Obs_1,Av_N_ph);
    RooPoisson  frac_2( "frac_2","fraction of 2 ph signal events"  ,N_Obs_2,Av_N_ph);
    RooPoisson  frac_3( "frac_3","fraction of 3 ph signal events"  ,N_Obs_3,Av_N_ph);
    RooPoisson  frac_4( "frac_4","fraction of 4 ph signal events"  ,N_Obs_4,Av_N_ph);
    RooPoisson  frac_5( "frac_5","fraction of 5 ph signal events"  ,N_Obs_5,Av_N_ph);
    RooPoisson  frac_6( "frac_6","fraction of 6 ph signal events"  ,N_Obs_6,Av_N_ph);
    RooPoisson  frac_7( "frac_7","fraction of 7 ph signal events"  ,N_Obs_7,Av_N_ph);
    RooPoisson  frac_8( "frac_8","fraction of 8 ph signal events"  ,N_Obs_8,Av_N_ph);
    RooPoisson  frac_9( "frac_9","fraction of 9 ph signal events"  ,N_Obs_9,Av_N_ph);
    RooPoisson  frac_10( "frac_10","fraction of 10 ph signal events"  ,N_Obs_10,Av_N_ph);
    
    RooFormulaVar  Frac_0( "Frac_0","F_{0}","frac_0*1",RooArgList(frac_0));
    RooFormulaVar  Frac_1( "Frac_1","F_{1}","frac_1*(1-q)" ,RooArgList(frac_1,q));
    RooFormulaVar  Frac_2( "Frac_2","F_{2}","frac_2*(1-q)*(1-q)+frac_1*q",RooArgList(frac_2,q,frac_1));
    RooFormulaVar  Frac_3( "Frac_3","F_{3}","frac_3*(1-q)*(1-q)*(1-q)+frac_2*2*q*(1-q)",RooArgList(frac_3,q,frac_2));
    RooFormulaVar  Frac_4( "Frac_4","F_{4}","frac_4*TMath::Power(1-q,4)+frac_3*3*q*TMath::Power(1-q,2)+frac_2*2*TMath::Power(q,2)",RooArgList(frac_4,q,frac_3,frac_2));
    RooFormulaVar  Frac_5( "Frac_5","F_{5}","frac_5*TMath::Power(1-q,5)+frac_4*4*q*TMath::Power(1-q,3)+frac_3*3*2*q*q*(1-q)",RooArgList(frac_5,q,frac_4,frac_3));
    RooFormulaVar  Frac_6( "Frac_6","F_{6}","frac_6*TMath::Power(1-q,6)+frac_5*5*q*TMath::Power(1-q,4)+frac_4*12*q*q*(1-q)*(1-q)+frac_3*3*2*q*q*q",RooArgList(frac_6,q,frac_5,frac_4,frac_3));
    RooFormulaVar  Frac_7( "Frac_7","F_{7}","frac_7*TMath::Power(1-q,7)+frac_6*6*q*(1-q)*(1-q)*(1-q)*(1-q)*(1-q)+frac_5*20*q*q*(1-q)*(1-q)*(1-q)+frac_4*4*3*2*q*q*q*(1-q)",RooArgList(frac_7,q,frac_6,frac_5,frac_4));
    RooFormulaVar  Frac_8( "Frac_8","F_{9}","frac_8*TMath::Power(1-q,8)+frac_7*7*q*TMath::Power(1-q,6)+frac_6*30*q*q*TMath::Power(1-q,4)+frac_5*5*4*3*q*q*q*(1-q)*(1-q)+frac_4*4*3*2*q*q*q*q",RooArgList(frac_8,q,frac_7,frac_6,frac_5,frac_4));
    RooFormulaVar  Frac_9( "Frac_9","F_{9}","frac_9*TMath::Power(1-q,9)+frac_8*8*q*TMath::Power(1-q,7)+frac_7*42*q*q*TMath::Power(1-q,5)+frac_6*6*5*4*q*q*q*TMath::Power(1-q,3)+frac_5*5*4*3*2*q*q*q*q*(1-q)*(1-q)*(1-q)",RooArgList(frac_9,q,frac_8,frac_7,frac_6,frac_5));
    RooFormulaVar  Frac_10( "Frac_10","F_{10}","frac_10*TMath::Power(1-q,10)+frac_9*9*q*TMath::Power(1-q,8)+frac_8*56*TMath::Power(q,2)*TMath::Power(1-q,6)+frac_7*7*6*5*TMath::Power(q,3)*TMath::Power(1-q,4)+frac_6*6*5*4*3*TMath::Power(q,4)*TMath::Power(1-q,2)+frac_5*5*4*3*2*q*q*q*q*q",RooArgList(frac_10,q,frac_9,frac_8,frac_7,frac_6,frac_5));
    
    RooRealVar  Frac_11("Frac_11"  ,"Frac_{11}", 0.006,0.0001,0.5);
    RooRealVar  Frac_12("Frac_12"  ,"Frac_{12}", 0.09,0.001,0.5);
    RooRealVar  Frac_13("Frac_13"  ,"Frac_{13}", 0.006,0.0001,0.5);
    RooRealVar  Frac_14("Frac_14"  ,"Frac_{14}", 0.006,0.0001,0.5);
    RooRealVar  Frac_15("Frac_15"  ,"Frac_{15}", 0.006,0.0001,0.5);
    
    
  }else if(fit_model==3){
    
    
    ////////////////////////////////////////////////////////////////////
    ///////////////////////////// FREE FRACTIONS////////////////////////
    ////////////////////////////////////////////////////////////////////
    RooRealVar  Frac_0("Frac_0"    ,"Frac_{0}", 0.25 ,0,1);
    RooRealVar  Frac_1("Frac_1"    ,"Frac_{1}", 0.20,0,0.5);
    RooRealVar  Frac_2("Frac_2"    ,"Frac_{2}", 0.10,0,0.5);
    RooRealVar  Frac_3("Frac_3"    ,"Frac_{3}", 0.07,0,0.5);
    RooRealVar  Frac_4("Frac_4"    ,"Frac_{4}", 0.04,0,0.5);
    RooRealVar  Frac_5("Frac_5"    ,"Frac_{5}", 0.025,0,0.5);
    RooRealVar  Frac_6("Frac_6"    ,"Frac_{6}", 0.09,0.001,0.5);
    RooRealVar  Frac_7("Frac_7"    ,"Frac_{7}", 0.006,0.001,0.5);
    RooRealVar  Frac_8("Frac_8"    ,"Frac_{8}", 0.04,0,0.5);
    RooRealVar  Frac_9("Frac_9"    ,"Frac_{9}", 0.025,0,0.5);
    RooRealVar  Frac_10("Frac_10"  ,"Frac_{10}", 0.09,0.001,0.5);
    RooRealVar  Frac_11("Frac_11"  ,"Frac_{11}", 0.006,0.0001,0.5);
    RooRealVar  Frac_12("Frac_12"  ,"Frac_{12}", 0.09,0.001,0.5);
    RooRealVar  Frac_13("Frac_13"  ,"Frac_{13}", 0.006,0.0001,0.5);
    RooRealVar  Frac_14("Frac_14"  ,"Frac_{14}", 0.006,0.0001,0.5);
    RooRealVar  Frac_15("Frac_15"  ,"Frac_{15}", 0.006,0.0001,0.5);

  }				    
 
     
  RooArgList pdfList(PDF_1,PDF_2,PDF_3,PDF_4,PDF_5,PDF_6,PDF_7,PDF_8);
  pdfList.add(RooArgList(PDF_9,PDF_10,PDF_11));//,PDF_12,PDF_13,PDF_14,PDF_15));
  RooArgList fracList(Frac_1,Frac_2,Frac_3,Frac_4,Frac_5,Frac_6,Frac_7,Frac_8);
  fracList.add(RooArgList(Frac_9,Frac_10,Frac_11));//,Frac_12,Frac_13,Frac_14,Frac_15));
	       
  RooAddPdf model("model","model",pdfList,fracList);
  
  // I m p o r t    d a t a
  // ------------------------
  TString input_basepath ="dati/";
  TFile *f_input = new TFile(input_basepath+"long_run_T50_out_temp.root");
  TTree *t_input = (TTree*)f_input->Get("times");
  TH1D *h_amps = new TH1D("h_amps","h_amps",1000,0,x_max);
  float t_min=T_central-T_range;
  float t_max=T_central+T_range;
  cout<<"pereppe "<<t_min<<"  "<<t_max<<endl;
  
  t_input->Project("h_amps","amp",Form("%f<=temp&&temp<%f&&((amp<40)||(amp>40&&(51.2<time&&time<51.63)))",t_min,t_max));
  //t_input->Project("h_amps","amp","cat==0||cat==1||cat==2||cat==3||cat==6");
  RooDataHist *ds = new RooDataHist("ds","ds",RooArgSet(x),Import(*h_amps)) ;
 
  // D o  F i t 
  // ------------------------
  
  RooPlot* xframe2 = x.frame(Title("Data")) ;
  ds->plotOn(xframe2,DataError(RooAbsData::SumW2)) ;
  RooFitResult* fit_results = model.fitTo(*ds,Save(),Extended(kFALSE),SumW2Error(kTRUE));
  
  
  TH2 *h_correlation = fit_results->correlationHist();
  //model.plotOn(xframe2,Components("PDF_bkg"),LineStyle(kDashed),LineColor(1)) ;
  model.plotOn(xframe2,Components("PDF_1"),LineStyle(kDashed),LineColor(2)) ;
  model.plotOn(xframe2,Components("PDF_2"),LineStyle(kDashed),LineColor(8)) ;
  model.plotOn(xframe2,Components("PDF_3"),LineStyle(kDashed),LineColor(4)) ;
  model.plotOn(xframe2,Components("PDF_4"),LineStyle(kDashed),LineColor(5)) ;
  model.plotOn(xframe2,Components("PDF_5"),LineStyle(kDashed),LineColor(6)) ;
  model.plotOn(xframe2,Components("PDF_6"),LineStyle(kDashed),LineColor(7)) ;
  model.plotOn(xframe2,Components("PDF_7"),LineStyle(kDashed),LineColor(1)) ;
  model.plotOn(xframe2,Components("PDF_8"),LineStyle(kDashed),LineColor(2)) ;
  model.plotOn(xframe2,Components("PDF_9"),LineStyle(kDashed),LineColor(8)) ;
  model.plotOn(xframe2,Components("PDF_10"),LineStyle(kDashed),LineColor(4)) ;
  model.plotOn(xframe2,Components("PDF_11"),LineStyle(kDashed),LineColor(5)) ;
  model.plotOn(xframe2,Components("PDF_12"),LineStyle(kDashed),LineColor(6)) ;
  model.plotOn(xframe2,Components("PDF_13"),LineStyle(kDashed),LineColor(7)) ;
  model.plotOn(xframe2,Components("PDF_14"),LineStyle(kDashed),LineColor(11)) ;
  model.plotOn(xframe2,Components("PDF_15"),LineStyle(kDashed),LineColor(2)) ;
  
  model.plotOn(xframe2) ;
  //model.paramOn(xframe2);
  
  
  
  // Construct a histogram with the residuals of the data w.r.t. the curve
  RooHist* hresid = xframe2->residHist() ;
  
  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist* hpull = xframe2->pullHist() ;
  
  // Create a new frame to draw the residual distribution and add the distribution to the frame
  RooPlot* xframe3 = x.frame(Title("Residual Distribution")) ;
  xframe3->addPlotable(hresid,"P") ;
  
  // Create a new frame to draw the pull distribution and add the distribution to the frame
  RooPlot* xframe4 = x.frame(Title("Pull Distribution")) ;
  xframe4->addPlotable(hpull,"P") ;
  
  
  
  
  // Draw all frames on a canvas
  fit_results->Print("v");
  //cout<<z_Frac.getVal()+s_Frac.getVal()+d_Frac.getVal()+tr_Frac.getVal()+f_Frac.getVal()+p_Frac.getVal()+e_Frac.getVal()<<endl;
  
  
  
  TF1 *line_0s = new TF1("line_0s","0",-100,1000);line_0s->SetLineColor(8);
  TF1 *line_1ps = new TF1("line_1ps","1",-100,1000);line_1ps->SetLineColor(4);
  TF1 *line_1ns = new TF1("line_1ns","-1",-100,1000);line_1ns->SetLineColor(4);
  TF1 *line_2ps = new TF1("line_2ps","2",-100,1000);line_2ps->SetLineColor(kOrange-2);
  TF1 *line_2ns = new TF1("line_2ns","-2",-100,1000);line_2ns->SetLineColor(kOrange-2);
  TF1 *line_3ps = new TF1("line_3ps","3",-100,1000);line_3ps->SetLineColor(2);
  TF1 *line_3ns = new TF1("line_3ns","-3",-100,1000);line_3ns->SetLineColor(2);
  
  TCanvas* c_Fit = new TCanvas("Fit results","Fit results",800,400) ;
  c_Fit->Divide(2,2) ;
  c_Fit->cd(4) ; gPad->SetLeftMargin(0.15) ; /*xframe->GetYaxis()->SetTitleOffset(1.6)*/ ; h_correlation->Draw("colz:text");
  c_Fit->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe2->GetYaxis()->SetTitleOffset(1.6) ; xframe2->Draw() ; 
  c_Fit->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe2->GetYaxis()->SetTitleOffset(1.6) ; gPad->SetLogy(1); xframe2->Draw() ;
  c_Fit->cd(3) ; gPad->SetLeftMargin(0.15) ; xframe4->GetYaxis()->SetTitleOffset(1.6) ; gPad->SetGridy(); xframe4->GetYaxis()->SetRangeUser(-5,5); xframe4->Draw() ;
  line_3ns->Draw("same");line_3ps->Draw("same");
  line_2ns->Draw("same");line_2ps->Draw("same");
  line_1ns->Draw("same");line_1ps->Draw("same");
  line_0s->Draw("same");
  model.paramOn(xframe2);
  TCanvas* c_2 = new TCanvas("Fit results 2","Fit results 2",800,400) ;
  c_2->Divide(1,2);
  c_2->cd(1); gPad->SetLeftMargin(0.15) ; xframe2->GetYaxis()->SetTitleOffset(1.6) ; /*gPad->SetLogy(1)*/ ; xframe2->Draw() ;
  c_2->cd(2); gPad->SetLeftMargin(0.15) ; xframe4->GetYaxis()->SetTitleOffset(1.6) ; gPad->SetGridy(); xframe4->GetYaxis()->SetRangeUser(-5,5); xframe4->Draw() ;
  line_3ns->Draw("same");line_3ps->Draw("same");
  line_2ns->Draw("same");line_2ps->Draw("same");
  line_1ns->Draw("same");line_1ps->Draw("same");
  line_0s->Draw("same");
  TCanvas* c_3 = new TCanvas("Fit results 3","Fit results 3",800,400) ;
  c_3->Divide(1,2);
  gStyle->SetOptFit(1111);
  c_3->cd(1); gPad->SetLeftMargin(0.15) ; xframe2->GetYaxis()->SetTitleOffset(1.6) ; gPad->SetLogy(1); xframe2->Draw() ;
  c_3->cd(2); gPad->SetLeftMargin(0.15) ; xframe4->GetYaxis()->SetTitleOffset(1.6) ; gPad->SetGridy(); xframe4->GetYaxis()->SetRangeUser(-5,5); xframe4->Draw() ;
  line_3ns->Draw("same");line_3ps->Draw("same");
  line_2ns->Draw("same");line_2ps->Draw("same");
  line_1ns->Draw("same");line_1ps->Draw("same");
  line_0s->Draw("same");
  
  
  
  
  double fracs[15];
  double err_fracs[15];
  double fractions[15];
  double err_fractions[15];
  double n[15];
  double err_n[15];
  for(int g=0;g<15;g++){
    n[g]=g+1;
    err_n[g]=0;
  }
  

  fracs[0]=Frac_1.getVal();
  fracs[1]=Frac_2.getVal();
  fracs[2]=Frac_3.getVal();
  fracs[3]=Frac_4.getVal();
  fracs[4]=Frac_5.getVal();
  fracs[5]=Frac_6.getVal();
  fracs[6]=Frac_7.getVal();
  fracs[7]=Frac_8.getVal();
  fracs[8]=Frac_9.getVal();
  fracs[9]=Frac_10.getVal();
  fracs[10]=Frac_11.getVal();
  fracs[11]=Frac_12.getVal();
  fracs[12]=Frac_13.getVal();
  fracs[13]=Frac_14.getVal();
  //fracs[14]=Frac_15.getVal();
  
  err_fracs[0]=0;//Frac_0.getError();
  err_fracs[1]=0;//Frac_1.getError();
  err_fracs[2]=0;//Frac_2.getError();
  err_fracs[3]=0;//Frac_3.getError();
  err_fracs[4]=0;//Frac_4.getError();
  err_fracs[5]=0;//Frac_5.getError();
  err_fracs[6]=0;//Frac_6.getError();
  err_fracs[7]=0;//Frac_7.getError();
  err_fracs[8]=0;//Frac_8.getError();
  err_fracs[9]=0;//Frac_9.getError();
  err_fracs[10]=0;//Frac_10.getError();
  err_fracs[11]=0;//Frac_11.getError();
  err_fracs[12]=0;//Frac_12.getError();
  err_fracs[13]=0;//Frac_13.getError();
  err_fracs[14]=0;//Frac_14.getError();
  
  
  fractions[0]=fracs[0];
  err_fractions[0]=err_fracs[0];
  double sum_fractions=fractions[0];
  double prod_fractions_n=(1-fractions[0]);
  for(int f=1;f<14;f++){
    fractions[f]=fracs[f];//*fractions[f-1]*(1-fracs[f-1])/fracs[f-1];
    err_fractions[f]=0;//fractions[f]*(err_fracs[f]/fracs[f]+err_fractions[f-1]/fractions[f-1]+(err_fracs[f-1])/(fracs[f-1]*(1-fracs[f-1])));
    sum_fractions=sum_fractions+fractions[f];
    //prod_fractions_n=prod_fractions_n*(1-fracs[f]);
  }
  fractions[14]=0;//1-sum_fractions;
  err_fractions[14]=0;
  
  TCanvas *c_frac = new TCanvas("c_frac","c_frac");
  
  TGraphErrors *FRAC = new TGraphErrors(15,n,fractions,err_n,err_fractions); FRAC->SetTitle("f_{L}");
  
  FRAC->SetMarkerStyle(20);
  FRAC->SetMarkerSize(2);
  FRAC->SetMarkerColor(2);
  
  FRAC->Draw("AP");
  
  /*
    TFile *f_output = new TFile("Output.root","recreate");
    f_output->cd();
    c_Fit->Write();
    c_2->Write();
    c_3->Write();
    c_frac->Write();
    FRAC->Write();
    f_output->cd();
    f_output->Close();
    delete f_output;
    
    
    double sigmas[15];
    double err_sigmas[15];
    
    
    sigmas[0]=sigma_0.getVal();
    sigmas[1]=sigma_1.getVal();
    sigmas[2]=sigma_2.getVal();
    sigmas[3]=sigma_3.getVal();
    sigmas[4]=sigma_4.getVal();
    sigmas[5]=sigma_5.getVal();
    sigmas[6]=sigma_6.getVal();
    sigmas[7]=sigma_7.getVal();
    sigmas[8]=sigma_8.getVal();
    sigmas[9]=sigma_9.getVal();
    sigmas[10]=sigma_10.getVal();
    sigmas[11]=sigma_11.getVal();
    sigmas[12]=sigma_12.getVal();
    sigmas[13]=sigma_13.getVal();
    sigmas[14]=sigma_14.getVal();
    
    err_sigmas[0]=sigma_0.getError();
    err_sigmas[1]=sigma_1.getError();
    err_sigmas[2]=sigma_2.getError();
    err_sigmas[3]=sigma_3.getError();
    err_sigmas[4]=sigma_4.getError();
    err_sigmas[5]=sigma_5.getError();
    err_sigmas[6]=sigma_6.getError();
    err_sigmas[7]=sigma_7.getError();
    err_sigmas[8]=sigma_8.getError();
    err_sigmas[9]=sigma_9.getError();
    err_sigmas[10]=sigma_10.getError();
    err_sigmas[11]=sigma_11.getError();
    err_sigmas[12]=sigma_12.getError();
    err_sigmas[13]=sigma_13.getError();
    err_sigmas[14]=sigma_14.getError();
    
    
    TCanvas *c_sigma = new TCanvas("c_sigma","c_sigma");
    TGraphErrors *SIGMA = new TGraphErrors(15,n,sigmas,err_n,err_sigmas); SIGMA->SetTitle("#sigma_{i}");
    
    SIGMA->SetMarkerStyle(20);
    SIGMA->SetMarkerSize(2);
    SIGMA->SetMarkerColor(2);
    
    SIGMA->Draw("AP");
  */
  
  Results.mean[0]=mean.getVal();
  Results.mean[1]=mean.getError();
  Results.delta[0]=Delta.getVal();
  Results.delta[1]=Delta.getError();
  Results.Av_N_ph[0]=Av_N_ph.getVal();
  Results.Av_N_ph[1]=Av_N_ph.getError();
  Results.CT_prob[0]=q.getVal();
  Results.CT_prob[1]=q.getError();
  Results.s0[0]=s0.getVal();
  Results.s0[1]=s0.getError();
  Results.s1[0]=s1.getVal();
  Results.s1[1]=s1.getError();
  
  TFile *f_output = new TFile(input_basepath+"Output_"+string_out+".root","update");
  f_output->cd();
  f_output->mkdir(Form("Temp_%f",T_central));
  f_output->cd(Form("Temp_%f",T_central));
  c_Fit->Write();
  c_2->Write();
  c_3->Write();
  c_frac->Write();
  FRAC->Write();
  f_output->cd();
  f_output->Close();
  delete f_output;

  return Results;
}




















































void run_scan_fit(){
  Float_t LOW_TEMP_RANGE=23.0;
  Float_t UP_TEMP_RANGE=27.0;
  Int_t N_scan=7;
  float T_range=(UP_TEMP_RANGE-LOW_TEMP_RANGE)/N_scan;

  float T[7];
  float err_T[7];
  float mean[7];
  float delta[7];
  float av_n_ph[7];
  float q_CT[7];
  float s0[7];
  float s1[7];
  float err_mean[7];
  float err_delta[7];
  float err_av_n_ph[7];
  float err_q_CT[7];
  float err_s0[7];
  float err_s1[7];
    
  for(int j=1;j<=N_scan;j++){
    Fit_results Results;
    float t_central=(LOW_TEMP_RANGE+0.5*T_range)+(j-1)*T_range;
    Results=Fit_Laser_Stability(1,t_central,0.5*T_range);

    T[j-1]=t_central;
    err_T[j-1]=0.5*T_range;
    mean[j-1]=Results.mean[0];
    delta[j-1]=Results.delta[0];
    av_n_ph[j-1]=Results.Av_N_ph[0];
    q_CT[j-1]=Results.CT_prob[0];
    s0[j-1]=Results.s0[0];
    s1[j-1]=Results.s1[0];
    err_mean[j-1]=Results.mean[1];
    err_delta[j-1]=Results.delta[1];
    err_av_n_ph[j-1]=Results.Av_N_ph[1];
    err_q_CT[j-1]=Results.CT_prob[1];
    err_s0[j-1]=Results.s0[1];
    err_s1[j-1]=Results.s1[1];
  }

  TGraphErrors *DELTA = new TGraphErrors(7,T,delta,err_T,err_delta);
  DELTA->SetTitle("#Delta");
  DELTA->SetName("Separation between two peaks");
  DELTA->SetMarkerStyle(20);
  DELTA->SetMarkerSize(2);
  DELTA->SetMarkerColor(4);
  DELTA->GetYaxis()->SetTitle("#Delta [ADC counts]");
  DELTA->GetXaxis()->SetTitle("Temperature [C°]");

  TGraphErrors *MEAN = new TGraphErrors(7,T,mean,err_T,err_mean);
  MEAN->SetTitle("#mu_{0}");
  MEAN->SetName("position of first peak");
  MEAN->SetMarkerStyle(20);
  MEAN->SetMarkerSize(2);
  MEAN->SetMarkerColor(4);
  MEAN->GetYaxis()->SetTitle("#mu_{0} [ADC counts]");
  MEAN->GetXaxis()->SetTitle("Temperature [C°]");
  
  TGraphErrors *AV_N_PH = new TGraphErrors(7,T,av_n_ph,err_T,err_av_n_ph);
  AV_N_PH->SetTitle("n_{#gamma}");
  AV_N_PH->SetName("position of first peak");
  AV_N_PH->SetMarkerStyle(20);
  AV_N_PH->SetMarkerSize(2);
  AV_N_PH->SetMarkerColor(4);
  AV_N_PH->GetYaxis()->SetTitle("n_{#gamma}");
  AV_N_PH->GetXaxis()->SetTitle("Temperature [C°]");
  
  TGraphErrors *Q_CT = new TGraphErrors(7,T,q_CT,err_T,err_q_CT);
  Q_CT->SetTitle("q_{CT}");
  Q_CT->SetName("CT probability");
  Q_CT->SetMarkerStyle(20);
  Q_CT->SetMarkerSize(2);
  Q_CT->SetMarkerColor(4);
  Q_CT->GetYaxis()->SetTitle("q_{CT}");
  Q_CT->GetXaxis()->SetTitle("Temperature [C°]");
  
  TGraphErrors *S0 = new TGraphErrors(7,T,s0,err_T,err_s0);
  S0->SetTitle("s_{0}");
  S0->SetName("#sigma=s_{0}+n#cdot s_{1}");
  S0->SetMarkerStyle(20);
  S0->SetMarkerSize(2);
  S0->SetMarkerColor(4);
  S0->GetYaxis()->SetTitle("s_{0}");
  S0->GetXaxis()->SetTitle("Temperature [C°]");
  
  TGraphErrors *S1 = new TGraphErrors(7,T,s1,err_T,err_s1);
  S1->SetTitle("s_{1}");
  S1->SetName("#sigma=s_{0}+n#cdot s_{1}");
  S1->SetMarkerStyle(20);
  S1->SetMarkerSize(2);
  S1->SetMarkerColor(4);
  S1->GetYaxis()->SetTitle("s_{1}");
  S1->GetXaxis()->SetTitle("Temperature [C°]");
  

  TCanvas *c = new TCanvas("c","c");
  c->Divide(2,3);
  c->cd(1);
  DELTA->Draw("AP");
  c->cd(3);
  MEAN->Draw("AP");
  c->cd(5);
  AV_N_PH->Draw("AP");
  c->cd(2);
  Q_CT->Draw("AP");
  c->cd(4);
  S0->Draw("AP");
  c->cd(6);
  S1->Draw("AP");
  
  
  
}














void make_plots_Scan(){
 
  TFile *_file0 = new TFile("Output_Free_fractions_T50_7.root");
  TFile *_file1 = new TFile("Output_Free_fractions_T50_6.root");
  TFile *_file2 = new TFile("Output_Free_fractions_T50_5.root");
  TFile *_file3 = new TFile("Output_Free_fractions_T50_4.root");
  TFile *_file4 = new TFile("Output_Free_fractions_T50_3.root");
  TFile *_file5 = new TFile("Output_Free_fractions_T50_2.root");
  TFile *_file6 = new TFile("Output_Free_fractions_T50_1.root");


  TGraphErrors *tge_res_1 = (TGraphErrors*)_file0->Get("Graph");
  TGraphErrors *tge_res_2 = (TGraphErrors*)_file1->Get("Graph");
  TGraphErrors *tge_res_3 = (TGraphErrors*)_file2->Get("Graph");
  TGraphErrors *tge_res_4 = (TGraphErrors*)_file3->Get("Graph");
  TGraphErrors *tge_res_5 = (TGraphErrors*)_file4->Get("Graph");
  TGraphErrors *tge_res_6 = (TGraphErrors*)_file5->Get("Graph");
  TGraphErrors *tge_res_7 = (TGraphErrors*)_file6->Get("Graph");
  tge_res_1->SetLineColor(1);
  tge_res_2->SetLineColor(2);
  tge_res_3->SetLineColor(4);
  tge_res_4->SetLineColor(8);
  tge_res_5->SetLineColor(7);
  tge_res_6->SetLineColor(6);
  tge_res_7->SetLineColor(5);
  tge_res_1->SetFillColor(0);
  tge_res_2->SetFillColor(0);
  tge_res_3->SetFillColor(0);
  tge_res_4->SetFillColor(0);
  tge_res_5->SetFillColor(0);
  tge_res_6->SetFillColor(0);
  tge_res_7->SetFillColor(0);
  tge_res_1->SetMarkerColor(1);
  tge_res_2->SetMarkerColor(2);
  tge_res_3->SetMarkerColor(4);
  tge_res_4->SetMarkerColor(8);
  tge_res_5->SetMarkerColor(7);
  tge_res_6->SetMarkerColor(6);
  tge_res_7->SetMarkerColor(5);
  tge_res_1->SetTitle("7");
  tge_res_2->SetTitle("6");
  tge_res_3->SetTitle("5");
  tge_res_4->SetTitle("4");
  tge_res_5->SetTitle("3");
  tge_res_6->SetTitle("2");
  tge_res_7->SetTitle("1");
  TCanvas *c_res = new TCanvas("c_res","c_res");
  c_res->cd();
  TMultiGraph *tmg_res = new TMultiGraph();
  tmg_res->Add(tge_res_1);
  tmg_res->Add(tge_res_2);
  tmg_res->Add(tge_res_3);
  tmg_res->Add(tge_res_4);
  tmg_res->Add(tge_res_5);
  tmg_res->Add(tge_res_6);
  tmg_res->Add(tge_res_7);
  tmg_res->Draw("AP");
  tmg_res->GetXaxis()->SetTitle("number of photons");
  tmg_res->GetYaxis()->SetTitle("fraction");
  //tmg_res->GetYaxis()->SetRangeUser(0.05,0.08);
  c_res->SetGridy();
  gPad->BuildLegend();
  gPad->Update();

  TFile *f = new TFile("output_Stability_scan.root","recreate");
  f->cd();
  
  c_res->Write();
 
  tmg_res->Write();
  f->cd();
  f->Close();
  delete f;

  
}


