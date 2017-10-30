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

void Fit_head_2(int fit_model){
  // S e t u p   m o d e l 
  // ---------------------
  
  // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
  //RooRealVar x("Amplitude_212","Amplitude_212",-150,300) ;
  TString string_out="";
  if(fit_model==1){string_out="Poisson_NoOCT";}
  else if (fit_model==2){string_out="Poisson_OCT";}
  else if (fit_model==3){string_out="Free_fractions";}
  
  gROOT->ProcessLine(".x ./myRooPdfs/RooExpGauss.cxx+") ;
  gROOT->ProcessLine(".x ./myRooPdfs/RooAsymGauss.cxx+") ;
  
  int x_max=1000;
  RooRealVar x("amps","amps",0,x_max) ;//
  RooRealVar t_min("t_min","t_min",289,285,292) ;//
  RooRealVar mean("mean","zero_value",44.5);//,-50,60) ;
  RooRealVar Delta("Delta","Delta",40,30,120);
  //RooRealVar Delta("Delta","Delta",50,20,120);
  //RooRealVar mean("mean","zero_value",20,-50,50) ;
  RooRealVar Av_N_ph("Av_N_ph","Average number of photons",7,0.,20);
  
  RooRealVar N_tot("N_tot","N Signal events",95608,50000,98000); N_tot.setConstant(kTRUE);
  RooRealVar N_bkg("N_bkg","Number of background events",2000,500,5000) ;
  
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
  */
  
  RooRealVar sigma_0("sigma_0","#sigma_{0}",10,1,20) ;
  RooRealVar s0("s0","s0",5.8559);//,5.855-0.2139,5.855+0.2139);//
  RooRealVar s1("s1","s1",0.93189);//,0.9318-0.0385,0.9318+0.0385);
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
  
  
  
  
  
  RooRealVar sigma_bkg("sigma_bkg","sigma_bkg",500,50,1000);
  RooRealVar mean_bkg("mean_bkg","mean_bkg",300,200,500);
  RooGaussian PDF_bkg("PDF_bkg","PDF background",x,mean_bkg,sigma_bkg) ;
  RooRealVar  Frac_sig("Frac_sig"  ,"fraction of signal events", 0.9 ,0.8,1);
  
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
  RooGaussian PDF_11("PDF_11","gaussian PDF background",x,mean_11,sigma_11) ;
  RooGaussian PDF_12("PDF_12","gaussian PDF background",x,mean_12,sigma_12) ;
  RooGaussian PDF_13("PDF_13","gaussian PDF background",x,mean_13,sigma_13) ;
  RooGaussian PDF_14("PDF_14","gaussian PDF background",x,mean_14,sigma_14) ;
  RooRealVar alpha_15("alpha_15","Observed number of photons",-0.5,-10,0);
  RooExpGauss PDF_15("PDF_15","gaussian PDF background",x,mean_15,sigma_15,alpha_15) ;
  
  
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
    RooPoisson  Frac_11( "Frac_11","fraction of 11 ph signal events"  ,N_Obs_11,Av_N_ph);
    RooPoisson  Frac_12( "Frac_12","fraction of 12 ph signal events"  ,N_Obs_12,Av_N_ph);
    RooPoisson  Frac_13( "Frac_13","fraction of 13 ph signal events"  ,N_Obs_13,Av_N_ph);
    RooPoisson  Frac_14( "Frac_14","fraction of 13 ph signal events"  ,N_Obs_14,Av_N_ph);
    RooRealVar  Frac_15( "Frac_15"  ,"Frac_{15}", 0.006,0.0001,0.5);
    
    RooArgList pdfList(PDF_1,PDF_2,PDF_3,PDF_4,PDF_5,PDF_6,PDF_7,PDF_8);
    pdfList.add(RooArgList(PDF_9,PDF_10,PDF_11,PDF_12,PDF_13,PDF_14,PDF_15));
    RooArgList fracList(Frac_1,Frac_2,Frac_3,Frac_4,Frac_5,Frac_6,Frac_7,Frac_8);
    fracList.add(RooArgList(Frac_9,Frac_10,Frac_11,Frac_12,Frac_13,Frac_14,Frac_15));
    
    
    RooAddPdf model("model","model",pdfList,fracList);
    
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
    
    RooRealVar q("q","q_{CT}", 0.1,0.0,1.0);
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
    
      
    RooArgList pdfList(PDF_1,PDF_2,PDF_3,PDF_4,PDF_5,PDF_6,PDF_7,PDF_8);
    pdfList.add(RooArgList(PDF_9,PDF_10,PDF_11,PDF_12,PDF_13,PDF_14,PDF_15));
    RooArgList fracList(Frac_1,Frac_2,Frac_3,Frac_4,Frac_5,Frac_6,Frac_7,Frac_8);
    fracList.add(RooArgList(Frac_9,Frac_10,Frac_11,Frac_12,Frac_13,Frac_14,Frac_15));
    
    
    RooAddPdf model("model","model",pdfList,fracList);
    
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

    RooArgList  pdfList(PDF_1,PDF_2,PDF_3,PDF_4,PDF_5,PDF_6,PDF_7,PDF_8);
    pdfList.add(RooArgList(PDF_9,PDF_10,PDF_11,PDF_12,PDF_13,PDF_14,PDF_15));
    RooArgList fracList(Frac_1,Frac_2,Frac_3,Frac_4,Frac_5,Frac_6,Frac_7,Frac_8);
    fracList.add(RooArgList(Frac_9,Frac_10,Frac_11,Frac_12,Frac_13,Frac_14,Frac_15));
    
    
    RooAddPdf model("model","model",pdfList,fracList);
  }				    
 
  
  
  // I m p o r t    d a t a
  // ------------------------
  TFile *f_input = new TFile("scan-att-Vb30.5-thr-20-T80-F15000.root");
  TTree *t_input = (TTree*)f_input->Get("times");
  TH1D *h_amps = new TH1D("h_amps","h_amps",1000,0,x_max);
  t_input->Project("h_amps","amps[0]");
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
  
  
  TFile *f_output = new TFile("Output_"+string_out+".root","recreate");
  f_output->cd();
  c_Fit->Write();
  c_2->Write();
  c_3->Write();
  c_frac->Write();
  FRAC->Write();
  f_output->cd();
  f_output->Close();
  delete f_output;
  
}




























































void Fit_test(int fit_model){
  // S e t u p   m o d e l 
  // ---------------------
  
  // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
  //RooRealVar x("Amplitude_212","Amplitude_212",-150,300) ;
  TString string_out="";
  if(fit_model==1){string_out="Poisson_NoOCT";}
  else if (fit_model==2){string_out="Poisson_OCT";}
  else if (fit_model==3){string_out="Free_fractions";}
  
  gROOT->ProcessLine(".x ./myRooPdfs/RooExpGauss.cxx+") ;
  gROOT->ProcessLine(".x ./myRooPdfs/RooAsymGauss.cxx+") ;
  
  int x_max=200;
  RooRealVar x("amps","amps",25,x_max) ;//
  RooRealVar mean("mean","zero_value",44.5,40,60) ;
  RooRealVar Delta("Delta","Delta",42.92) ;//,10,120);
  RooRealVar Av_N_ph("Av_N_ph","Average number of photons",7,0.,20);
  
  RooRealVar N_tot("N_tot","N Signal events",95608,50000,98000); N_tot.setConstant(kTRUE);
  RooRealVar N_bkg("N_bkg","Number of background events",2000,500,5000) ;
  
  RooRealVar N_Obs_0("N_Obs_0","Observed number of photons",0);
  RooRealVar N_Obs_1("N_Obs_1","Observed number of photons",1);
  RooRealVar N_Obs_2("N_Obs_2","Observed number of photons",2);
  RooRealVar N_Obs_3("N_Obs_3","Observed number of photons",3);
  RooRealVar N_Obs_4("N_Obs_4","Observed number of photons",4);
 
 
    RooRealVar sigma_0("sigma_0","#sigma_{0}",10,1,20) ;
    RooRealVar sigma_1("sigma_1","#sigma_{1}",10,1,20) ;
    RooRealVar sigma_2("sigma_2","#sigma_{2}",5,1,25) ;
    RooRealVar sigma_2_CT("sigma_2_CT","#sigma_{2}",5,1,20) ;
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

  /*  
  RooRealVar s0("s0","s0",5.8559,5.855-0.2139,5.855+0.2139);//
  RooRealVar s1("s1","s1",0.93189,0.9318-0.0385,0.9318+0.0385);
  RooFormulaVar sigma_1("sigma_1","sigma_1","s0+0*s1",RooArgList(s0,s1));
  RooFormulaVar sigma_2("sigma_2","sigma_2","s0+1*s1",RooArgList(s0,s1));
  RooRealVar sigma_2_CT("sigma_2_CT","#sigma_{2}",10,1,20) ;
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
  */
 
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
    
  
  
  
  
  RooRealVar alpha("alpha","Observed number of photons",-0.5,-10,0);
 
  RooGaussian PDF_1("PDF_1","gaussian PDF background",x,mean_1,sigma_1) ;
  
  //RooGaussian PDF_2_C("PDF_2_C","gaussian PDF background",x,mean_2,sigma_2) ;
  RooExpGauss PDF_2_C("PDF_2_C","gaussian PDF background",x,mean_2,sigma_2,alpha) ;
  RooRealVar Delta_2_CT("Delta_2_CT","Delta_2_CT",20,10,30);
  RooFormulaVar mean_2_CT("mean_2_CT","mean_2_CT","mean_2-Delta_2_CT",RooArgList(mean_2,Delta_2_CT));
  RooGaussian PDF_2_CT("PDF_2_CT","gaussian PDF background",x,mean_2_CT,sigma_2_CT) ;
  RooRealVar  Frac_2_CT("Frac_2_CT"    ,"Frac_{2}", 0.80,0,1);
  RooAddPdf PDF_2("PDF_2","PDF_2",RooArgList(PDF_2_C,PDF_2_CT),Frac_2_CT);
  
  RooExpGauss PDF_3("PDF_3","gaussian PDF background",x,mean_3,sigma_3,alpha) ;
  if(fit_model){
    
    
    ////////////////////////////////////////////////////////////////////
    ///////////////////////////// FREE FRACTIONS////////////////////////
    ////////////////////////////////////////////////////////////////////
    RooRealVar  Frac_1("Frac_1"    ,"Frac_{1}", 0.20,0,0.8);
    RooRealVar  Frac_2("Frac_2"    ,"Frac_{2}", 0.10,0,1);

    RooRealVar  Frac_3("Frac_3"    ,"Frac_{3}", 0.07,0,0.5);
    RooRealVar  Frac_4("Frac_4"    ,"Frac_{4}", 0.04,0,0.5);

    RooArgList  pdfList(PDF_1,PDF_2);
    RooArgList fracList(Frac_1);
    
    
    RooAddPdf model("model","model",pdfList,fracList,kTRUE);
  }				    
 
  
  
  // I m p o r t    d a t a
  // ------------------------
  TFile *f_input = new TFile("scan-att-Vb30.5-thr-20-T80-F15000.root");
  TTree *t_input = (TTree*)f_input->Get("times");
  TH1D *h_amps = new TH1D("h_amps","h_amps",500,0,x_max);
  t_input->Project("h_amps","amps[0]");
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
  model.plotOn(xframe2,Components("PDF_2_CT"),LineStyle(kDashed),LineColor(7)) ;
  model.plotOn(xframe2,Components("PDF_2_C"),LineStyle(kDashed),LineColor(1)) ;
  //model.plotOn(xframe2,Components("PDF_3"),LineStyle(kDashed),LineColor(4)) ;
 
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
  
  
  
  
}



