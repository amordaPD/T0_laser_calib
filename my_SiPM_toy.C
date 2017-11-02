////////Code to study SiPM behavior


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



void my_SiPM_toy(){

  float digit_to_time_factor= (1.0/3.2);
  float av_norm_wf[1024];
  for(int hh =0; hh<1024;hh++){
    av_norm_wf[hh]=0;
  }
  Int_t counter_good_wf=0;
  TH1F *h_av_norm_wf= new TH1F("h_av_norm_wf","h_av_norm_wf",1024,0,1024*(digit_to_time_factor));
  
  TH1D *h_t0 = new TH1D("h_t0","h_t0",1024,0,1024);
  TH2F *toy_WF = new TH2F("toy_WF","toy_WF",1024,0,1024*digit_to_time_factor,1000,-5,250);
  TH2F *toy_WF_norm = new TH2F("toy_WF_norm","toy_WF_norm",1024,0,1024*digit_to_time_factor,250,-0.2,1.1);
  TH2F *sig_sub_WF = new TH2F("sig_sub_WF","sig_sub_WF",800 /*300*/,-100*digit_to_time_factor,700*digit_to_time_factor,500,-5,500);
  TH2F *sig_sub_WF_APCT = new TH2F("sig_sub_WF_APCT","sig_sub_WF_APCT",700 /*300*/,0*digit_to_time_factor,700*digit_to_time_factor,500,-5,300);
  
  std::vector<float> WF_sub_shape;
  std::vector<float> WF_shape;
  std::vector<float> TR_shape;

  Float_t wf_integral=0;
  Int_t N_peaks=-9;
  Float_t Amp_main_peak=-99;

  TFile *ff = new TFile("my_SiPM_output.root","recreate");
  
  vector<float> peak_height;
  vector<float> peak_pos;

  TTree *tree_wf = new TTree("tree_wf","waveform analysis");
  tree_wf->Branch("wf_integral",&wf_integral,"wf_integral/F");
  tree_wf->Branch("N_peaks",&N_peaks,"N_peaks/I");
  tree_wf->Branch("Amp_main_peak",&Amp_main_peak,"Amp_main_peak/F");
  tree_wf->Branch("peak_height",&peak_height,"peak_height[N_peaks]/F");
  tree_wf->Branch("peak_pos",&peak_pos,"peak_pos[N_peaks]/F");
  

  
  /*
    TFile *f_input= new TFile("Toy_wf.root");
    TTree *t_input = (TTree*)f_input->Get("tree_wf");
    t_input->SetBranchAddress("wf",&(WF_shape[0]));
  */


  TFile *f_template = new TFile("waveforms-thr--20-T50-F2000-waves0-2_sig_wf_Template.root");
  TH1D *h_SigWfTemp = (TH1D*)f_template->Get("h_av_norm_wf");
  Float_t max_bin = h_SigWfTemp->GetMaximumBin();


  /*
  vector<double> SigWfTemp(800);
  for(int h=0;h<800;h++){
    SigWfTemp[h]=h_SigWfTemp->GetBinContent(h);
  }
  double rho = 1.0; //default value
  TKDE * kde = new TKDE(800, &SigWfTemp[0],50 ,300, "KernelType:Gaussian", rho);
  */

   
  float _WF_shape[8][1024];
  TFile *f_input = new TFile("dati/waveforms-thr--20-T50-F2000-waves0-2.root");
  TTree *t_input = (TTree*)f_input->Get("DT5743");
  t_input->SetBranchAddress("adcval",&(_WF_shape[0]));


  Int_t nentries = t_input->GetEntries();
  cout<<"input wf : "<<nentries<<endl;
  for(Int_t i=0; i<nentries;i++) {//£
    Float_t integral_WF=0;
    Int_t N_extra_peaks_WF=0;
    ////////////////////////////////////////////////////////////////
    /// INITIALIZE TRIGGER AND SiPM WAVEFUNCTIONS///////////////////
    ////////////////////////////////////////////////////////////////
    WF_shape.clear();
    TR_shape.clear();
    WF_sub_shape.clear();
    
    t_input->GetEvent(i);
    for(int j=0;j<1024;j++){
      WF_shape.push_back(-_WF_shape[1][j]);
      TR_shape.push_back(_WF_shape[0][j]);
    }


    
    ////////////////////////////////////////////////////////////////
    ////////COMPUTE TRIGGER TIME WITH CONSTANT FRACTION METHOD//////
    ////////////////////////////////////////////////////////////////
    
    double trigger_time=-10;
    double thr=0.5;
    int samples=100;
    // remove offset to trigger signal
    if (TR_shape.size()>samples) {
      float avg = accumulate(TR_shape.begin(), TR_shape.begin()+samples,0);
      avg = avg/double(samples);
      float minval = *(min_element(TR_shape.begin(), TR_shape.end()));
      //     cout << avg ;
      for (unsigned int ii=0; ii<TR_shape.size(); ++ii) TR_shape[ii] -= avg;
      float maxval1 = accumulate(TR_shape.end()-120, TR_shape.end()-20,0)/100.;
      float maxval2 = accumulate(TR_shape.begin(), TR_shape.begin()+100,0)/100.;
      float maxval = maxval1*thr;
      if (maxval2>maxval1) maxval = maxval2*thr;
      unsigned int indx=0;
      for (unsigned int ii=5; ii<TR_shape.size()-5; ++ii) {
	indx = TR_shape.size()-ii+1;
	if (TR_shape[TR_shape.size()-ii] < maxval) break;
      }
      double cft = -1;
      if (indx>1) {
	cft =  ((indx-1)/(maxval-TR_shape[indx-1]) + (indx)/(TR_shape[indx]-maxval))/(1./(maxval-TR_shape[indx-1]) + 1./(TR_shape[indx]-maxval));
      }
      trigger_time =cft;
    }else {
      trigger_time= -1.;
    }

    /////////////////////////////////////////
    /// BASELINE ESTIMATION ////////////////
    ////////////////////////////////////////
    int n_samples = 150;
    float avg_wf = accumulate(WF_shape.begin()+20, WF_shape.begin()+20+n_samples,0);
    avg_wf=avg_wf/n_samples;
    
    int t0 = -99;
    bool _find_t0 = true;
    int i_max =-99;
    float max=-99;

    /////////////////////////////
    /// WAVEFORM ESTIMATION  ////
    /////////////////////////////
    
    for(int j=0;j<1024;j++){//§
      
      /////BASELINE SUBTRACTION
      WF_shape[j]-=avg_wf;

      integral_WF+=WF_shape[j];
      /////SEARCH FOR WF MAXIMUM AT NOMINAL SIGNAL TIME WINDOW
      //////// TRY WITH THAT https://root.cern.ch/root/html/tutorials/math/exampleTKDE.C.html
      if(TMath::Abs((float(j)-trigger_time+200)*(digit_to_time_factor)-70)<10){
	if(_find_t0&&WF_shape[j]>25&&(
		      WF_shape[j]>WF_shape[j-2]&&
		      WF_shape[j]>WF_shape[j+2]&&
		      WF_shape[j]>WF_shape[j-1]&&
		      WF_shape[j]>WF_shape[j+1])) {t0=j; _find_t0=false;}
      }

    }// end   for(int j=0;j<1024;j++){ //§
    wf_integral=integral_WF/1000000.0;
    if(t0>0) {Amp_main_peak=WF_shape[t0];}else{Amp_main_peak=t0;}
    
    ////if he founds a maximum in the right position will fill the scope histogram and the subtracted wavefunction
    if(!_find_t0){//$
      for(int j=0;j<800;j++){//§§
	sig_sub_WF->Fill((-100+j)*(digit_to_time_factor),WF_shape[t0-100+j]-WF_shape[t0]*(h_SigWfTemp->GetBinContent(max_bin-100+j)));
	if(j>=100) WF_sub_shape.push_back(WF_shape[t0-100+j]-WF_shape[t0]*(h_SigWfTemp->GetBinContent(max_bin-100+j)));
      }//end for(int j=0;j<1024;j++){//§§

      //look for extra maxima
      for(int f=4;f<WF_sub_shape.size();f++){//@
	if(WF_sub_shape[f]>WF_sub_shape[f-1]&&
	   WF_sub_shape[f]>WF_sub_shape[f-2]&&
	   WF_sub_shape[f]>WF_sub_shape[f+1]&&
	   WF_sub_shape[f]>WF_sub_shape[f+2]&&
	   WF_sub_shape[f]>WF_sub_shape[f-3]&&
	   WF_sub_shape[f]>WF_sub_shape[f+3]&&
	   WF_sub_shape[f]>WF_sub_shape[f-4]&&
	   WF_sub_shape[f]>WF_sub_shape[f+4]&&
	   WF_sub_shape[f]>25&&
	   (WF_sub_shape[f]-WF_sub_shape[f-2]>10)
	   ) {

	  N_extra_peaks_WF++;
	  sig_sub_WF_APCT->Fill(f*(digit_to_time_factor),WF_sub_shape[f]);
	  //break;
	  //maxima_apct.push_back(WF_sub_shape[f]);
	  //pos_maxima_apct.push_back(f*(digit_to_time_factor));
	}
	
      }//end for(int f=0;f<WF_sub_shape.size();f++){//@

      if(N_extra_peaks_WF>0) {N_peaks=N_extra_peaks_WF;}else {N_peaks=-99;}
      
    }//end if(!_find_t0){//$
    tree_wf->Fill();
  }//end for (Int_t i=0; i<nentries;i++) {//£


  sig_sub_WF->Draw("colz");
  gPad->SetLogz();
  TCanvas *c_APCT = new TCanvas("c_APCT","c_APCT");
  c_APCT->cd();
  sig_sub_WF_APCT->Draw("colz");
  ff->cd();
  tree_wf->Write();
  ff->Close();
  delete ff;
}


void my_SiPM_produce_signal_wf_template(){

  float digit_to_time_factor= (1.0/3.2);
  float av_norm_wf[1024];
  for(int hh =0; hh<1024;hh++){
    av_norm_wf[hh]=0;
  }
  Int_t counter_good_wf=0;

  std::vector<float> WF_shape;
  std::vector<float> TR_shape;


  TH1F *h_av_norm_wf= new TH1F("h_av_norm_wf","h_av_norm_wf",1024,0,1024*(digit_to_time_factor));
  TH2F *toy_WF = new TH2F("toy_WF","toy_WF",1024,0,1024*digit_to_time_factor,1000,-5,250);
  TH2F *toy_WF_norm = new TH2F("toy_WF_norm","toy_WF_norm",1024,0,1024*digit_to_time_factor,250,-0.2,1.1);
  
  float _WF_shape[8][1024];
  TFile *f_input = new TFile("dati/waveforms-thr--20-T50-F2000-waves0-2.root");
  TTree *t_input = (TTree*)f_input->Get("DT5743");
  t_input->SetBranchAddress("adcval",&(_WF_shape[0]));


  Int_t nentries = t_input->GetEntries();
  cout<<"input wf : "<<nentries<<endl;
  for(Int_t i=0; i<nentries;i++) {

    ////////////////////////////////////////////////////////////////
    /// INITIALIZE TRIGGER AND SiPM WAVEFUNCTIONS///////////////////
    ////////////////////////////////////////////////////////////////
    WF_shape.clear();
    TR_shape.clear();
    t_input->GetEvent(i);
    float max=-99;
    int i_max=-10;
    for(int j=0;j<1024;j++){
      WF_shape.push_back(-_WF_shape[1][j]);
      TR_shape.push_back(_WF_shape[0][j]);
    }


    
    ////////////////////////////////////////////////////////////////
    ////////COMPUTE TRIGGER TIME WITH CONSTANT FRACTION METHOD//////
    ////////////////////////////////////////////////////////////////
    
    double trigger_time=-10;
    double thr=0.5;
    int samples=100;
    // remove offset to trigger signal
    if (TR_shape.size()>samples) {
      float avg = accumulate(TR_shape.begin(), TR_shape.begin()+samples,0);
      avg = avg/double(samples);
      float minval = *(min_element(TR_shape.begin(), TR_shape.end()));
      //     cout << avg ;
      for (unsigned int ii=0; ii<TR_shape.size(); ++ii) TR_shape[ii] -= avg;
      float maxval1 = accumulate(TR_shape.end()-120, TR_shape.end()-20,0)/100.;
      float maxval2 = accumulate(TR_shape.begin(), TR_shape.begin()+100,0)/100.;
      float maxval = maxval1*thr;
      if (maxval2>maxval1) maxval = maxval2*thr;
      unsigned int indx=0;
      for (unsigned int ii=5; ii<TR_shape.size()-5; ++ii) {
	indx = TR_shape.size()-ii+1;
	if (TR_shape[TR_shape.size()-ii] < maxval) break;
      }
      double cft = -1;
      if (indx>1) {
	cft =  ((indx-1)/(maxval-TR_shape[indx-1]) + (indx)/(TR_shape[indx]-maxval))/(1./(maxval-TR_shape[indx-1]) + 1./(TR_shape[indx]-maxval));
      }
      trigger_time =cft;
    }else {
      trigger_time= -1.;
    }

    /////////////////////////////////////////
    /// BASELINE ESTIMATION ////////////////
    ////////////////////////////////////////
    int n_samples = 150;
    float avg_wf = accumulate(WF_shape.begin()+20, WF_shape.begin()+20+n_samples,0);
    avg_wf=avg_wf/n_samples;
    
    int t0 = 0;
    bool _find_t0 = true;
    int i_max =-99;
    float max=-99;
    
    
    for(int j=0;j<1024;j++){
      WF_shape[j]-=avg_wf;
      if(WF_shape[j]>max) {i_max=j; max=WF_shape[j];}
      
    }

    /////////////////////////////////////////
    ///COMPUTATION OF NORMALIZED WAVEFORMS///
    /////////////////////////////////////////

    //check if the absolute maximum is in the expected window around the nominal arrival time wrt trigger signal
    if(TMath::Abs((float(i_max)-trigger_time+200)*(digit_to_time_factor)-70)<2) {counter_good_wf++;}
    //fill the scope histogram for all the events (toy_WF) and only the good ones (toy_wf_norm)
    for(int j=0;j<1024;j++){     
      if(TMath::Abs((float(i_max)-trigger_time+200)*(digit_to_time_factor)-70)<2) {
	toy_WF_norm->Fill((float(j)-trigger_time+200)*(digit_to_time_factor),WF_shape[j]/max);
	if(i_max-200+j<1024) av_norm_wf[j]+=WF_shape[i_max-200+j]/max;
      }
      toy_WF->Fill((float(j)-trigger_time+200)*(digit_to_time_factor),WF_shape[j]);
    }

  }//end for (Int_t i=0; i<nentries;i++) {
  
  for(int kk=0;kk<1024;kk++){
    h_av_norm_wf->SetBinContent(kk+1,av_norm_wf[kk]/counter_good_wf);
  }

  cout<<"number of good wf : "<<counter_good_wf<<endl;
  
  TCanvas *c = new TCanvas("c","c");
  c->cd();
  c->Divide(2,2);
  c->cd(1);
  h_av_norm_wf->Draw();
  c->cd(2);
  toy_WF_norm->DrawNormalized("colz");
  //h_av_norm_wf->Draw("same");
  gPad->SetLogz();
  c->cd(3);
  toy_WF->Draw("colz");
  gPad->SetLogz();
  

  TFile *f_out = new TFile("waveforms-thr--20-T50-F2000-waves0-2_sig_wf_Template.root","recreate");
  f_out->cd();
  h_av_norm_wf->Write();
  toy_WF_norm->Write();
  toy_WF->Write();
  //c->Write();
  f_out->cd();
  f_out->Close();
  delete f_out;/*
  delete h_av_norm_wf;
  delete toy_WF_norm;
  delete toy_WF;
  //delete c;
  */
}










void my_SiPM_generate_toy_waveforms(int n_events=50000){

  
  float WF_shape[1024];
  float signal_shape[1024];

  TFile *f_toy_wf = new TFile("Toy_wf.root","recreate");
  TTree *tree_toy_wf = new TTree("tree_wf","tree waveform");
  tree_toy_wf->Branch("wf",WF_shape,"wf[1024]/F");
  
  float A=20;
  float x0=200;
  float tau=150;
  float amp_tot=-10;
  float mu_noise=0;
  float sigma_noise=0.5;


  ///////functions for SiPM processes simulation
  // for simulation of signal rise range
  TF1 *fRise = new TF1("fRise","1",2,8);
  // for simulation of number of prompt CT
  TF1 *fCT = new TF1("fCT","exp(-(x-1))",1,10);
  // for simulation of secondary pulses (no matter if after pulses or delayed cross-talk)
  // possibilmente cambiarla in una gaussian molto stretta attorno ad una data posizione per essere
  TF1 *fap = new TF1("fap","1",100,110);//"exp(-x/150)",0,1024-200); 
  // for simulation of secondary pulses amplitudes (pure afterpulses)
  TF1 *fampap = new TF1("faampp","exp(-x/0.3)",0,1);
  // for simulation of pure electric noise
  TF1 *f_gauss_noise = new TF1("f_gauss_noise","(1./(sqrt(2*3.14)))*exp(-x*x/(2))",-2,2);
  // variables with generated quantities
  int n_CT=-9;
  int n_rise=-9;
  int n_rise_AP=-9;
  int x_AP=-9;
  float xap=-10;
  float amp_ap=-10;



  
  TH1D *h_CT_gen = new TH1D("h_CT_gen","h_CT_gen",10,1,11);
  TH1D *h_AP_gen = new TH1D("h_AP_gen","h_AP_gen",100,0,1024-200);

  TH2F *toy_WF = new TH2F("toy_WF","toy_WF",1024,0,1024,1000,-5,250);
  int  coeff_wf_sig;
  int coeff_wf_ap;
  int  coeff_wf_sig_rise;
  int coeff_wf_ap_rise;


  TH2F *toy_WF_deco = new TH2F("toy_WF_deco","toy_WF_deco",1024,0,1024,1000,-5,250);

  for(int i=0;i<n_events; i++){
    // get the number of optical cross talk
    n_CT=int(fCT->GetRandom());
    h_CT_gen->Fill(n_CT);
    amp_tot=A*n_CT;
    // get the position of the afterpulse
    x_AP=fap->GetRandom();
    h_AP_gen->Fill(x_AP);
    xap=x0+x_AP;
    // get the signal rise range
    n_rise=int(fRise->GetRandom());
    n_rise_AP=int(fRise->GetRandom());
    // get the ap amplitude
    amp_ap=1*n_CT;//fampap->GetRandom();
    
    for(int j=0;j<1024;j++){
      if(j<x0){ coeff_wf_sig =0;} else {coeff_wf_sig =1;}
      if(j<xap){ coeff_wf_ap =0;} else {coeff_wf_ap =1;}
      if(j<x0-n_rise||j>=x0){ coeff_wf_sig_rise =0;} else if(j>x0-n_rise&&j<x0){coeff_wf_sig_rise =1;}
      if(j<xap-n_rise_AP||j>=xap){ coeff_wf_ap_rise =0;} else if(j>xap-n_rise_AP&&j<xap){coeff_wf_ap_rise =1;}
      float wf_val=f_gauss_noise->GetRandom()+
	coeff_wf_sig_rise*(j*amp_tot/n_rise+amp_tot*(1-float(x0)/n_rise))+
	coeff_wf_sig*amp_tot*exp(-(j-x0)/tau)+
	coeff_wf_ap_rise*(j*A/n_rise_AP+A*(1-float(xap)/n_rise_AP))+
	coeff_wf_ap*amp_ap*A*exp(-(j-xap)/tau)	;
      WF_shape[j]=wf_val;
      signal_shape[j]=coeff_wf_sig*exp(-(j-x0)/tau);
      toy_WF->Fill(j,wf_val);
    }//end     for(j=0;j<1024;j++){
    

    tree_toy_wf->Fill();
    
  }// end   for(int i=0;i<n_events; i++){
  f_toy_wf->cd();
  tree_toy_wf->Write();
  tree_toy_wf->Print();
  f_toy_wf->Close();
  delete f_toy_wf;
  
  TCanvas *c_0 = new TCanvas("c_0","c_0");
  c_0->cd();
  toy_WF->Draw("colz");
  /*
    TCanvas *c = new TCanvas("c","c");
    c->Divide(2,2);
    c->cd(1);
    h_CT_gen->Draw();
    c->cd(2);
    h_AP_gen->Draw();
  */
}//end my_SiPM_toy
