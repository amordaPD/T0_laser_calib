
////////////////////////////////////////////////
////////////// Author: A. Mordà ////////////////
////////////////////////////////////////////////
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
  TH2F* xframe2_t_amp_0;
  RooPlot* xframe2_amp_0;
  RooPlot* xframe2_fit_0;
  RooPlot* xframe2_pull_0;
  TH2 *h_correlation_0;
  TH2F* xframe2_t_amp_1;
  RooPlot* xframe2_amp_1;
  RooPlot* xframe2_fit_1;
  RooPlot* xframe2_pull_1;
  TH2 *h_correlation_1;
  TH2F* xframe2_t_amp_01;
  RooPlot* xframe2_amp_01;
  RooPlot* xframe2_fit_01;
  RooPlot* xframe2_pull_01;
  TH2 *h_correlation_01;
  
};

struct Fit_results_basic{
  float T0[2];
  float T_Res_H[2];
  float T_Res_L[2];
  float frac_L_0[2];
  float delta_H_0[2];
  float delta_H_KR[2];
  float frac_H_KR[2];
  
  RooPlot* xframe2_fit_0;
  RooPlot* xframe2_fit_0_log;
  RooPlot* xframe2_pull_0;
  TH2 *h_correlation_0;
  TH1D* h_MC_nominal_lens;
  TH1D* h_MC_ring_lens;
  
};

void run_data_production(TString input_raw_name="", int my_slot=-99){

  TString file_origin="PD";
  if(my_slot>0){file_origin=="KEK";}

  /*cout<<"PRODUCING FLAT NTUPLE INPUT"<<endl;
  cscoper_DAQ("",input_raw_name,"");//produce_flat_ntuples
  */
  cout<<"PRODUCING HISTOGRAMS FOR FITTING"<<endl;
  for(int i=1;i<=4;i++){
    TString file_kind='';
    if(i==1){file_kind="recreate";}else{file_kind="update";}
    make_data_histos_column("",input_raw_name,"",file_kind,my_slot,i);
  }
  
  cout<<"PRODUCING PMT SPECTRA"<<endl;
  make_pmt_plots("",input_raw_name,file_origin);
  
}



int cscoper_DAQ(TString origin_path, TString fname, TString output_path){

  
  float _times_DAQ[32];
  float _amps_DAQ[32];
  int _chans_DAQ[32];
  
  if(origin_path==""){origin_path="D:/Padova/TOP/Time_calibration_codes/Dati_PD/1.raw_input/";}
  if(output_path==""){output_path="D:/Padova/TOP/Time_calibration_codes/Dati_PD/2.flat_input/";}
  
  float time_DAQ(-99);
  float amp_DAQ(-99);
  int channel_DAQ(-99);



  TFile *f_input = new TFile(origin_path+fname+".root");
  TTree *_tree_DAQ=(TTree*)f_input->Get("times");
  _tree_DAQ->SetBranchAddress("times",&(_times_DAQ[0]));
  _tree_DAQ->SetBranchAddress("amps",&(_amps_DAQ[0]));
  _tree_DAQ->SetBranchAddress("chans",&(_chans_DAQ[0]));


  
  TFile *f_data = new TFile(output_path+fname+"_DAQ_flat.root","recreate");
  TTree *tree_DAQ = new TTree("times","times");
  tree_DAQ->Branch("time",&time_DAQ,"time/F");
  tree_DAQ->Branch("amp",&amp_DAQ,"amp/F");
  tree_DAQ->Branch("channel",&channel_DAQ,"channel/I");
   
  Int_t _nevents = _tree_DAQ->GetEntries();
  
  Int_t _curevent;
  cout<<"Input \".root\" file opened"<<endl;


  for( _curevent =0; _curevent <_nevents; _curevent++) {
    _tree_DAQ->GetEntry(_curevent);
    for(int _ichan=0; _ichan<32;_ichan++){
      if(0<_times_DAQ[_ichan]&&_times_DAQ[_ichan]<1000){
	time_DAQ=_times_DAQ[_ichan];
	channel_DAQ=_chans_DAQ[_ichan];
	amp_DAQ=_amps_DAQ[_ichan];
	tree_DAQ->Fill();
      }
    }
    
  }
  f_data->cd();
  tree_DAQ->Write();
  f_data->Close();
  delete f_data;

  return 1;
}



void  make_data_histos_column(TString input_path, TString file_name, TString output_path, TString out_file_status, int my_slot, int my_column){
  if (input_path=="") input_path="Dati_PD/2.flat_input/";
  if (output_path=="") output_path="Dati_PD/3.column_data/";
  if(out_file_status!="recreate"&& out_file_status!="update"){
    cout<<"ERROR: incorrect kind of output file: recreate or update, output file is being recreated from scretch"<<endl;
    out_file_status="recreate";
  }
  
  //TString file_name="quartz-2pmt-waves-T77_out";
  //if(my_slot>0) file_name="run3730_f00000-f00049_allch-7nsWidth_digits_col1-4";
  TString tree_input="times"; //For PD  data
  if(my_slot>0) tree_input="laser"; //For KEK data
  
  TString inp_f =input_path+file_name+"_DAQ_flat.root";
  TFile *file_input = new TFile(inp_f);
  cout<<"input file opened : "<<inp_f<<endl;
  TTree *t_input = (TTree*)file_input->Get(tree_input);

  TString data_origin="PD_";
  if(my_slot>0) data_origin="KEK_";
  TString slotID = "";
  if(my_slot>0) slotID=Form("_slot_%i",my_slot);
  
  TString columnID =Form("col_%i",my_column);

  TString of=output_path+data_origin+file_name+slotID+"_data_histos_test.root"; //*******TO BE removed '_test'
  TFile *f_data = new TFile(of,out_file_status);
  cout<<"output file created (updating) : "<<of<<endl;

  f_data->cd();
  f_data->mkdir(Form("column_%i",my_column));
  f_data->cd(Form("column_%i",my_column));



  
  TFile *file_input_MC = new TFile("ana_laser_s01_0reso_500k.root");
  TTree *tree_MC = (TTree*)file_input_MC->Get("laser");
  TFile *file_input_MC_ring = new TFile("ana_laser_s01_0reso_ring_500k.root");
  TTree *tree_MC_ring = (TTree*)file_input_MC_ring->Get("laser");
  TH1D *h_yields = new TH1D("h_yields","event yields",8,0,9);


  
  int my_pixelID=-9;
  my_pixelID=my_column;
  int my_pixelID_data=-9;
  float upper_time;
  float lower_time;
  int KEK_PMT=13;
  int my_row=-9;
  if(my_slot<=0){upper_time=100; lower_time=0;}else{upper_time=-10; lower_time=-16;}
  int index_pmt=0;
  for(int g=1; g<=8;g++){
    if(g>4)index_pmt=1;
    cout<<"doing "<<g<<"th row"<<endl;
    my_pixelID=my_column;
    if(g<=4){my_row=g;}else{my_row=g-4;}
    my_pixelID_data=(my_column)*4+20*index_pmt-g;///////actualy the readout channel for PD testbench system data
    
    Float_t upper_bound_hist=2;
    Int_t n_bins = 200;
    TH1D *h_temp = new TH1D("h_temp","h_temp",1000,lower_time,upper_time);
    TCut cut ;
    if(my_slot<=0) cut = Form("amp>15&&0<time&&time<100&&channel==%i",my_pixelID_data); //For PD
    if(my_slot>0) cut = Form("%f<time&&time<%f&&pmt==13&&column==%i&&row==%i",lower_time,upper_time,my_column+49,my_row); //For KEK
    
    t_input->Project("h_temp","time",cut);
    h_yields->SetBinContent(g,t_input->GetEntries(cut));
    //h_temp->Draw();
    Float_t max_bin = h_temp->GetMaximumBin();
    TAxis *xaxis = h_temp->GetXaxis(); 
    Double_t max_pos = xaxis->GetBinCenter(max_bin);
    //cout<< max_pos <<endl;
    TH1D *h_time = new TH1D("h_time","Time [ns]",n_bins,-1.5,upper_bound_hist);
    t_input->Project("h_time",Form("time-%f",max_pos),cut);
    TH1D *h_amp = new TH1D("h_amp","Max Amplitude [ADC counts]",n_bins,0,300);
    t_input->Project("h_amp","amp",cut);
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
    /*    tree_MC->Project("h_MC_f1",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==1&&pixel==%i",my_pixelID));
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
    */

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
    /*    tree_MC_ring->Project("h_MC_ring_f1",Form("propTime-%f",max_pos_MC),Form("propTime<1&&fiberNo==1&&pixel==%i",my_pixelID));
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
    */
    
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
  
    
    
    f_data->mkdir(Form("column_%i/histos_%i_%i",my_column,my_column,g));
    f_data->cd(Form("column_%i/histos_%i_%i",my_column,my_column,g));
    h_amp->Write();
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
    delete h_amp;
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
  f_data->cd(Form("column_%i",my_column));
  h_yields->Write();
  f_data->cd();
  f_data->Close();
  delete f_data;
}




void make_pmt_plots(TString input_path, TString input_filebasename, TString data_origin){

  bool save = true;
  if (input_path=="") input_path="Dati_PD/3.column_data/";
  
  TCanvas *pmt_up = new TCanvas("pmt_up_"+input_filebasename,"pmt_up_"+input_filebasename);
  pmt_up->Divide(4,4);
  TCanvas *pmt_down = new TCanvas("pmt_down_"+input_filebasename,"pmt_down_"+input_filebasename);
  pmt_down->Divide(4,4);
  TCanvas *pmts = new TCanvas("pmts_"+input_filebasename,"pmts_"+input_filebasename);
  pmts->Divide(4,8);
  TCanvas *pmts_amp = new TCanvas("pmts_amp_"+input_filebasename,"pmts_amp_"+input_filebasename);
  pmts_amp->Divide(4,8);

  TH1D* h[5][9];
  TH1D* h_amp[5][9];
  TH2D *h_yield_map = new TH2D("h_yield_map","event yield maps",4,0,5,8,0,9);

  TString f_prefix = "PD_";
  if(data_origin=="KEK"){f_prefix = "KEK_";}
  TString input_filename=f_prefix+input_filebasename;
  input_filename=input_path+input_filename+"_data_histos";
  TFile *f = new TFile(input_filename+".root");
  
  for(int i=1;i<=4;i++){
    int pmt_index=1;
    int pmt_row=-1;
    TH1D *h_yields = (TH1D*)f->Get(Form("column_%i/h_yields",i));
    for(int g=1; g<=8;g++){
      if(g<=4){pmt_row=g;} else{pmt_row=g-4;}
      h[i][g]=(TH1D*)f->Get(Form("column_%i/histos_%i_%i/h_time",i,i,g));
      //***h_amp[i][g]=(TH1D*)f->Get(Form("histos_%i_%i/h_amp",i,g));
      int canvasID=i+4*(4-pmt_row);
      if(g<=4){pmt_down->cd(canvasID);
	h[i][g]->Draw("same");
	pmts->cd(16+canvasID);
	h[i][g]->Draw("same");
	pmts_amp->cd(16+canvasID);
	//***h_amp[i][g]->Draw("same");
      }else {pmt_up->cd(canvasID);
	h[i][g]->Draw("same");
	pmts->cd(canvasID);
	h[i][g]->Draw("same");
	pmts_amp->cd(canvasID);
	//***h_amp[i][g]->Draw("same");
      }
      
      h_yield_map->SetBinContent(i,g,h_yields->GetBinContent(g));
    }
    delete h_yields;
  }
  TCanvas *c_map = new TCanvas("pmts_occupancy_"+input_filebasename,"pmts_occupancy_"+input_filebasename);
  h_yield_map->DrawNormalized("colz");
  
  if(save){
    TFile *f_out = new TFile(input_path+input_filebasename+"pmt_plots.root","recreate");
    f_out->cd();
    c_map->Write();
    pmts_amp->Write();
    pmts->Write();
    pmt_up->Write();
    pmt_down->Write();
    f_out->cd();
    f_out->Close();
    delete f_out;
  }
  delete h_yield_map;
}






Fit_results_basic basic_fit_data(TString input_basepath, TString input_basefilename, int fit_model_d, int fit_model_r, int column_number, int row_number){
  if (input_basepath=="") input_basepath="Dati_PD/3.column_data/";
  
  //gROOT->ProcessLine(".L Fit_testbench.cpp"); loop_fit_PD_column(1); > ee.log
  /////////////////////////////////////////////////////////////////
  ////// HERE Starts the fit part /////////////////////////////////
  /////////////////////////////////////////////////////////////////
  Fit_results_basic Results;
  bool change_model = true;
  gROOT->ProcessLine(".x myRooPdfs/RooExpGauss.cxx+") ;
  gROOT->ProcessLine(".x myRooPdfs/RooAsymGauss.cxx+") ;

  TString input_filename="";
  TString input_column=Form("_col_%i",column_number);
  input_filename=input_basepath+input_basefilename+"_data_histos"+input_column;
  TFile *f_input = new TFile(input_filename+".root"); //FOR PD data
 
  //TFile *f_input = new TFile(Form("KEK_data_histos_slot_8_col_%i_00049_allchs.root",column_number));
  Int_t pixelID=column_number+64*(row_number-1);
  
  TH1D* h_time = f_input->Get(Form("histos_%i_%i/h_time",column_number,row_number));
 
  TH1D* h_MC_tot = f_input->Get(Form("histos_%i_%i/h_MC_tot",column_number,row_number));

  TH1D* h_MC_tot_ring = f_input->Get(Form("histos_%i_%i/h_MC_ring_tot",column_number,row_number));
 
  TCanvas *can0 = f_input->Get(Form("histos_%i_%i/c",column_number,row_number));
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
  int  bkg_Chebychev_polynomial_degree=1;//set to n to have a n+1 degree Chebychev Polynomial!!!!!!!!!
  bool add_background_component=false;
  int  amplitude_cut = -40;
  
  bool suppress_negligible_first_peak=false;
  bool do_simultaneous_fit=false;
  bool add_third_signal_pos0=false;
  //if(row_number>4){add_third_signal_pos0=true;}
  bool simulate_CB_tail=false;
  
  
  if(!print_prefit_info){MN_output_print_level_prefit=-1;}else{MN_output_print_level_prefit=MN_output_print_level;}
  bool draw_results;
  if(_draw_results=="draw"){draw_results=true;}else if(_draw_results=="blind"){draw_results=false;}else{draw_results=false;}


  
  double my_low_x=-0.7;//-0.5;//1;
  double my_up_x=0.8;//0.8;
  TAxis *xaxis_MC = h_MC_tot->GetXaxis();
  xaxis_MC->SetRange(0,(TMath::Abs(my_low_x)/(my_up_x-my_low_x))*h_MC_tot->GetNbinsX()-1);
  Float_t max_bin_MC = h_MC_tot->GetMaximumBin();
  Double_t max_first_pos_MC = xaxis_MC->GetBinCenter(max_bin_MC);
  //  cout<<"ciao "<<max_first_pos_MC<<endl;
  
  RooRealVar x("Time","Time [ns]",my_low_x,my_up_x) ;
  x.setRange("Fit_Range_L",-0.55,-0.2) ;
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
  
  low_delta_H_0=0.00;
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
   up_alpha_0=0.7;
  
  low_beta_0=0.5;
   up_beta_0=1.0;

  
   starting_delta_H_0=0.5;
   starting_delta_T_0=0.3;
   starting_sigma_L_0=0.08100;
   starting_sigma_H_0=0.08100;
   starting_sigma_T_0=0.108100;


   
   starting_alpha_0=0.125;
   starting_beta_0=0.9;
   
   starting_mean_H_0=0.0;//starting_mean_L_0+starting_delta_H_0;
   low_mean_H_0=starting_mean_H_0-0.075;//315;
   up_mean_H_0=starting_mean_H_0+0.075;//0.315;

   
   starting_mean_L_0=0.0;
   low_mean_L_0=starting_mean_L_0-0.315;
   up_mean_L_0=starting_mean_L_0+0.315;

 
   starting_alpha_CB_L=-1.5;
  starting_n_CB_L=6;

  low_alpha_CB_L=-6.5;
   up_alpha_CB_L=-0.050;
   
  low_n_CB_L=0;
   up_n_CB_L=200;


  starting_alpha_CB_H=-1.5;
  starting_n_CB_H=6;

  low_alpha_CB_H=-5;
   up_alpha_CB_H=-0.050;
   
  low_n_CB_H=0;
   up_n_CB_H=200;

   starting_alpha_CB_T=-1.5;
  starting_n_CB_T=6;

  low_alpha_CB_T=-5;
   up_alpha_CB_T=-0.050;
   
  low_n_CB_T=0;
   up_n_CB_T=20;

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
  //RooRealVar mean_T_0("mean_T_0","mean_{T}^{0}",1,0.7,1.3);
  



  RooRealVar sigma_L_0("sigma_L_0","width of L gaussian background",starting_sigma_L_0,low_sigma_L_0,up_sigma_L_0);
  /*RooRealVar sigma_H_ScF("sigma_H_ScF","width of H gaussian background",1,0.5,1.5);
    RooFormulaVar sigma_H_0("sigma_H_0", "sigma_H_0","sigma_L_0*sigma_H_ScF",RooArgList(sigma_L_0,sigma_H_ScF));*/
  RooRealVar sigma_H_0("sigma_H_0","width of H gaussian background",starting_sigma_H_0,low_sigma_H_0,up_sigma_H_0);
  RooRealVar sigma_T_0("sigma_T_0","width of T gaussian background",starting_sigma_T_0,low_sigma_T_0,up_sigma_T_0);
  
  RooRealVar sigma_L_T_0("sigma_L_T_0","width of L gaussian background",starting_sigma_L_0,low_sigma_L_0,up_sigma_L_0);
  RooRealVar sigma_H_T_0("sigma_H_T_0","width of H gaussian background",starting_sigma_H_0,low_sigma_H_0,up_sigma_H_0);
  RooRealVar sigma_T_T_0("sigma_T_T_0","width of T gaussian background",starting_sigma_T_0,low_sigma_T_0,up_sigma_T_0);

  RooRealVar  sigma_KR_T_0("sigma_KR_T_0","sigma_KR_T_0",0.08,0.01,0.1);
  RooRealVar  delta_mean_KR_T_0("delta_mean_KR_T_0","delta_mean_KR_T_0",0.100,0.05,0.3);
  RooRealVar  alpha_KR_C_0("alpha_KR_C_0","alpha_KR_C_0",0.8,0.1,1);



  
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////// FIT MODEL FOR DIRECT PATH ///////////////////////////////
  //////////////////////////////////////////////////////////////////////////////// 

  


  
  if(fit_model_d==0){
    RooGaussian PDF_L_0("PDF_L_0","gaussian L_0",x,mean_L_0,sigma_L_0) ;
  }else if(fit_model_d==1){
    //RooGaussian PDF_L_0("PDF_L_0","gaussian L_0",x,mean_L_0,sigma_L_0) ;
    RooCBShape PDF_L_0("PDF_L_0","gaussian L_0",x,mean_L_0,sigma_L_0,alpha_CB_L,n_CB_L) ;
  }else if(fit_model_d==2){
    RooExpGauss PDF_L_0("PDF_L_0","gaussian L_0",x,mean_L_0,sigma_L_0,alpha_CB_L) ;
  }else if(fit_model_d==3){
    RooAsymGauss PDF_L_0("PDF_L_0","gaussian L_0",x,mean_L_0,sigma_L_0,sigma_L_T_0) ;
  }else if(fit_model_d==4){
    RooGaussian PDF_L_C_0("PDF_H_C_0","gaussian L_0",x,mean_L_0,sigma_L_0) ;
    RooFormulaVar mean_L_T_0("mean_L_T_0","mean_L_T_0","mean_L_0+delta_mean_KR_T_0",RooArgList(mean_L_0,delta_mean_KR_T_0));
    RooExpGauss PDF_L_T_0("PDF_H_T_0","gaussian L_0",x,mean_L_T_0,sigma_KR_T_0,alpha_CB_T) ;
    RooArgList  pdfList_sig_L_0(PDF_L_C_0,PDF_L_T_0);
    RooArgList  fracList_sig_L_0(alpha_KR_C_0);
    RooAddPdf   PDF_L_0("PDF_L_0","PDF_L_0",pdfList_sig_L_0,fracList_sig_L_0,kTRUE);
  }




  
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////// FIT MODEL FOR REFLECTED PATH ///////////////////////////////
  //////////////////////////////////////////////////////////////////////////////// 

  


  
  if(fit_model_r==0){
    RooGaussian PDF_H_0("PDF_H_0","gaussian H_0",x,mean_H_0,sigma_H_0) ;//,alpha_CB_H,n_CB_H) ;
    RooGaussian PDF_T_0("PDF_T_0","gaussian T_0",x,mean_T_0,sigma_T_0) ;
  }else if(fit_model_r==1){
    RooCBShape PDF_H_0("PDF_H_0","gaussian H_0",x,mean_H_0,sigma_H_0,alpha_CB_H,n_CB_H) ;
    RooCBShape PDF_T_0("PDF_T_0","gaussian T_0",x,mean_T_0,sigma_T_0,alpha_CB_T,n_CB_T) ;
  }else if(fit_model_r==2){
    RooExpGauss PDF_H_0("PDF_H_0","gaussian H_0",x,mean_H_0,sigma_H_0,alpha_CB_H) ;
    RooExpGauss PDF_T_0("PDF_T_0","gaussian T_0",x,mean_T_0,sigma_T_0,alpha_CB_T) ;
  }else if(fit_model_r==3){
    RooAsymGauss PDF_H_0("PDF_H_0","gaussian H_0",x,mean_H_0,sigma_H_0,sigma_H_T_0) ;
    RooAsymGauss PDF_T_0("PDF_T_0","gaussian T_0",x,mean_T_0,sigma_T_0,sigma_T_T_0) ;
  }else if(fit_model_r==4){
    RooGaussian PDF_H_C_0("PDF_H_C_0","gaussian L_0",x,mean_H_0,sigma_H_0) ;
    RooFormulaVar mean_H_T_0("mean_H_T_0","mean_H_T_0","mean_H_0+delta_mean_KR_T_0",RooArgList(mean_H_0,delta_mean_KR_T_0));
    RooExpGauss PDF_H_T_0("PDF_H_T_0","gaussian L_0",x,mean_H_T_0,sigma_H_T_0,alpha_CB_T) ;
    RooArgList  pdfList_sig_H_0(PDF_H_C_0,PDF_H_T_0);
    RooArgList  fracList_sig_H_0(alpha_KR_C_0);
    RooAddPdf   PDF_H_0("PDF_H_0","PDF_H_0",pdfList_sig_H_0,fracList_sig_H_0,kTRUE);

    RooExpGauss PDF_T_0("PDF_T_0","gaussian T_0",x,mean_T_0,sigma_T_0,alpha_CB_T) ;

    /*
      RooGaussian PDF_T_C_0("PDF_T_C_0","gaussian L_0",x,mean_T_0,sigma_T_0) ;
      RooRealVar  sigma_T_T_0("sigma_T_T_0","sigma_T_T_0",0.2,0.15,0.5);
      RooRealVar  delta_mean_T_T_0("delta_mean_T_T_0","delta_mean_T_T_0",0.02,0.0,1);
      RooFormulaVar mean_T_T_0("mean_T_T_0","mean_T_T_0","mean_T_0+delta_mean_T_T_0",RooArgList(mean_T_0,delta_mean_T_T_0));
      RooGaussian PDF_T_T_0("PDF_T_T_0","gaussian L_0",x,mean_T_T_0,sigma_T_T_0) ;
      RooRealVar  alpha_T_C_0("alpha_T_C_0","alpha_T_C_0",0.9,0.7,1);
      RooArgList  pdfList_sig_T_0(PDF_T_C_0,PDF_T_T_0);
      RooArgList  fracList_sig_T_0(alpha_KR_C_0);
      RooAddPdf   PDF_T_0("PDF_T_0","PDF_T_0",pdfList_sig_T_0,fracList_sig_T_0,kTRUE);
    */
    
  }



  
  ////////////////////////////////////////////////////////////////////
  //////////// BUILDING FULL MODEL /////////////////////////////
  //////////////////////////////////////////////////////////////






  
  RooRealVar alpha_0("alpha_0","alpha_0",starting_alpha_0,low_alpha_0,up_alpha_0);
  //  alpha_0.setConstant(kTRUE);
  RooRealVar beta_0("beta_0","beta_0",starting_beta_0,low_beta_0,up_beta_0);
  RooFormulaVar Frac_L_0("Frac_L_0","Frac_L_0","alpha_0",RooArgList(alpha_0));
  RooFormulaVar Frac_H_0("Frac_H_0","Frac_H_0","beta_0-alpha_0*beta_0",RooArgList(beta_0,alpha_0));
  RooFormulaVar Frac_T_0("Frac_T_0","Frac_T_0","1-alpha_0-beta_0+alpha_0*beta_0",RooArgList(alpha_0,beta_0));
  RooArgList  pdfList_sig_0(PDF_L_0,PDF_H_0); if(add_third_signal_pos0) {pdfList_sig_0.add(PDF_T_0);}//
  RooArgList  fracList_sig_0(alpha_0); if(add_third_signal_pos0) {fracList_sig_0.add(beta_0);}
  
  if(change_model&&row_number<=3){
    RooArgList  pdfList_sig_0(PDF_H_0); if(add_third_signal_pos0) {pdfList_sig_0.add(PDF_T_0);}//
    RooArgList  fracList_sig_0(); if(add_third_signal_pos0) {fracList_sig_0.add(beta_0);}
    RooAddPdf PDF_sig_0("PDF_sig_0","PDF_sig_0",pdfList_sig_0,fracList_sig_0,kTRUE);
  }else{
    RooArgList  pdfList_sig_0(PDF_L_0,PDF_H_0); if(add_third_signal_pos0) {pdfList_sig_0.add(PDF_T_0);}//
    RooArgList  fracList_sig_0(alpha_0); if(add_third_signal_pos0) {fracList_sig_0.add(beta_0);}
    RooAddPdf PDF_sig_0("PDF_sig_0","PDF_sig_0",pdfList_sig_0,fracList_sig_0,kTRUE);
  }


  
  RooRealVar a0_0("a0_0", "", 0.0, -10, 10);
  RooRealVar a1_0("a1_0", "", 0.0, -20, 20);
  RooRealVar a2_0("a2_0", "", 0.00, -20, 20);
  RooArgList  coeffList_sig_0(a0_0);
  if(bkg_Chebychev_polynomial_degree>=1){
    coeffList_sig_0.add(a1_0);    
    if(bkg_Chebychev_polynomial_degree>=2){
      coeffList_sig_0.add(a2_0);
    }
  }
  
  alpha_CB_L.setVal(-0.4);
  alpha_CB_H.setVal(-0.4);
  //alpha_CB_L.setConstant(kTRUE);
  //alpha_CB_H.setConstant(kTRUE);
  
  
  RooChebychev PDF_B_0("PDF_B_0","PDF_B_0",x,coeffList_sig_0);
  //RooGaussian PDF_B_0("PDF_B_0","gaussian T_0",x,mean_T_0,sigma_T_0) ;
  RooRealVar  Frac_sig_0("Frac_sig_0","fraction of sig events", 0.9, 0.7,1.0);
  if(!add_background_component){
    a0_0.setConstant(kTRUE);
    a1_0.setConstant(kTRUE);
    a2_0.setConstant(kTRUE);
    Frac_sig_0.setVal(1);
    Frac_sig_0.setConstant(kTRUE);
  }
  RooArgList  pdfList_0(PDF_sig_0,PDF_B_0);
  RooArgList  fracList_0(Frac_sig_0);
  RooAddPdf   model_0("model_0","model_0",pdfList_0,fracList_0,kTRUE);
  RooAddPdf   model_0_b("model_0_b","model_0_b",pdfList_0,fracList_0,kTRUE);

  /*
  if(change_model&&row_number<=3){
    cout<<"the two peaks are too close (Delta_H<0.15): repeating the fit with only one signal "<<endl; 
    Delta_H_0.setVal(0.00);
    //mean_L_0.setVal(0.00);
    //sigma_L_0.setVal(0.00);
    alpha_0.setVal(0.00);
    alpha_0.setConstant(kTRUE);
    //mean_L_0.setConstant(kTRUE);
    //sigma_L_0.setConstant(kTRUE);
    //Delta_H_0.setConstant(kTRUE);
  }
  */
  if(do_prefit){

    /*
  mean_H_0.setConstant(kTRUE);
  PDF_L_0.fitTo(ds_0_H,Save(),Minimizer(Type_minim_pf,Algo_minim_pf),Strategy(2),SumW2Error(kFALSE),PrintLevel(MN_output_print_level_prefit),PrintEvalErrors(-1),Range("Fit_Range_L"),Verbose(kFALSE));//);
  Delta_H_0.setConstant(kTRUE);
  sigma_L_0.setConstant(kTRUE);
  alpha_CB_L.setConstant(kTRUE);
  mean_H_0.setConstant(kFALSE);
    */
    


    
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
 
  if(change_model&&row_number<=3){
    cout<<"the two peaks are too close (Delta_H<0.15): repeating the fit with only one signal "<<endl; 
    Delta_H_0.setVal(0.00);
    alpha_0.setVal(0.00);
    sigma_L_0.setVal(0.100);
    Delta_H_0.setConstant(kTRUE);
    sigma_L_0.setConstant(kTRUE);
    alpha_CB_L.setConstant(kTRUE);
    n_CB_L.setConstant(kTRUE);
    alpha_0.setConstant(kTRUE);
    
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
    
    
  }
  
  
  fit_results_0->Print("v");
  RooPlot* xframe2_0 = x.frame(Title("Fit")) ;
  ds_0.plotOn(xframe2_0);//,DataError(RooAbsData::SumW2)) ;
  TH2 *h_correlation_0 = fit_results_0->correlationHist(); h_correlation_0->SetTitle("correlation matrix ");h_correlation_0->GetYaxis()->SetLabelSize(0.1); h_correlation_0->GetYaxis()->SetLabelFont(70);h_correlation_0->GetXaxis()->SetLabelSize(0.075); h_correlation_0->GetXaxis()->SetLabelFont(70); h_correlation_0->SetMarkerSize(2);
  model_0.plotOn(xframe2_0,Range("Fit_Range_0"),Components("PDF_L_C_0"),LineStyle(2),LineColor(1)) ;
  model_0.plotOn(xframe2_0,Range("Fit_Range_0"),Components("PDF_L_T_0"),LineStyle(2),LineColor(1)) ;
  model_0.plotOn(xframe2_0,Range("Fit_Range_0"),Components("PDF_H_C_0"),LineStyle(3),LineColor(8)) ;
  model_0.plotOn(xframe2_0,Range("Fit_Range_0"),Components("PDF_H_T_0"),LineStyle(3),LineColor(8)) ;
  model_0.plotOn(xframe2_0,Range("Fit_Range_0"),Components("PDF_L_0"),LineStyle(1),LineColor(1)) ;
  model_0.plotOn(xframe2_0,Range("Fit_Range_0"),Components("PDF_H_0"),LineStyle(1),LineColor(8)) ;
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
  model_0.plotOn(xframe2_0_log,Range("Fit_Range_0"),Components("PDF_L_C_0"),LineStyle(kDashed),LineColor(1)) ;
  model_0.plotOn(xframe2_0_log,Range("Fit_Range_0"),Components("PDF_L_T_0"),LineStyle(kDashed),LineColor(1)) ;
  model_0.plotOn(xframe2_0_log,Range("Fit_Range_0"),Components("PDF_H_C_0"),LineStyle(kDashed),LineColor(8)) ;
  model_0.plotOn(xframe2_0_log,Range("Fit_Range_0"),Components("PDF_H_T_0"),LineStyle(kDashed),LineColor(8)) ;
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
    
   /*
    TCanvas* c_Fit_0 = new TCanvas(Form("c_Fit_0_col%i_row%i",column_number,row_number),Form("c_Fit_0_col%i_row%i",column_number,row_number));
    xframe2_0->Draw() ;
    TFile *f_out = new TFile("file_output.root","update");
    f_out->cd();
    c_Fit_0->Write();
    f_out->Close();
    delete f_out;
    
    */
  
  }


  Results.xframe2_fit_0=xframe2_0;
  Results.xframe2_fit_0_log=xframe2_0_log;
  Results.xframe2_pull_0=xframe4_0;
  Results.h_correlation_0=h_correlation_0;
  Results.T_Res_H[0]=sigma_H_0.getVal();
  Results.T_Res_H[1]=sigma_H_0.getError();
  Results.T0[0]=mean_H_0.getVal();
  Results.T0[1]=mean_H_0.getError();
  if(!change_model||row_number>3){
    Results.T_Res_L[0]=sigma_L_0.getVal();
    Results.T_Res_L[1]=sigma_L_0.getError();
    Results.delta_H_0[0]=Delta_H_0.getVal();
    Results.delta_H_0[1]=Delta_H_0.getError();
    Results.frac_L_0[0]=alpha_0.getVal();
    Results.frac_L_0[1]=alpha_0.getError();
  }else{
    Results.T_Res_L[0]=sigma_H_0.getVal();
    Results.T_Res_L[1]=sigma_H_0.getError();
    Results.delta_H_0[0]=0.0;
    Results.delta_H_0[1]=0.0;
    Results.frac_L_0[0]=0.0;
    Results.frac_L_0[1]=0.0;
  }

  if(fit_model_r==4){
    Results.delta_H_KR[0]=delta_mean_KR_T_0.getVal();
    Results.delta_H_KR[1]=delta_mean_KR_T_0.getError();
    Results.frac_H_KR[0]=alpha_KR_C_0.getVal();
    Results.frac_H_KR[1]=alpha_KR_C_0.getError();
  }else{
    Results.delta_H_KR[0]=0.0;
    Results.delta_H_KR[1]=0.0;
    Results.frac_H_KR[0]=0.0;
    Results.frac_H_KR[1]=0.0;
  }
    

  return Results;






  
}


void loop_fit_column(TString input_basepath, TString input_basefilename, TString output_basepath, int column_number,int fit_model_ID){
  /*
To Run:
gROOT->ProcessLine(".L Fit_testbench.cpp"); loop_fit_column("","file","",3,4)
   */


  if (input_basepath=="") input_basepath="Dati_PD/3.column_data/";
  if (output_basepath=="") output_basepath="Dati_PD/4.fit_results/";
  TCanvas *cc = new TCanvas("cc","cc",0,0,1124,500);
  cc->Divide(8,3);
  float x[8];
  float err_x[8];
  float T_0[8];
  float err_T_0[8];
  float sigma_H_0[8];
  float err_sigma_H_0[8];
  float sigma_L_0[8];
  float err_sigma_L_0[8];
  float frac_L[8];
  float err_frac_L[8];
  float delta_H[8];
  float err_delta_H[8];
  float frac_H_KR[8];
  float err_frac_H_KR[8];
  float delta_H_KR[8];
  float err_delta_H_KR[8];
  
  TF1 *line_0s = new TF1("line_0s","0",-100,100);line_0s->SetLineColor(8);
  TF1 *line_1ps = new TF1("line_1ps","1",-100,100);line_1ps->SetLineColor(4);
  TF1 *line_1ns = new TF1("line_1ns","-1",-100,100);line_1ns->SetLineColor(4);
  TF1 *line_2ps = new TF1("line_2ps","2",-100,100);line_2ps->SetLineColor(kOrange-2);
  TF1 *line_2ns = new TF1("line_2ns","-2",-100,100);line_2ns->SetLineColor(kOrange-2);
  TF1 *line_3ps = new TF1("line_3ps","3",-100,100);line_3ps->SetLineColor(2);
  TF1 *line_3ns = new TF1("line_3ns","-3",-100,100);line_3ns->SetLineColor(2);

  TString output_filename="";
  TString out_column=Form("_col_%i.root",column_number);
  output_filename=input_basefilename+"_fit_results"+out_column;
  TFile *f_result = new TFile(output_basepath+output_filename,"recreate");
  for(int h=1;h<=8;h++){
    
    x[h-1]=h; err_x[h-1]=0;
    Fit_results_basic my_Results;                        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    my_Results = basic_fit_data("",input_basefilename,0,fit_model_ID,column_number,h);   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cc->cd(h);
    my_Results.xframe2_fit_0->Draw() ;
    cc->cd(8+h);gPad->SetLeftMargin(0.15) ; my_Results.xframe2_pull_0->GetYaxis()->SetTitleOffset(1.6) ; gPad->SetGridy(); my_Results.xframe2_pull_0->GetYaxis()->SetRangeUser(-5,5); my_Results.xframe2_pull_0->Draw() ;
    line_3ns->Draw("same");line_3ps->Draw("same");
    line_2ns->Draw("same");line_2ps->Draw("same");
    line_1ns->Draw("same");line_1ps->Draw("same");
    line_0s->Draw("same");line_0s->Draw("same");
    cc->cd(16+h);
    my_Results.h_correlation_0->Draw("colz:text");
    T_0[h-1]=my_Results.T0[0];
    err_T_0[h-1]=my_Results.T0[1];
    sigma_H_0[h-1]=my_Results.T_Res_H[0];
    sigma_L_0[h-1]=my_Results.T_Res_L[0];
    err_sigma_H_0[h-1]=my_Results.T_Res_H[1];
    err_sigma_L_0[h-1]=my_Results.T_Res_L[1];
    frac_L[h-1]=my_Results.frac_L_0[0];
    err_frac_L[h-1]=my_Results.frac_L_0[1];
    delta_H[h-1]=my_Results.delta_H_0[0];
    err_delta_H[h-1]=my_Results.delta_H_0[1];
    frac_H_KR[h-1]=my_Results.frac_H_KR[0];
    err_frac_H_KR[h-1]=my_Results.frac_H_KR[1];
    delta_H_KR[h-1]=my_Results.delta_H_KR[0];
    err_delta_H_KR[h-1]=my_Results.delta_H_KR[1];
    
    
 
  }


    TF1 *line = new TF1("line","0.100",-100,100);
    line->SetTitle("100 ps");



    TGraphErrors *SIGMA_H_0 = new TGraphErrors(8,x,sigma_H_0,err_x,err_sigma_H_0);SIGMA_H_0->SetTitle("#delta t_{H}");SIGMA_H_0->SetName("tge_Res_T_H");
    TGraphErrors *SIGMA_L_0 = new TGraphErrors(8,x,sigma_L_0,err_x,err_sigma_L_0);SIGMA_L_0->SetTitle("#delta t_{L}");SIGMA_L_0->SetName("tge_Res_L_H");

    SIGMA_L_0->SetMarkerStyle(20);
    SIGMA_H_0->SetMarkerStyle(20);
    SIGMA_L_0->SetMarkerSize(2);
    SIGMA_H_0->SetMarkerSize(2);
    SIGMA_L_0->SetMarkerColor(2);
    SIGMA_H_0->SetMarkerColor(4);
        
    TMultiGraph *mg_SIGMA_0 = new TMultiGraph();
    mg_SIGMA_0->SetTitle("Time-resolution (#delta t) vs Channel");
    mg_SIGMA_0->Add(SIGMA_L_0);
    mg_SIGMA_0->Add(SIGMA_H_0);
    
    TCanvas* c_Res_pos0 = new TCanvas("c_Res_pos0","c_Res_pos0",0,0,1124,700);
    mg_SIGMA_0->Draw("AP");
    line->Draw("same");
    gPad->Update();
    mg_SIGMA_0->GetXaxis()->SetTitle("row");
    mg_SIGMA_0->GetYaxis()->SetTitle("#delta t [ns]");
    gPad->Update();
    gPad->BuildLegend();
    
    TGraphErrors *tge_T_0 = new TGraphErrors(8,x,T_0,err_x,err_T_0); tge_T_0->SetTitle("T_{0}");tge_T_0->SetName("tge_T0");

    tge_T_0->SetMarkerStyle(20);
    tge_T_0->SetMarkerSize(2);
    tge_T_0->SetMarkerColor(2);

    TCanvas* c_T0 = new TCanvas("c_T0","c_T0",0,0,1124,700);
    tge_T_0->Draw("AP");
    gPad->Update();
    tge_T_0->GetXaxis()->SetTitle("row");
    tge_T_0->GetYaxis()->SetTitle("T_{0} [ns]");
    gPad->Update();
    gPad->BuildLegend();
  
    TGraphErrors *FRAC_L_0 = new TGraphErrors(8,x,frac_L,err_x,err_frac_L); FRAC_L_0->SetTitle("f_{L}");FRAC_L_0->SetName("tge_f_L");

    FRAC_L_0->SetMarkerStyle(20);
    FRAC_L_0->SetMarkerSize(2);
    FRAC_L_0->SetMarkerColor(2);
  
        
    TMultiGraph *mg_FRAC_0 = new TMultiGraph();
    mg_FRAC_0->SetTitle("fraction of low time peak (f_{L}) vs Channel");
    mg_FRAC_0->Add(FRAC_L_0);
   



    TGraphErrors *DELTA_H_0 = new TGraphErrors(8,x,delta_H,err_x,err_delta_H); DELTA_H_0->SetTitle("#Delta_{H}");DELTA_H_0->SetName("tge_Delta_H");

    DELTA_H_0->SetMarkerStyle(20);
    DELTA_H_0->SetMarkerSize(2);
    DELTA_H_0->SetMarkerColor(4);
  
    TMultiGraph *mg_DELTA_0 = new TMultiGraph();
    mg_DELTA_0->SetTitle("Time separation between main peaks (#Delta_{H}) vs Channel");
    mg_DELTA_0->Add(DELTA_H_0);



    TGraphErrors *DELTA_H_KR_0 = new TGraphErrors(8,x,delta_H_KR,err_x,err_delta_H_KR); DELTA_H_KR_0->SetTitle("#Delta_{H}^{KR}");DELTA_H_KR_0->SetName("tge_Delta_KR");

    DELTA_H_KR_0->SetMarkerStyle(20);
    DELTA_H_KR_0->SetMarkerSize(2);
    DELTA_H_KR_0->SetMarkerColor(4);
  
    TMultiGraph *mg_DELTA_H_KR_0 = new TMultiGraph();
    mg_DELTA_H_KR_0->SetTitle("Time separation between prompt and kinematic recoil  vs Channel");
    mg_DELTA_H_KR_0->Add(DELTA_H_KR_0);
    
    TGraphErrors *FRAC_H_KR_0 = new TGraphErrors(8,x,frac_H_KR,err_x,err_frac_H_KR); FRAC_H_KR_0->SetTitle("f_{KR}^{core}");FRAC_H_KR_0->SetName("tge_f_KR");

    FRAC_H_KR_0->SetMarkerStyle(20);
    FRAC_H_KR_0->SetMarkerSize(2);
    FRAC_H_KR_0->SetMarkerColor(4);
  
    TMultiGraph *mg_FRAC_H_KR_0 = new TMultiGraph();
    mg_FRAC_H_KR_0->SetTitle("Time separation between prompt and kinematic recoil  vs Channel");
    mg_FRAC_H_KR_0->Add(FRAC_H_KR_0);


    
    TCanvas* c_frac_pos0 = new TCanvas("c_frac_pos0","c_frac_pos0",0,0,1124,700);
    mg_FRAC_0->Draw("AP");
    gPad->Update();
    mg_FRAC_0->GetXaxis()->SetTitle("row");
    mg_FRAC_0->GetYaxis()->SetTitle("f_{L}");
    gPad->Update();
    gPad->BuildLegend();



    TCanvas* c_delta_pos0 = new TCanvas("c_delta_pos0","c_delta_pos0",0,0,1124,700);
    mg_DELTA_0->Draw("AP");
    gPad->Update();
    mg_DELTA_0->GetXaxis()->SetTitle("row");
    mg_DELTA_0->GetYaxis()->SetTitle("#Delta_{H}");
    gPad->Update();
    gPad->BuildLegend();
    
    TCanvas* c_delta_KR = new TCanvas("c_delta_KR","c_delta_KR",0,0,1124,700);
    mg_DELTA_H_KR_0->Draw("AP");
    gPad->Update();
    mg_DELTA_H_KR_0->GetXaxis()->SetTitle("row");
    mg_DELTA_H_KR_0->GetYaxis()->SetTitle("#Delta_{H}^{KR}");
    gPad->Update();
    gPad->BuildLegend();
    TCanvas* c_frac_KR = new TCanvas("c_frac_KR","c_frac_KR",0,0,1124,700);
    mg_FRAC_H_KR_0->Draw("AP");
    gPad->Update();
    mg_FRAC_H_KR_0->GetXaxis()->SetTitle("row");
    mg_FRAC_H_KR_0->GetYaxis()->SetTitle("f_{H}^{KR}");
    gPad->Update();
    gPad->BuildLegend();
    
    f_result->cd();
    cc->Write();
    c_frac_pos0->Write();
    c_delta_pos0->Write();
    c_Res_pos0->Write();
    c_frac_KR->Write();
    c_delta_KR->Write();
    mg_SIGMA_0->Write();
    mg_FRAC_0->Write();
    mg_FRAC_H_KR_0->Write();
    mg_DELTA_H_KR_0->Write();
    tge_T_0->Write();
    SIGMA_H_0->Write();
    SIGMA_L_0->Write();
    FRAC_L_0->Write();
    DELTA_H_KR_0->Write();
    FRAC_H_KR_0->Write();
    
    f_result->Close();
    delete f_result;
}









void loop_fit_row(TString input_basepath, TString input_basefilename, TString output_basepath, int row_number,int fit_model_ID){
  /*
To Run:
gROOT->ProcessLine(".L Fit_testbench.cpp"); loop_fit_column("","file","",3,4)
   */

  
  if (input_basepath=="") input_basepath="Dati_PD/3.column_data/";
  if (output_basepath=="") output_basepath="Dati_PD/4.fit_results/";
  
  TCanvas *cc = new TCanvas("cc","cc",0,0,1124,500);
  cc->Divide(4,3);
  float x[4];
  float err_x[4];
  float T_0[4];
  float err_T_0[4];
  float sigma_H_0[4];
  float err_sigma_H_0[4];
  float sigma_L_0[4];
  float err_sigma_L_0[4];
  float frac_L[4];
  float err_frac_L[4];
  float delta_H[4];
  float err_delta_H[4];
  float frac_H_KR[4];
  float err_frac_H_KR[4];
  float delta_H_KR[4];
  float err_delta_H_KR[4];
  
  TF1 *line_0s = new TF1("line_0s","0",-100,100);line_0s->SetLineColor(8);
  TF1 *line_1ps = new TF1("line_1ps","1",-100,100);line_1ps->SetLineColor(4);
  TF1 *line_1ns = new TF1("line_1ns","-1",-100,100);line_1ns->SetLineColor(4);
  TF1 *line_2ps = new TF1("line_2ps","2",-100,100);line_2ps->SetLineColor(kOrange-2);
  TF1 *line_2ns = new TF1("line_2ns","-2",-100,100);line_2ns->SetLineColor(kOrange-2);
  TF1 *line_3ps = new TF1("line_3ps","3",-100,100);line_3ps->SetLineColor(2);
  TF1 *line_3ns = new TF1("line_3ns","-3",-100,100);line_3ns->SetLineColor(2);

  TString output_filename="";
  TString out_row=Form("_row_%i.root",row_number);
  TString out_fitID =Form("_modelID_%i",fit_model_ID);
  output_filename=input_basefilename+"_fit_results"+out_fitID+out_row;
  TFile *f_result = new TFile(output_basepath+output_filename,"recreate");
  for(int h=1;h<=4;h++){
    x[h-1]=h; err_x[h-1]=0;
    Fit_results_basic my_Results;                        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    my_Results = basic_fit_data("",input_basefilename,fit_model_ID,h,row_number);   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cc->cd(h);
    my_Results.xframe2_fit_0->Draw() ;
    cc->cd(4+h);gPad->SetLeftMargin(0.15) ; my_Results.xframe2_pull_0->GetYaxis()->SetTitleOffset(1.6) ; gPad->SetGridy(); my_Results.xframe2_pull_0->GetYaxis()->SetRangeUser(-5,5); my_Results.xframe2_pull_0->Draw() ;
    line_3ns->Draw("same");line_3ps->Draw("same");
    line_2ns->Draw("same");line_2ps->Draw("same");
    line_1ns->Draw("same");line_1ps->Draw("same");
    line_0s->Draw("same");line_0s->Draw("same");
    cc->cd(8+h);
    my_Results.h_correlation_0->Draw("colz:text");
    T_0[h-1]=my_Results.T0[0];
    err_T_0[h-1]=my_Results.T0[1];
    sigma_H_0[h-1]=my_Results.T_Res_H[0];
    sigma_L_0[h-1]=my_Results.T_Res_L[0];
    err_sigma_H_0[h-1]=my_Results.T_Res_H[1];
    err_sigma_L_0[h-1]=my_Results.T_Res_L[1];
    frac_L[h-1]=my_Results.frac_L_0[0];
    err_frac_L[h-1]=my_Results.frac_L_0[1];
    delta_H[h-1]=my_Results.delta_H_0[0];
    err_delta_H[h-1]=my_Results.delta_H_0[1];
    frac_H_KR[h-1]=my_Results.frac_H_KR[0];
    err_frac_H_KR[h-1]=my_Results.frac_H_KR[1];
    delta_H_KR[h-1]=my_Results.delta_H_KR[0];
    err_delta_H_KR[h-1]=my_Results.delta_H_KR[1];
    
    
 
  }


    TF1 *line = new TF1("line","0.100",-100,100);
    line->SetTitle("100 ps");

    TGraphErrors *SIGMA_H_0 = new TGraphErrors(4,x,sigma_H_0,err_x,err_sigma_H_0);SIGMA_H_0->SetTitle("#delta t_{H}");SIGMA_H_0->SetName("tge_Res_T_H");
    TGraphErrors *SIGMA_L_0 = new TGraphErrors(4,x,sigma_L_0,err_x,err_sigma_L_0);SIGMA_L_0->SetTitle("#delta t_{L}");SIGMA_L_0->SetName("tge_Res_L_H");

    SIGMA_L_0->SetMarkerStyle(20);
    SIGMA_H_0->SetMarkerStyle(20);
    SIGMA_L_0->SetMarkerSize(2);
    SIGMA_H_0->SetMarkerSize(2);
    SIGMA_L_0->SetMarkerColor(2);
    SIGMA_H_0->SetMarkerColor(4);
        
    TMultiGraph *mg_SIGMA_0 = new TMultiGraph();
    mg_SIGMA_0->SetTitle("Time-resolution (#delta t) vs Channel");
    mg_SIGMA_0->Add(SIGMA_L_0);
    mg_SIGMA_0->Add(SIGMA_H_0);
    
    TCanvas* c_Res_pos0 = new TCanvas("c_Res_pos0","c_Res_pos0",0,0,1124,700);
    mg_SIGMA_0->Draw("AP");
    line->Draw("same");
    gPad->Update();
    mg_SIGMA_0->GetXaxis()->SetTitle("row");
    mg_SIGMA_0->GetYaxis()->SetTitle("#delta t [ns]");
    gPad->Update();
    gPad->BuildLegend();
    
    TGraphErrors *tge_T_0 = new TGraphErrors(4,x,T_0,err_x,err_T_0); tge_T_0->SetTitle("T_{0}");tge_T_0->SetName("tge_T0");

    tge_T_0->SetMarkerStyle(20);
    tge_T_0->SetMarkerSize(2);
    tge_T_0->SetMarkerColor(2);

    TCanvas* c_T0 = new TCanvas("c_T0","c_T0",0,0,1124,700);
    tge_T_0->Draw("AP");
    gPad->Update();
    tge_T_0->GetXaxis()->SetTitle("row");
    tge_T_0->GetYaxis()->SetTitle("T_{0} [ns]");
    gPad->Update();
    gPad->BuildLegend();
  
    TGraphErrors *FRAC_L_0 = new TGraphErrors(4,x,frac_L,err_x,err_frac_L); FRAC_L_0->SetTitle("f_{L}");FRAC_L_0->SetName("tge_f_L");

    FRAC_L_0->SetMarkerStyle(20);
    FRAC_L_0->SetMarkerSize(2);
    FRAC_L_0->SetMarkerColor(2);
  
        
    TMultiGraph *mg_FRAC_0 = new TMultiGraph();
    mg_FRAC_0->SetTitle("fraction of low time peak (f_{L}) vs Channel");
    mg_FRAC_0->Add(FRAC_L_0);
   



    TGraphErrors *DELTA_H_0 = new TGraphErrors(4,x,delta_H,err_x,err_delta_H); DELTA_H_0->SetTitle("#Delta_{H}");DELTA_H_0->SetName("tge_Delta_H");

    DELTA_H_0->SetMarkerStyle(20);
    DELTA_H_0->SetMarkerSize(2);
    DELTA_H_0->SetMarkerColor(4);
  
    TMultiGraph *mg_DELTA_0 = new TMultiGraph();
    mg_DELTA_0->SetTitle("Time separation between main peaks (#Delta_{H}) vs Channel");
    mg_DELTA_0->Add(DELTA_H_0);



    TGraphErrors *DELTA_H_KR_0 = new TGraphErrors(4,x,delta_H_KR,err_x,err_delta_H_KR); DELTA_H_KR_0->SetTitle("#Delta_{H}^{KR}");DELTA_H_KR_0->SetName("tge_Delta_KR");

    DELTA_H_KR_0->SetMarkerStyle(20);
    DELTA_H_KR_0->SetMarkerSize(2);
    DELTA_H_KR_0->SetMarkerColor(4);
  
    TMultiGraph *mg_DELTA_H_KR_0 = new TMultiGraph();
    mg_DELTA_H_KR_0->SetTitle("Time separation between prompt and kinematic recoil  vs Channel");
    mg_DELTA_H_KR_0->Add(DELTA_H_KR_0);
    
    TGraphErrors *FRAC_H_KR_0 = new TGraphErrors(4,x,frac_H_KR,err_x,err_frac_H_KR); FRAC_H_KR_0->SetTitle("f_{KR}^{core}");FRAC_H_KR_0->SetName("tge_f_KR");

    FRAC_H_KR_0->SetMarkerStyle(20);
    FRAC_H_KR_0->SetMarkerSize(2);
    FRAC_H_KR_0->SetMarkerColor(4);
  
    TMultiGraph *mg_FRAC_H_KR_0 = new TMultiGraph();
    mg_FRAC_H_KR_0->SetTitle("Time separation between prompt and kinematic recoil  vs Channel");
    mg_FRAC_H_KR_0->Add(FRAC_H_KR_0);



    
    TCanvas* c_frac_pos0 = new TCanvas("c_frac_pos0","c_frac_pos0",0,0,1124,700);
    mg_FRAC_0->Draw("AP");
    gPad->Update();
    mg_FRAC_0->GetXaxis()->SetTitle("row");
    mg_FRAC_0->GetYaxis()->SetTitle("f_{L}");
    gPad->Update();
    gPad->BuildLegend();



    TCanvas* c_delta_pos0 = new TCanvas("c_delta_pos0","c_delta_pos0",0,0,1124,700);
    mg_DELTA_0->Draw("AP");
    gPad->Update();
    mg_DELTA_0->GetXaxis()->SetTitle("row");
    mg_DELTA_0->GetYaxis()->SetTitle("#Delta_{H}");
    gPad->Update();
    gPad->BuildLegend();
    
    TCanvas* c_delta_KR = new TCanvas("c_delta_KR","c_delta_KR",0,0,1124,700);
    mg_DELTA_H_KR_0->Draw("AP");
    gPad->Update();
    mg_DELTA_H_KR_0->GetXaxis()->SetTitle("row");
    mg_DELTA_H_KR_0->GetYaxis()->SetTitle("#Delta_{H}^{KR}");
    gPad->Update();
    gPad->BuildLegend();
    TCanvas* c_frac_KR = new TCanvas("c_frac_KR","c_frac_KR",0,0,1124,700);
    mg_FRAC_H_KR_0->Draw("AP");
    gPad->Update();
    mg_FRAC_H_KR_0->GetXaxis()->SetTitle("row");
    mg_FRAC_H_KR_0->GetYaxis()->SetTitle("f_{H}^{KR}");
    gPad->Update();
    gPad->BuildLegend();
    
    f_result->cd();
    cc->Write();
    c_frac_pos0->Write();
    c_delta_pos0->Write();
    c_Res_pos0->Write();
    c_frac_KR->Write();
    c_delta_KR->Write();
    mg_SIGMA_0->Write();
    mg_FRAC_0->Write();
    mg_FRAC_H_KR_0->Write();
    mg_DELTA_H_KR_0->Write();
    tge_T_0->Write();
    SIGMA_H_0->Write();
    SIGMA_L_0->Write();
    FRAC_L_0->Write();
    DELTA_H_KR_0->Write();
    FRAC_H_KR_0->Write();
    
    f_result->Close();
    delete f_result;
}


















/*
TGraphErrors *tge_res_1 = (TGraphErrors*)_file0->Get("tge_Res_T_H")
TGraphErrors *tge_res_2 = (TGraphErrors*)_file1->Get("tge_Res_T_H")
TGraphErrors *tge_res_3 = (TGraphErrors*)_file2->Get("tge_Res_T_H")
TGraphErrors *tge_res_4 = (TGraphErrors*)_file3->Get("tge_Res_T_H")
tge_res_1->SetLineColor(1)
tge_res_2->SetLineColor(2)
tge_res_3->SetLineColor(4)
tge_res_4->SetLineColor(8)
tge_res_1->SetMarkerColor(1)
tge_res_2->SetMarkerColor(2)
tge_res_3->SetMarkerColor(4)
tge_res_4->SetMarkerColor(8)
tge_res_1->SetTitle("#sigma ")
tge_res_2->SetTitle("#sigma ")
tge_res_3->SetTitle("#sigma ")
tge_res_4->SetTitle("#sigma ")
TMultiGraph *tmg = new TMultiGraph()
tmg->Add(tge_res_1)
tmg->Add(tge_res_2)
tmg->Add(tge_res_3)
tmg->Add(tge_res_4)
tmg->Draw("AP")
gPad->BuildLegend()
gPad->Update()
TGraphErrors *tge_T0_1 = (TGraphErrors*)_file0->Get("tge_T0")
TGraphErrors *tge_T0_2 = (TGraphErrors*)_file1->Get("tge_T0")
TGraphErrors *tge_T0_3 = (TGraphErrors*)_file2->Get("tge_T0")
TGraphErrors *tge_T0_4 = (TGraphErrors*)_file3->Get("tge_T0")
tge_T0_1->SetLineColor(1)
tge_T0_2->SetLineColor(2)
tge_T0_3->SetLineColor(4)
tge_T0_4->SetLineColor(8)
tge_T0_1->SetMarkerColor(1)
tge_T0_2->SetMarkerColor(2)
tge_T0_3->SetMarkerColor(4)
tge_T0_4->SetMarkerColor(8)
tge_T0_1->SetTitle("T_0 ")
tge_T0_2->SetTitle("T_0 ")
tge_T0_3->SetTitle("T_0 ")
tge_T0_4->SetTitle("T_0 ")
TMultiGraph *tmg = new TMultiGraph()
tmg->Add(tge_T0_1)
tmg->Add(tge_T0_2)
tmg->Add(tge_T0_3)
tmg->Add(tge_T0_4)
tmg->Draw("AP")
gPad->BuildLegend()
gPad->Update()
TGraphErrors *tge_1 = (TGraphErrors*)_file0->Get("tge_Res_T_H")
TGraphErrors *tge_2 = (TGraphErrors*)_file1->Get("tge_Res_T_H")
TGraphErrors *tge_3 = (TGraphErrors*)_file2->Get("tge_Res_T_H")
TGraphErrors *tge_4 = (TGraphErrors*)_file3->Get("tge_Res_T_H")
tge_1->SetLineColor(1)
tge_2->SetLineColor(2)
tge_3->SetLineColor(4)
tge_4->SetLineColor(8)
tge_1->SetMarkerColor(1)
tge_2->SetMarkerColor(2)
tge_3->SetMarkerColor(4)
tge_4->SetMarkerColor(8)
tge_1->SetTitle("#sigma ")
tge_2->SetTitle("#sigma ")
tge_3->SetTitle("#sigma ")
tge_4->SetTitle("#sigma ")
TMultiGraph *tmg = new TMultiGraph()
tmg->Add(tge_1)
tmg->Add(tge_2)
tmg->Add(tge_3)
tmg->Add(tge_4)
tmg->Draw("AP")
gPad->BuildLegend()
gPad->Update()
*/
void make_plots_signal_shape(TString tune, TString row){
  TString basepath="Dati_PD/4.fit_results/";
  TFile *_file0 = new TFile(basepath+"PD_quarzo-att6-2f-f1-thr--15-T"+tune+"-F2000_DAQ_flat_fit_results_modelID_4_row_"+row+".root");
  TFile *_file1 = new TFile(basepath+"PD_quarzo-att6-2f-f1-thr--15-T"+tune+"-F2000_DAQ_flat_fit_results_modelID_2_row_"+row+".root");
  TFile *_file2 = new TFile(basepath+"PD_quarzo-att6-2f-f1-thr--15-T"+tune+"-F2000_DAQ_flat_fit_results_modelID_1_row_"+row+".root");





  TGraphErrors *tge_res_1 = (TGraphErrors*)_file0->Get("tge_Res_T_H");
  TGraphErrors *tge_res_2 = (TGraphErrors*)_file1->Get("tge_Res_T_H");
  TGraphErrors *tge_res_3 = (TGraphErrors*)_file2->Get("tge_Res_T_H");
  tge_res_1->SetLineColor(1);
  tge_res_2->SetLineColor(2);
  tge_res_3->SetLineColor(4);
  tge_res_1->SetMarkerColor(1);
  tge_res_2->SetMarkerColor(2);
  tge_res_3->SetMarkerColor(4);
  tge_res_1->SetTitle("#sigma Gaussian core + GaussExp tail");
  tge_res_2->SetTitle("#sigma GaussExp");
  tge_res_3->SetTitle("#sigma CB");
  TCanvas *c_res = new TCanvas("c_res","c_res");
  c_res->cd();
  TMultiGraph *tmg_res = new TMultiGraph();
  tmg_res->Add(tge_res_1);
  tmg_res->Add(tge_res_2);
  tmg_res->Add(tge_res_3);
  tmg_res->Draw("AP");
  tmg_res->GetXaxis()->SetTitle("column number");
  tmg_res->GetYaxis()->SetTitle("#sigma_{T} [ns]");
  tmg_res->GetYaxis()->SetRangeUser(0.05,0.08);
  c_res->SetGridy();
  gPad->BuildLegend();
  gPad->Update();
 
  TGraphErrors *tge_T0_1 = (TGraphErrors*)_file0->Get("tge_T0");
  TGraphErrors *tge_T0_2 = (TGraphErrors*)_file1->Get("tge_T0");
  TGraphErrors *tge_T0_3 = (TGraphErrors*)_file2->Get("tge_T0");
  tge_T0_1->SetLineColor(1);
  tge_T0_2->SetLineColor(2);
  tge_T0_3->SetLineColor(4);
  tge_T0_1->SetMarkerColor(1);
  tge_T0_2->SetMarkerColor(2);
  tge_T0_3->SetMarkerColor(4);
  tge_T0_1->SetTitle("T_0 Gaussian core + GaussExp tail");
  tge_T0_2->SetTitle("T_0 GaussExp");
  tge_T0_3->SetTitle("T_0 CB");
  TCanvas *c_T0 = new TCanvas("c_T0","c_T0");
  c_T0->cd();
  TMultiGraph *tmg_T0 = new TMultiGraph();
  tmg_T0->Add(tge_T0_1);
  tmg_T0->Add(tge_T0_2);
  tmg_T0->Add(tge_T0_3);
  tmg_T0->Draw("AP");
  tmg_T0->GetXaxis()->SetTitle("column number");
  tmg_T0->GetYaxis()->SetTitle("T_{0} [ns]");
  tmg_T0->GetYaxis()->SetRangeUser(-0.1,0.1);
  c_T0->SetGridy();
  gPad->BuildLegend();
  gPad->Update();

  TFile *f = new TFile("Dati_PD/plots_presentations/output_T"+tune+"_row_"+row+".root","recreate");
  f->cd();
  c_T0->Write();
  c_res->Write();
  tmg_T0->Write();
  tmg_res->Write();
  f->cd();
  f->Close();
  delete f;

  
}


void make_plots_Tune(TString row){
  TString basepath="Dati_PD/4.fit_results/";
  TFile *_file0 = new TFile(basepath+"PD_quarzo-att6-2f-f1-thr--15-T80-F2000_DAQ_flat_fit_results_modelID_2_row_"+row+".root");
  TFile *_file1 = new TFile(basepath+"PD_quarzo-att6-2f-f1-thr--15-T75-F2000_DAQ_flat_fit_results_modelID_2_row_"+row+".root");
  TFile *_file2 = new TFile(basepath+"PD_quarzo-att6-2f-f1-thr--15-T70-F2000_DAQ_flat_fit_results_modelID_2_row_"+row+".root");
  TFile *_file3 = new TFile(basepath+"PD_quarzo-att6-2f-f1-thr--15-T65-F2000_DAQ_flat_fit_results_modelID_2_row_"+row+".root");
  TFile *_file4 = new TFile(basepath+"PD_quarzo-att6-2f-f1-thr--15-T60-F2000_DAQ_flat_fit_results_modelID_2_row_"+row+".root");
  TFile *_file5 = new TFile(basepath+"PD_quarzo-att6-2f-f1-thr--15-T55-F2000_DAQ_flat_fit_results_modelID_2_row_"+row+".root");
  TFile *_file6 = new TFile(basepath+"PD_quarzo-att6-2f-f1-thr--15-T50-F2000_DAQ_flat_fit_results_modelID_2_row_"+row+".root");





  TGraphErrors *tge_res_1 = (TGraphErrors*)_file0->Get("tge_Res_T_H");
  TGraphErrors *tge_res_2 = (TGraphErrors*)_file1->Get("tge_Res_T_H");
  TGraphErrors *tge_res_3 = (TGraphErrors*)_file2->Get("tge_Res_T_H");
  TGraphErrors *tge_res_4 = (TGraphErrors*)_file3->Get("tge_Res_T_H");
  TGraphErrors *tge_res_5 = (TGraphErrors*)_file4->Get("tge_Res_T_H");
  TGraphErrors *tge_res_6 = (TGraphErrors*)_file5->Get("tge_Res_T_H");
  TGraphErrors *tge_res_7 = (TGraphErrors*)_file6->Get("tge_Res_T_H");
  tge_res_1->SetLineColor(1);
  tge_res_2->SetLineColor(2);
  tge_res_3->SetLineColor(4);
  tge_res_4->SetLineColor(8);
  tge_res_5->SetLineColor(7);
  tge_res_6->SetLineColor(6);
  tge_res_7->SetLineColor(5);
  tge_res_1->SetMarkerColor(1);
  tge_res_2->SetMarkerColor(2);
  tge_res_3->SetMarkerColor(4);
  tge_res_4->SetMarkerColor(8);
  tge_res_5->SetMarkerColor(7);
  tge_res_6->SetMarkerColor(6);
  tge_res_7->SetMarkerColor(5);
  tge_res_1->SetTitle("T80");
  tge_res_2->SetTitle("T75");
  tge_res_3->SetTitle("T70");
  tge_res_4->SetTitle("T65");
  tge_res_5->SetTitle("T60");
  tge_res_6->SetTitle("T55");
  tge_res_7->SetTitle("T50");
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
  tmg_res->GetXaxis()->SetTitle("column number");
  tmg_res->GetYaxis()->SetTitle("#sigma_{T} [ns]");
  //tmg_res->GetYaxis()->SetRangeUser(0.05,0.08);
  c_res->SetGridy();
  gPad->BuildLegend();
  gPad->Update();
 
  TGraphErrors *tge_T0_1 = (TGraphErrors*)_file0->Get("tge_T0");
  TGraphErrors *tge_T0_2 = (TGraphErrors*)_file1->Get("tge_T0");
  TGraphErrors *tge_T0_3 = (TGraphErrors*)_file2->Get("tge_T0");
  tge_T0_1->SetLineColor(1);
  tge_T0_2->SetLineColor(2);
  tge_T0_3->SetLineColor(4);
  tge_T0_1->SetMarkerColor(1);
  tge_T0_2->SetMarkerColor(2);
  tge_T0_3->SetMarkerColor(4);
  tge_T0_1->SetTitle("T_0 Gaussian core + GaussExp tail");
  tge_T0_2->SetTitle("T_0 GaussExp");
  tge_T0_3->SetTitle("T_0 CB");
  TCanvas *c_T0 = new TCanvas("c_T0","c_T0");
  c_T0->cd();
  TMultiGraph *tmg_T0 = new TMultiGraph();
  tmg_T0->Add(tge_T0_1);
  tmg_T0->Add(tge_T0_2);
  tmg_T0->Add(tge_T0_3);
  tmg_T0->Draw("AP");
  tmg_T0->GetXaxis()->SetTitle("column number");
  tmg_T0->GetYaxis()->SetTitle("T_{0} [ns]");
  tmg_T0->GetYaxis()->SetRangeUser(-0.1,0.1);
  c_T0->SetGridy();
  gPad->BuildLegend();
  gPad->Update();

  TFile *f = new TFile("Dati_PD/plots_presentations/output_Tune_scan_row_"+row+".root","recreate");
  f->cd();
  c_T0->Write();
  c_res->Write();
  tmg_T0->Write();
  tmg_res->Write();
  f->cd();
  f->Close();
  delete f;

  
}
