////////////////////////////////////////////////////
//// Added lines for branches exercises merging ////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
/// This is used to analyze Padova data, as well ///
///                as KEKCC data                 ///
////////////////////////////////////////////////////
/////////////// Author: A. Mordà ///////////////////
////////////////////////////////////////////////////

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

void run_data_production(TString input_raw_name="", bool flat_ntuple=false){
  TString file_origin="PD";

  cout<<"PRODUCING FLAT NTUPLE INPUT"<<endl;
  if(flat_ntuple) cscoper_DAQ("",input_raw_name,"");//produce_flat_ntuples
  
  cout<<"PRODUCING HISTOGRAMS FOR FITTING"<<endl;
  for(int i=1;i<=4;i++){
    TString file_kind='';
    if(i==1){file_kind="recreate";}else{file_kind="update";}
    make_PD_data_histos_column("",input_raw_name,"",file_kind,i);
  }
  
  cout<<"PRODUCING PMT SPECTRA"<<endl;
  make_pmt_plots("",input_raw_name,3);
  
}



int cscoper_DAQ(TString origin_path, TString fname, TString output_path){

  
  float _times_DAQ[32];
  float _amps_DAQ[32];
  int _chans_DAQ[32];
  float _temp_DAQ;
  int _nchs_DAQ;
  UInt_t _evtnum_DAQ;
  
  if(origin_path==""){origin_path="D:/Padova/TOP/Time_calibration_codes/Dati_PD/1.raw_input/";}
  if(output_path==""){output_path="D:/Padova/TOP/Time_calibration_codes/Dati_PD/2.flat_input/";}
  
  float time_DAQ(-99);
  float amp_DAQ(-99);
  int channel_DAQ(-99);
  float temp_DAQ(-99);
  int nchs_DAQ(-99);
  UInt_t evtnum_DAQ(-99);
  Double_t numevts(-99);



  TFile *f_input = new TFile(origin_path+fname+".root");
  TTree *_tree_DAQ=(TTree*)f_input->Get("times");
  _tree_DAQ->SetBranchAddress("times",&(_times_DAQ[0]));
  _tree_DAQ->SetBranchAddress("amps",&(_amps_DAQ[0]));
  _tree_DAQ->SetBranchAddress("chans",&(_chans_DAQ[0]));
  _tree_DAQ->SetBranchAddress("nchans",&_nchs_DAQ);
  _tree_DAQ->SetBranchAddress("evtnum",&_evtnum_DAQ);
  _tree_DAQ->SetBranchAddress("temp",&_temp_DAQ);


  
  TFile *f_data = new TFile(output_path+fname+"_DAQ_flat.root","recreate");
  TTree *tree_DAQ = new TTree("times","times");
  tree_DAQ->Branch("time",&time_DAQ,"time/F");
  tree_DAQ->Branch("amp",&amp_DAQ,"amp/F");
  tree_DAQ->Branch("channel",&channel_DAQ,"channel/I");
  tree_DAQ->Branch("nchannels",&nchs_DAQ,"nchannels/I");
  tree_DAQ->Branch("evtnum",&evtnum_DAQ,"evtnum/I");
  tree_DAQ->Branch("temp",&temp_DAQ,"temp/F");
  
   
  Int_t _nevents = _tree_DAQ->GetEntries();
  
  Int_t _curevent;
  cout<<"Input \".root\" file opened"<<endl;


  for( _curevent =0; _curevent <_nevents; _curevent++) {
    _tree_DAQ->GetEntry(_curevent);
    for(int _ichan=0; _ichan<_nchs_DAQ;_ichan++){
      if(60<_times_DAQ[_ichan]&&_times_DAQ[_ichan]<72){
      //if(0<_times_DAQ[_ichan]&&_times_DAQ[_ichan]<100){
	time_DAQ=_times_DAQ[_ichan];
	channel_DAQ=_chans_DAQ[_ichan];
	amp_DAQ=_amps_DAQ[_ichan];
	nchs_DAQ=_nchs_DAQ;
	evtnum_DAQ=_evtnum_DAQ;
	temp_DAQ=_temp_DAQ;
	tree_DAQ->Fill();
	if(_evtnum_DAQ>numevts) {numevts=_evtnum_DAQ;}
      }
    }
    
  }
  cout<<"n trigger : "<<numevts<<endl;
  TH1D *h_ntriggers = new TH1D("h_ntriggers","h_ntriggers",1,0,1);
  h_ntriggers->SetBinContent(1,numevts);
  Double_t n_triggers = h_ntriggers->GetBinContent(1);
  cout<<"n triggers : "<<n_triggers<<endl;
  f_data->cd();
  tree_DAQ->Write();
  h_ntriggers->Write();
  delete h_ntriggers;
  f_data->cd();
  f_data->Close();
  delete f_data;
  return 1;
}



void  make_PD_data_histos_column(TString input_path, TString file_name, TString output_path, TString out_file_status, int my_column){
  if (input_path=="") input_path="Dati_PD/2.flat_input/";
  if (output_path=="") output_path="Dati_PD/3.column_data/";
  if(out_file_status!="recreate"&& out_file_status!="update"){
    cout<<"ERROR: incorrect kind of output file: recreate or update, output file is being recreated from scretch"<<endl;
    out_file_status="recreate";
  }

  /////// INITIALIZING INPUT FILE
  TString inp_f =input_path+file_name+"_DAQ_flat.root";
  TFile *file_input = new TFile(inp_f);
  cout<<"input file opened : "<<inp_f<<endl;
  TString tree_input="times"; 
  TTree *t_input = (TTree*)file_input->Get(tree_input);

  /////// INITIALIZING OUTPUT FILE FOR DATA
  TString columnID =Form("col_%i",my_column);
  TString of=output_path+file_name+"_data_histos.root"; 
  TFile *f_data = new TFile(of,out_file_status);
  cout<<"output file created (updating) : "<<of<<endl;
  f_data->cd();
  f_data->mkdir(Form("column_%i",my_column));
  f_data->cd(Form("column_%i",my_column));

  TH1D *h_ntrig = file_input->Get("h_ntriggers");
  
  Double_t n_triggers = h_ntrig->GetBinContent(1);
  cout<<"n triggers : "<<n_triggers<<endl;

  TH1D *h_yields = new TH1D("h_yields","event yields",8,0,9);
  
  float upper_time;
  float lower_time;
  lower_time=60;
  upper_time=72; 
  //upper_time=60; lower_time=72;
  int my_pixelID_data=-9; ////////For PD data
  
  int my_pixelID=-9;
  my_pixelID=my_column;
  int KEK_PMT=13;
  int my_row=-9;
  int index_pmt=0;
  for(int g=1; g<=8;g++){
    cout<<"doing "<<g<<"th row"<<endl;
    if(g>4)index_pmt=1;
    my_pixelID=my_column;
    if(g<=4){my_row=g;}else{my_row=g-4;}
    my_pixelID_data=(my_column)*4+20*index_pmt-g;///////actually the readout channel for PD testbench system data
    
 
    TCut cut = Form("amp>15&&60<time&&time<72&&channel==%i",my_pixelID_data); //For PD
    TCut cut_16ch = Form("nchannels==16"); //For PD
    TCut cut_32ch = Form("nchannels==32"); //For PD


    /////// THIS IS THE HISTOGAM USED TO RESCALE THE TIME DISTRIBUTION TO ZERO
    TH1D *h_temp = new TH1D("h_temp",Form("time_col_%i_row_%i",my_column,g),500,lower_time,upper_time);
    t_input->Project("h_temp","time",cut);
    h_yields->SetBinContent(g,t_input->GetEntries(cut));
    Float_t max_bin = h_temp->GetMaximumBin();
    TAxis *xaxis = h_temp->GetXaxis(); 
    Double_t max_pos = xaxis->GetBinCenter(max_bin);
    //h_temp->Scale(1.0/n_triggers);

    //////// THIS IS THE HISTOGRAM THAT WILL BE FITTED
    Float_t upper_bound_hist=2;
    Int_t n_bins = 200;
    TH1D *h_time = new TH1D("h_time","Time [ns]",n_bins,-1.5,upper_bound_hist);
    t_input->Project("h_time",Form("time-%f",max_pos),cut);
    TH1D *h_time_16 = new TH1D("h_time_16","Time [ns]",n_bins,-1.5,upper_bound_hist);
    t_input->Project("h_time_16",Form("time-%f",max_pos),cut&&cut_16ch);
    TH1D *h_time_32 = new TH1D("h_time_32","Time [ns]",n_bins,-1.5,upper_bound_hist);
    t_input->Project("h_time_32",Form("time-%f",max_pos),cut&&cut_32ch);
    
    /*
    float mean_tmp = h_time->GetMean();
    if(mean_tmp>0&&my_row>4){
      h_time->Reset();
      TH1D *h_tmp_1 = new TH1D("h_tmp_1","Time [ns]",n_bins,0,upper_bound_hist);
      t_input->Project("h_tmp_1",Form("time-%f-%f",max_pos,mean_tmp),cut);
      Float_t max_bin_1 = h_tmp_1->GetMaximumBin();
      TAxis *xaxis_1 = h_tmp_1->GetXaxis(); 
      Double_t max_pos_1 = xaxis_1->GetBinCenter(max_bin_1);
      t_input->Project("h_time",Form("time-%f",max_pos_1),cut);
      delete h_tmp_1;
    }
    */

    //// THIS IS THE AMPLITUDE HISTOGRAM
    TString variable_amp="amp";
    Float_t up_limit_amp =300;
    TH1D *h_amp = new TH1D("h_amp","Max Amplitude [ADC counts]",n_bins,0,300);
    t_input->Project("h_amp",variable_amp,cut);
    TH1D *h_amp_16 = new TH1D("h_amp_16","Max Amplitude [ADC counts]",n_bins,0,300);
    t_input->Project("h_amp_16",variable_amp,cut&&cut_16ch);
    TH1D *h_amp_32 = new TH1D("h_amp_32","Max Amplitude [ADC counts]",n_bins,0,300);
    t_input->Project("h_amp_32",variable_amp,cut&&cut_32ch);
   
    f_data->mkdir(Form("column_%i/histos_%i_%i",my_column,my_column,g));
    f_data->cd(Form("column_%i/histos_%i_%i",my_column,my_column,g));
    h_amp->SetName(Form("amp_col_%i_row_%i",my_column,g));
    h_amp->Write();
    h_temp->SetName(Form("time_col_%i_row_%i",my_column,g));
    h_temp->Write();
    h_time->Write();
    h_amp_16->Write();
    h_time_16->Write();
    h_amp_32->Write();
    h_time_32->Write();
    f_data->cd();
    delete h_temp;
    delete h_time;
    delete h_amp;
    delete h_time_16;
    delete h_amp_16;
    delete h_time_32;
    delete h_amp_32;
    delete h_ntrig;
  }
  f_data->cd(Form("column_%i",my_column));
  h_yields->Write();
  f_data->cd();
  f_data->Close();
  delete f_data;
}

void  make_KEK_data_histos_column(TString input_path, TString file_name, TString output_path){

  if (input_path=="") input_path="Dati_PD/2.flat_input/";
  if (output_path=="") output_path="Dati_PD/3.column_data/";
  
  TFile *f_input = new TFile(input_path+file_name+".root");
  TTree *t_input = (TTree*)f_input->Get("laser");
  TFile *f_output = new TFile(output_path+file_name+"_data_histos.root","recreate");
  TH1D* h[16][64][8];
  TH1D* h_amp[16][64][8];
  float low_time = 65;
  float up_time = 70;

  
  for(int j=0;j<16;j++){
    for(int jj=0;jj<64;jj++){
      for(int jjj=0;jjj<8;jjj++){
	h[j][jj][jjj] = new TH1D(Form("time_slot_%i_col_%i_row_%i",j+1,jj+1,jjj+1),Form("time_slot_%i_col_%i_row_%i",j+1,jj+1,jjj+1),200,low_time,up_time);
	h_amp[j][jj][jjj] = new TH1D(Form("amp_slot_%i_col_%i_row_%i",j+1,jj+1,jjj+1),Form("amp_slot_%i_col_%i_row_%i",j+1,jj+1,jjj+1),1000,0,1000);
      }
    }
  }
  cout<<"initialized histograms"<<endl;
    
  float time=0;
  int valuepeak=0;
  int column=-9;
  int row=-9;
  int slot=-9;
  int quality=-9;
  t_input->SetBranchAddress("time",&time);
  t_input->SetBranchAddress("valuepeak",&valuepeak);
  t_input->SetBranchAddress("column",&column);
  t_input->SetBranchAddress("row",&row);
  t_input->SetBranchAddress("slot",&slot);
  t_input->SetBranchAddress("quality",&quality);

  int nentries=t_input->GetEntries();
  for (int i=0;i<nentries;i++){
    if(i%1000000==0)cout<<i<<endl;
    t_input->GetEntry(i);
    if(quality!=1) continue;
    h[slot-1][column-1][row-1]->Fill(time);
    h_amp[slot-1][column-1][row-1]->Fill(valuepeak);
  }
  cout<<"filled histograms"<<endl;
  cout<<"now writing histograms"<<endl;
  f_output->cd();
  for(int j=0;j<16;j++){
    f_output->mkdir(Form("slot_%i",j+1));
    f_output->cd(Form("slot_%i",j+1));
    for(int jj=0;jj<64;jj++){
      TH1D *h_yields = new TH1D("h_yields","event yields",8,0,9);
      for(int jjj=0;jjj<8;jjj++){
	f_output->mkdir(Form("slot_%i/column_%i/histos_%i_%i",j+1,jj+1,jj+1,jjj+1));
	f_output->cd(Form("slot_%i/column_%i/histos_%i_%i",j+1,jj+1,jj+1,jjj+1));
	h[j][jj][jjj]->Write();
	h_amp[j][jj][jjj]->Write();
	h_yields->SetBinContent(jjj+1,h[j][jj][jjj]->GetEntries());
       }
      f_output->cd(Form("slot_%i/column_%i/",j+1,jj+1));
      h_yields->Write();
     }
  }

  f_output->Close();

  delete f_output;
  
}

void  make_MC_histos_column(TString output_path, int pmt_column, int pmt_pos){
  if (output_path=="") output_path="Dati_PD/3.column_data/";



  int my_column=-9;
  my_column=pmt_column+(pmt_pos-1)*4;
  ///////INITIALIZING MC FILES
  TFile *file_input_MC = new TFile("Dati_PD/1.raw_input/ana_laser_s01_0reso_500k.root");
  TTree *tree_MC = (TTree*)file_input_MC->Get("laser");
  TFile *file_input_MC_ring = new TFile("Dati_PD/1.raw_input/ana_laser_s01_0reso_ring_500k.root");
  TTree *tree_MC_ring = (TTree*)file_input_MC_ring->Get("laser");
  cout<<"input files initialized "<<endl;
  /////INITIALIZING output file
  TFile *f_data = new TFile(output_path+"MC_inputs_test.root","recreate");
  cout<<"output file initialized "<<endl;



  
  cout<<"initializing histograms"<<endl;
  TH1D* h_nominal[64][8];
  TH1D* h_nominal_f1[64][8];
  TH1D* h_nominal_f2[64][8];
  TH1D* h_nominal_f3[64][8];
  TH1D* h_nominal_f4[64][8];
  TH1D* h_nominal_f5[64][8];
  TH1D* h_nominal_f6[64][8];
  TH1D* h_nominal_f7[64][8];
  TH1D* h_nominal_f8[64][8];
  TH1D* h_nominal_f9[64][8];
  
  TH1D* h_ring[64][8];
  TH1D* h_ring_f1[64][8];
  TH1D* h_ring_f2[64][8];
  TH1D* h_ring_f3[64][8];
  TH1D* h_ring_f4[64][8];
  TH1D* h_ring_f5[64][8];
  TH1D* h_ring_f6[64][8];
  TH1D* h_ring_f7[64][8];
  TH1D* h_ring_f8[64][8];
  TH1D* h_ring_f9[64][8];
  for(int jj=0;jj<64;jj++){
    for(int jjj=0;jjj<8;jjj++){
      h_nominal[jj][jjj] = new TH1D(Form("time_nom_col_%i_row_%i",jj+1,jjj+1),Form("time_nom_col_%i_row_%i",jj+1,jjj+1),200,0.3,1);
      h_nominal_f1[jj][jjj] = new TH1D(Form("time_nom_f1_col_%i_row_%i",jj+1,jjj+1),Form("time_nom_f1_col_%i_row_%i",jj+1,jjj+1),200,0.3,1);
      h_nominal_f2[jj][jjj] = new TH1D(Form("time_nom_f2_col_%i_row_%i",jj+1,jjj+1),Form("time_nom_f2_col_%i_row_%i",jj+1,jjj+1),200,0.3,1);
      h_nominal_f3[jj][jjj] = new TH1D(Form("time_nom_f3_col_%i_row_%i",jj+1,jjj+1),Form("time_nom_f3_col_%i_row_%i",jj+1,jjj+1),200,0.3,1);
      h_nominal_f4[jj][jjj] = new TH1D(Form("time_nom_f4_col_%i_row_%i",jj+1,jjj+1),Form("time_nom_f4_col_%i_row_%i",jj+1,jjj+1),200,0.3,1);
      h_nominal_f5[jj][jjj] = new TH1D(Form("time_nom_f5_col_%i_row_%i",jj+1,jjj+1),Form("time_nom_f5_col_%i_row_%i",jj+1,jjj+1),200,0.3,1);
      h_nominal_f6[jj][jjj] = new TH1D(Form("time_nom_f6_col_%i_row_%i",jj+1,jjj+1),Form("time_nom_f6_col_%i_row_%i",jj+1,jjj+1),200,0.3,1);
      h_nominal_f7[jj][jjj] = new TH1D(Form("time_nom_f7_col_%i_row_%i",jj+1,jjj+1),Form("time_nom_f7_col_%i_row_%i",jj+1,jjj+1),200,0.3,1);
      h_nominal_f8[jj][jjj] = new TH1D(Form("time_nom_f8_col_%i_row_%i",jj+1,jjj+1),Form("time_nom_f8_col_%i_row_%i",jj+1,jjj+1),200,0.3,1);
      h_nominal_f9[jj][jjj] = new TH1D(Form("time_nom_f9_col_%i_row_%i",jj+1,jjj+1),Form("time_nom_f9_col_%i_row_%i",jj+1,jjj+1),200,0.3,1);
      h_ring[jj][jjj] = new TH1D(Form("time_ring_col_%i_row_%i",jj+1,jjj+1),Form("time_ring_col_%i_row_%i",jj+1,jjj+1),200,0.3,1);
      h_ring_f1[jj][jjj] = new TH1D(Form("time_ring_f1_col_%i_row_%i",jj+1,jjj+1),Form("time_ring_f1_col_%i_row_%i",jj+1,jjj+1),200,0.3,1);
      h_ring_f2[jj][jjj] = new TH1D(Form("time_ring_f2_col_%i_row_%i",jj+1,jjj+1),Form("time_ring_f2_col_%i_row_%i",jj+1,jjj+1),200,0.3,1);
      h_ring_f3[jj][jjj] = new TH1D(Form("time_ring_f3_col_%i_row_%i",jj+1,jjj+1),Form("time_ring_f3_col_%i_row_%i",jj+1,jjj+1),200,0.3,1);
      h_ring_f4[jj][jjj] = new TH1D(Form("time_ring_f4_col_%i_row_%i",jj+1,jjj+1),Form("time_ring_f4_col_%i_row_%i",jj+1,jjj+1),200,0.3,1);
      h_ring_f5[jj][jjj] = new TH1D(Form("time_ring_f5_col_%i_row_%i",jj+1,jjj+1),Form("time_ring_f5_col_%i_row_%i",jj+1,jjj+1),200,0.3,1);
      h_ring_f6[jj][jjj] = new TH1D(Form("time_ring_f6_col_%i_row_%i",jj+1,jjj+1),Form("time_ring_f6_col_%i_row_%i",jj+1,jjj+1),200,0.3,1);
      h_ring_f7[jj][jjj] = new TH1D(Form("time_ring_f7_col_%i_row_%i",jj+1,jjj+1),Form("time_ring_f7_col_%i_row_%i",jj+1,jjj+1),200,0.3,1);
      h_ring_f8[jj][jjj] = new TH1D(Form("time_ring_f8_col_%i_row_%i",jj+1,jjj+1),Form("time_ring_f8_col_%i_row_%i",jj+1,jjj+1),200,0.3,1);
      h_ring_f9[jj][jjj] = new TH1D(Form("time_ring_f9_col_%i_row_%i",jj+1,jjj+1),Form("time_ring_f9_col_%i_row_%i",jj+1,jjj+1),200,0.3,1);
    }
  }  
  cout<<"initialized histograms"<<endl;
  double time=0;
  int column=-9;
  int row=-9;
  int fiberNo=-9;
  tree_MC->SetBranchAddress("propTime",&time);
  tree_MC->SetBranchAddress("column",&column);
  tree_MC->SetBranchAddress("row",&row);
  tree_MC->SetBranchAddress("fiberNo",&fiberNo);

  int nentries=tree_MC->GetEntries();
  for (int i=0;i<nentries;i++){
    if(i%1000000==0)cout<<i<<endl;
    tree_MC->GetEntry(i);
    h_nominal[column-1][row-1]->Fill(time);
    if(fiberNo==1){h_nominal_f1[column-1][row-1]->Fill(time);}
    else if(fiberNo==2){h_nominal_f2[column-1][row-1]->Fill(time);}
    else if(fiberNo==3){h_nominal_f3[column-1][row-1]->Fill(time);}
    else if(fiberNo==4){h_nominal_f4[column-1][row-1]->Fill(time);}
    else if(fiberNo==5){h_nominal_f5[column-1][row-1]->Fill(time);}
    else if(fiberNo==6){h_nominal_f6[column-1][row-1]->Fill(time);}
    else if(fiberNo==7){h_nominal_f7[column-1][row-1]->Fill(time);}
    else if(fiberNo==8){h_nominal_f8[column-1][row-1]->Fill(time);}
    else if(fiberNo==9){h_nominal_f9[column-1][row-1]->Fill(time);}
  }
  cout<<"filled histograms nominal model"<<endl;
  
  tree_MC_ring->SetBranchAddress("propTime",&time);
  tree_MC_ring->SetBranchAddress("column",&column);
  tree_MC_ring->SetBranchAddress("row",&row);
  tree_MC_ring->SetBranchAddress("fiberNo",&fiberNo);
  int nentries_ring=tree_MC_ring->GetEntries();
  for (int i=0;i<nentries_ring;i++){
    if(i%1000000==0)cout<<i<<endl;
    tree_MC_ring->GetEntry(i);
    h_ring[column-1][row-1]->Fill(time);
    if(fiberNo==1){h_ring_f1[column-1][row-1]->Fill(time);}
    else if(fiberNo==2){h_ring_f2[column-1][row-1]->Fill(time);}
    else if(fiberNo==3){h_ring_f3[column-1][row-1]->Fill(time);}
    else if(fiberNo==4){h_ring_f4[column-1][row-1]->Fill(time);}
    else if(fiberNo==5){h_ring_f5[column-1][row-1]->Fill(time);}
    else if(fiberNo==6){h_ring_f6[column-1][row-1]->Fill(time);}
    else if(fiberNo==7){h_ring_f7[column-1][row-1]->Fill(time);}
    else if(fiberNo==8){h_ring_f8[column-1][row-1]->Fill(time);}
    else if(fiberNo==9){h_ring_f9[column-1][row-1]->Fill(time);}
  }
  cout<<"filled histograms ring model"<<endl;

  f_data->cd();
  for(int jj=0;jj<64;jj++){
    for(int jjj=0;jjj<8;jjj++){
      f_data->mkdir(Form("column_%i/histos_%i_%i",jj+1,jj+1,jjj+1));
      f_data->cd(Form("column_%i/histos_%i_%i",jj+1,jj+1,jjj+1));  
      h_nominal[jj][jjj]->Write(); 
      h_nominal_f1[jj][jjj]->Write();
      h_nominal_f2[jj][jjj]->Write();
      h_nominal_f3[jj][jjj]->Write();
      h_nominal_f4[jj][jjj]->Write();
      h_nominal_f5[jj][jjj]->Write();
      h_nominal_f6[jj][jjj]->Write();
      h_nominal_f7[jj][jjj]->Write();
      h_nominal_f8[jj][jjj]->Write();
      h_nominal_f9[jj][jjj]->Write();
      h_ring[jj][jjj]->Write();
      h_ring_f1[jj][jjj]->Write();
      h_ring_f2[jj][jjj]->Write();
      h_ring_f3[jj][jjj]->Write();
      h_ring_f4[jj][jjj]->Write();
      h_ring_f5[jj][jjj]->Write();
      h_ring_f6[jj][jjj]->Write();
      h_ring_f7[jj][jjj]->Write();
      h_ring_f8[jj][jjj]->Write();
      h_ring_f9[jj][jjj]->Write();
    }
  }
  f_data->cd();
  f_data->Close();
  delete f_data;




}




void make_pmt_plots(TString input_path, TString input_filebasename, int pmt_pos, int slotID=-1){

  bool save = true;
  if (input_path=="") input_path="Dati_PD/3.column_data/";
  TString channel_cut="";//can be "_16" or "_32"
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



  /////// READ INPUT DATA
  TString input_filename=input_filebasename;
  input_filename=input_path+input_filename+"_data_histos";
  TFile *f = new TFile(input_filename+".root");
  ////// READ MC INPUTS
  TFile *f_MC = new TFile(input_path+"MC_inputs.root");
  int MC_row=-99;
  TH1D* h_MC[5][9];
  TH1D* h_MC_ring[5][9];
  TString slot_ID= "";
  TString slot_ID_undsc= "";
  TString slot_ID_bs= "";
  int shift_position=0; //when looking at KEK data to shift the position of the PMT while keeping the PD column numbering scheme
  
  if(slotID>0) {
    slot_ID = Form("slot_%i",slotID);
    slot_ID_undsc = "_";
    slot_ID_bs= "/";
    shift_position=1;
  }
  
  for(int i=1;i<=4;i++){
    int pmt_index=1;
    int pmt_row=-1;
    TH1D *h_yields = (TH1D*)f->Get(slot_ID+slot_ID_bs+Form("column_%i/h_yields",i+shift_position*(pmt_pos-1)*4));
    for(int g=1; g<=8;g++){
      if(g<=4){pmt_row=g;} else{pmt_row=g-4;}
      int row_id_mc=9-g;
      int MC_column=0;
      TString histo_name_base_slot=slot_ID+slot_ID_bs;
      TString histo_name_base=Form("column_%i/histos_%i_%i/",i+shift_position*(pmt_pos-1)*4,i+shift_position*(pmt_pos-1)*4,g);
      TString histo_name_prefix_time = "time_";
      TString histo_name_prefix_amp = "amp_";
      TString histo_name_slot=slot_ID+slot_ID_undsc;
      TString histo_name_colrow= Form("col_%i_row_%i",i+shift_position*(pmt_pos-1)*4,g);
      TString histoname_time =histo_name_base_slot+histo_name_base+histo_name_prefix_time+histo_name_slot+histo_name_colrow+channel_cut;
      TString histoname_amp  =histo_name_base_slot+histo_name_base+histo_name_prefix_amp+histo_name_slot+histo_name_colrow+channel_cut;
      
      cout<<histoname_time<<endl;
      MC_column=i+(pmt_pos-1)*4;
      h[i][g]=(TH1D*)f->Get(histoname_time);
      //§ h_MC[i][g]=(TH1D*)f_MC->Get(Form("column_%i/histos_%i_%i/h_MC_tot",MC_column,MC_column,row_id_mc));
      //§ h_MC_ring[i][g]=(TH1D*)f_MC->Get(Form("column_%i/histos_%i_%i/h_MC_ring_tot",MC_column,MC_column,row_id_mc));
      h_amp[i][g]=(TH1D*)f->Get(histoname_amp);
      int canvasID=i+4*(4-pmt_row);
      if(g<=4){
	pmt_down->cd(canvasID);
	h[i][g]->Draw();
	//§ h_MC[i][g]->Scale(h[i][g]->GetMaximum()/h_MC[i][g]->GetMaximum());
	//§ h_MC_ring[i][g]->Scale(h[i][g]->GetMaximum()/h_MC_ring[i][g]->GetMaximum());
	//§ h_MC_ring[i][g]->SetLineColor(2);
	//§ h_MC_ring[i][g]->SetMarkerSize(0);
	//§ h_MC_ring[i][g]->SetMarkerColor(2);
	//§ h_MC[i][g]->Draw("same");
	//§ h_MC_ring[i][g]->Draw("same");
	gPad->BuildLegend();
	pmts->cd(16+canvasID);
	h[i][g]->Draw();
	pmts_amp->cd(16+canvasID);
	h_amp[i][g]->Draw("same");
      }else {
	pmt_up->cd(canvasID);
	h[i][g]->Draw();
	//§ h_MC[i][g]->Scale(h[i][g]->GetMaximum()/h_MC[i][g]->GetMaximum());
	//§ h_MC_ring[i][g]->Scale(h[i][g]->GetMaximum()/h_MC_ring[i][g]->GetMaximum());
	//§ h_MC[i][g]->Draw("same");
	//§ h_MC_ring[i][g]->SetLineColor(2);
	//§ h_MC_ring[i][g]->SetMarkerSize(0);
	//§ h_MC_ring[i][g]->SetMarkerColor(2);
	//§ h_MC_ring[i][g]->Draw("same");	
	gPad->BuildLegend();
	pmts->cd(canvasID);
	h[i][g]->Draw();
	pmts_amp->cd(canvasID);
	h_amp[i][g]->Draw("same");
      }
      h_yield_map->SetBinContent(i,g,h_yields->GetBinContent(g));
    }
    delete h_yields;
  }
  TCanvas *c_map = new TCanvas("pmts_occupancy_"+input_filebasename,"pmts_occupancy_"+input_filebasename);
  h_yield_map->DrawNormalized("colz");
  if(save){
    TFile *f_out = new TFile(input_path+"/pmt_plots/"+input_filebasename+Form("_pmt_pos%i_plots.root",pmt_pos),"recreate");
    f_out->cd();
    c_map->Write();
    pmts_amp->Write();
    pmts->Write();
    pmt_up->Write();
    pmt_down->Write();
    h_yield_map->Write();
    f_out->cd();
    f_out->Close();
    delete f_out;
  }
  delete h_yield_map;
}



perform_yields_shapes_comparison(TString input_path, TString file_in_0, TString file_in_1, TString file_in_2, TString file_in_3, TString file_in_4=""){

  /*
    file_0 is for the target comparison
    One histograms matrix per channel, enumerating from 1 to keep things more readeable
    ex. perform_yields_shapes_comparison("",
    "quarzo-06nov17-PMTpos3-3f-thr--20-T50-F1982-1",
    "quarzo-06nov17-PMTpos3-f1-thr--20-T50-F1982-1",
    "quarzo-06nov17-PMTpos3-f2-thr--20-T50-F1982-1",
    "quarzo-06nov17-PMTpos3-f3-thr--20-T50-F1982-1",
    "quarzo-06nov17-PMTpos3-3f-thr--20-T50-F1982-1"
    )

    perform_yields_shapes_comparison("","quarzo-06nov17-PMTpos3-3f-thr--20-T50-F1982-1","quarzo-06nov17-PMTpos3-f1-thr--20-T50-F1982-1","quarzo-06nov17-PMTpos3-f2-thr--20-T50-F1982-1",    "quarzo-06nov17-PMTpos3-f3-thr--20-T50-F1982-1","quarzo-06nov17-PMTpos3-3f-thr--20-T50-F1982-1")
  */
  
  if (input_path=="") input_path="Dati_PD/3.column_data/";
  TFile *f_input_0 = new TFile(input_path+file_in_0+"_data_histos.root");
  TFile *f_input_1 = new TFile(input_path+file_in_1+"_data_histos.root");
  TFile *f_input_2 = new TFile(input_path+file_in_2+"_data_histos.root");
  TFile *f_input_3 = new TFile(input_path+file_in_3+"_data_histos.root");
  TFile *f_input_4 = new TFile(input_path+file_in_4+"_data_histos.root");
  TFile *f_ntr_input_0 = new TFile(input_path+"../2.flat_input/"+file_in_0+"_DAQ_flat.root");
  TFile *f_ntr_input_1 = new TFile(input_path+"../2.flat_input/"+file_in_1+"_DAQ_flat.root");
  TFile *f_ntr_input_2 = new TFile(input_path+"../2.flat_input/"+file_in_2+"_DAQ_flat.root");
  TFile *f_ntr_input_3 = new TFile(input_path+"../2.flat_input/"+file_in_3+"_DAQ_flat.root");
  TFile *f_ntr_input_4 = new TFile(input_path+"../2.flat_input/"+file_in_4+"_DAQ_flat.root");
  TH1 *h_ntr_0 = f_ntr_input_0->Get("h_ntriggers");
  TH1 *h_ntr_1 = f_ntr_input_1->Get("h_ntriggers");
  TH1 *h_ntr_2 = f_ntr_input_2->Get("h_ntriggers");
  TH1 *h_ntr_3 = f_ntr_input_3->Get("h_ntriggers");
  TH1 *h_ntr_4 = f_ntr_input_4->Get("h_ntriggers");
  Double_t ntrigs[5];
  ntrigs[0] = h_ntr_0->GetBinContent(1);
  ntrigs[1] = h_ntr_1->GetBinContent(1);
  ntrigs[2] = h_ntr_2->GetBinContent(1);
  ntrigs[3] = h_ntr_3->GetBinContent(1);
  ntrigs[4] = h_ntr_4->GetBinContent(1);


  TH1D *h_0[5][9];
  TH1D *h[5][5][9];

  TCanvas *pmts = new TCanvas("pmts","pmts");
  pmts->Divide(4,8);
  TCanvas *pmts_up = new TCanvas("pmts_up","pmts_up");
  pmts_up->Divide(4,4);
  TCanvas *pmts_down = new TCanvas("pmts_down","pmts_down");
  pmts_down->Divide(4,4);



  
  for(int row=1;row<=8;row++){
    for(int column=1;column<=4;column++){
      h_0[column][row]  = (TH1D*)f_input_0->Get(Form("column_%i/histos_%i_%i/time_col_%i_row_%i",column,column,row,column,row));
      h[1][column][row] = (TH1D*)f_input_1->Get(Form("column_%i/histos_%i_%i/time_col_%i_row_%i",column,column,row,column,row));
      h[2][column][row] = (TH1D*)f_input_2->Get(Form("column_%i/histos_%i_%i/time_col_%i_row_%i",column,column,row,column,row));
      h[3][column][row] = (TH1D*)f_input_3->Get(Form("column_%i/histos_%i_%i/time_col_%i_row_%i",column,column,row,column,row));
      h[4][column][row] = (TH1D*)f_input_4->Get(Form("column_%i/histos_%i_%i/time_col_%i_row_%i",column,column,row,column,row));
      h[0][column][row] = (TH1D*)h[1][column][row]->Clone();
    }
  }
  for(int row=1;row<=8;row++){
    for(int column=1;column<=4;column++){
      h_0[column][row]->SetName(Form("full_row_%i__column_%i",row,column));
      h_0[column][row]->SetTitle(Form("full_row_%i__column_%i",row,column));
      h_0[column][row]->GetXaxis()->SetTitle("Time [ns]");
      h_0[column][row]->GetXaxis()->SetLabelFont(70);
      h_0[column][row]->GetXaxis()->SetLabelSize(0.05);
      h_0[column][row]->GetXaxis()->SetTitleFont(70);
      h_0[column][row]->GetXaxis()->SetTitleSize(0.05);
      h_0[column][row]->GetYaxis()->SetTitle("A.U");
      h_0[column][row]->GetYaxis()->SetLabelFont(70);
      h_0[column][row]->GetYaxis()->SetLabelSize(0.05);
      h_0[column][row]->GetYaxis()->SetTitleFont(70);
      h_0[column][row]->GetYaxis()->SetTitleSize(0.05);
      for(int ll = 0; ll<=4; ll++){
	h[ll][column][row]->SetName(Form("full_row_%i__column_%i__fib%i",row,column,ll));
	h[ll][column][row]->SetTitle(Form("full_row_%i__column_%i__fib%i",row,column,ll));
	h[ll][column][row]->GetXaxis()->SetTitle("Time [ns]");
	h[ll][column][row]->GetXaxis()->SetLabelFont(70);
	h[ll][column][row]->GetXaxis()->SetLabelSize(0.05);
	h[ll][column][row]->GetYaxis()->SetTitle("A.U");
	h[ll][column][row]->GetYaxis()->SetLabelFont(70);
	h[ll][column][row]->GetYaxis()->SetLabelSize(0.05);
	//cout<<"h["<<ll<<"] integral : "<<h[ll][column][row]->GetIntegral()<<"  , h["<<ll<<"] entries : " <<h[ll][column][row]->GetEntries()<<endl;
      }
      
    }
  }
  TH1D *h_spread = new TH1D("h_spread","h_spread",50,-0.5,0.5);
  TH3D *h3_spread = new TH3D("h3_spread","h3_spread",4,1,5,8,1,9,50,-0.5,0.5);
  TH2D *h2_spread = new TH2D("h2_spread","#frac{n. obs - n. exp}{n. obs}",4,1,5,8,1,9);
  TH2D *h22_spread = new TH2D("h22_spread","#frac{n. exp_{tot} - n. exp_{main fiber}}{n. exp_{tot}}",4,1,5,8,1,9);
  h2_spread->GetXaxis()->SetTitle("column");
  h2_spread->GetYaxis()->SetTitle("row");
  h22_spread->GetXaxis()->SetTitle("column");
  h22_spread->GetYaxis()->SetTitle("row");
  
  //////////Comparison of events yields
  for(int row=1;row<=8;row++){
    for(int column=1;column<=4;column++){
      Int_t obs_number = h_0[column][row]->GetEntries();
      
      float      exp_number_1 = (h[1][column][row]->GetEntries())*(float(ntrigs[0])/(ntrigs[1]));
      float      exp_number_2 = (h[2][column][row]->GetEntries())*(float(ntrigs[0])/(ntrigs[2]));
      float      exp_number_3 = (h[3][column][row]->GetEntries())*(float(ntrigs[0])/(ntrigs[3]));
      
      Int_t exp_number =  exp_number_1+exp_number_2+exp_number_3;

      /*
      cout<<"------------------"<<endl;
      cout<<"------------------"<<endl;
      cout<<"------------------"<<endl;
      cout<<exp_number_1<<"  =  "<<h[1][column][row]->GetEntries()<<" * ("<<float(ntrigs[0])<<")/("<<ntrigs[1]<<")"<<endl;
      cout<<exp_number_2<<"  =  "<<h[2][column][row]->GetEntries()<<" * ("<<float(ntrigs[0])<<")/("<<ntrigs[2]<<")"<<endl;
      cout<<exp_number_3<<"  =  "<<h[3][column][row]->GetEntries()<<" * ("<<float(ntrigs[0])<<")/("<<ntrigs[3]<<")"<<endl;
      cout<<exp_number<< " vs "<<obs_number<<endl;
      */
      h_spread->Fill((obs_number-exp_number)/float(obs_number));
      h3_spread->Fill(column,row,(obs_number-exp_number)/float(obs_number));
      h2_spread->SetBinContent(column,row,(obs_number-exp_number)/float(obs_number));
      h22_spread->SetBinContent(column,row,(exp_number-exp_number_2)/float(exp_number));

      /*
      cout<<"prima del rescaling"<<endl;

      cout<<"h[0][column][row]->GetEntries()        "<<h[0][column][row]->GetEntries()<<endl;
      cout<<"h[0][column][row]->GetIntegral()       "<<h[0][column][row]->GetIntegral()<<endl;
      cout<<"h[0][column][row]->Integral()          "<<h[0][column][row]->Integral()<<endl;

      cout<<"h[2][column][row]->GetEntries()        "<<h[2][column][row]->GetEntries()<<endl;
      cout<<"h[2][column][row]->GetIntegral()       "<<h[2][column][row]->GetIntegral()<<endl;
      cout<<"h[2][column][row]->Integral()          "<<h[2][column][row]->Integral()<<endl;
      */
      
      h[0][column][row]->Scale(double(exp_number_1)/(h[0][column][row]->Integral()));
      h[1][column][row]->Scale(double(exp_number_1)/(h[1][column][row]->Integral()));
      h[2][column][row]->Scale(double(exp_number_2)/(h[2][column][row]->Integral()));
      h[3][column][row]->Scale(double(exp_number_3)/(h[3][column][row]->Integral()));
      
      /*
      cout<<"prima del rescaling"<<endl;
      
      cout<<"h[0][column][row]->GetEntries()        "<<h[0][column][row]->GetEntries()<<endl;
      cout<<"h[0][column][row]->GetIntegral()       "<<h[0][column][row]->GetIntegral()<<endl;
      cout<<"h[0][column][row]->Integral()          "<<h[0][column][row]->Integral()<<endl;
      
      cout<<"h[2][column][row]->GetEntries()        "<<h[2][column][row]->GetEntries()<<endl;
      cout<<"h[2][column][row]->GetIntegral()       "<<h[2][column][row]->GetIntegral()<<endl;
      cout<<"h[2][column][row]->Integral()          "<<h[2][column][row]->Integral()<<endl;
      */
    
      for(int ii=2;ii<=3;ii++){
	h[0][column][row]->Add(h[ii][column][row]);
      }


      /*
      cout<<"------------------"<<endl;
      cout<<" test non normalization "<<endl;
      cout<<"h_0[column][row]->GetEntries()        "<<h_0[column][row]->GetEntries()<<endl;
      cout<<"h_0[column][row]->GetIntegral()       "<<h_0[column][row]->GetIntegral()<<endl;
      cout<<"h_0[column][row]->Integral()          "<<h_0[column][row]->Integral()<<endl;
      //h_0[column][row]->Scale(exp_number); 
      cout<<"h[0][column][row]->GetEntries()        "<<h[0][column][row]->GetEntries()<<endl;
      cout<<"h[0][column][row]->GetIntegral()       "<<h[0][column][row]->GetIntegral()<<endl;
      cout<<"h[0][column][row]->Integral()          "<<h[0][column][row]->Integral()<<endl;
      */

      
    }
  }
  
  TCanvas *c_spread = new TCanvas("c_spread","c_spread");
  c_spread->Divide(2,1);
  /*
  c_spread->cd(1);
  h_spread->DrawNormalized();
  c_spread->cd(2);
  h3_spread->Draw("colz");
  */
  c_spread->cd(1);
  h2_spread->Draw("colz");
  c_spread->cd(2);
  h22_spread->Draw("colz");

 
  
  /*
    int canvasID=-9;
    for(int row=1;row<=8;row++){
    for(int column=1;column<=4;column++){
    canvasID=column+4*(8-row);
    pmts->cd(canvasID);
    h[4][column][row]->SetLineColor(8);
      h[4][column][row]->DrawNormalized();
      h[3][column][row]->SetLineColor(4);
      h[3][column][row]->DrawNormalized("same");
      h[2][column][row]->SetLineColor(2);
      h[2][column][row]->DrawNormalized("same");
      h[1][column][row]->SetLineColor(1);
      h[1][column][row]->DrawNormalized("same");
    }
  }
  */
  
  int canvasID=-9;
  for(int row=1;row<=8;row++){
    for(int column=1;column<=4;column++){
      canvasID=column+4*(8-row);
      pmts->cd(canvasID);
      h_0[column][row]->SetLineColor(8);
      h[0][column][row]->SetLineColor(kMagenta+3);
      h[1][column][row]->SetLineColor(1);
      h[2][column][row]->SetLineColor(kOrange+7);
      h[3][column][row]->SetLineColor(4);
      
      h_0[column][row]->Draw();
      h[0][column][row]->Draw("same");
      h[1][column][row]->Draw("same");
      h[2][column][row]->Draw("same");
      h[3][column][row]->Draw("same");
      if(row<=4){
	pmts_down->cd(-4*row+column+16);
	h_0[column][row]->Draw();
	h[0][column][row]->Draw("same");
	h[1][column][row]->Draw("same");
	h[2][column][row]->Draw("same");
	h[3][column][row]->Draw("same");
	gPad->SetLogy();
      }
      else{
	pmts_up->cd(-4*(row-4)+column+16);
	h_0[column][row]->Draw();
	h[0][column][row]->Draw("same");
	h[1][column][row]->Draw("same");
	h[2][column][row]->Draw("same");
	h[3][column][row]->Draw("same");
	gPad->SetLogy();
      }
    }
  }



  TFile *f_out = new TFile("output_comparison.root","recreate");
  f_out->cd();
  c_spread->Write();
  pmts_up->Write();
  pmts_down->Write();
  pmts->Write();
  f_out->mkdir("histograms");
  f_out->cd("histograms");
  for(int row=1;row<=8;row++){
    for(int column=1;column<=4;column++){
      h_0[column][row]->Write();
      for(int ii=0;ii<=4;ii++){
	h[ii][column][row]->Write();
      }
    }
  }
  f_out->Close();
  delete f_out;
    

  
  
}

Fit_results_basic basic_fit_data(TString input_basepath, TString input_basefilename, int fit_model_d, int fit_model_r, int slotID, int pmt_pos, int column_number, int row_number){
  Fit_results_basic Results;
  if(fit_model_d>4||fit_model_r>4){cout<<"invalid model specified, execution stopped"<<endl; return Results;}

  
  if (input_basepath=="") input_basepath="Dati_PD/3.column_data/";
  
  //gROOT->ProcessLine(".L Fit_testbench.cpp"); loop_fit_PD_column(1); > ee.log
  /////////////////////////////////////////////////////////////////
  ////// HERE Starts the fit part /////////////////////////////////
  /////////////////////////////////////////////////////////////////
  bool raw_spectra = true;
  bool change_model = true;
  gROOT->ProcessLine(".x myRooPdfs/RooExpGauss.cxx+") ;
  gROOT->ProcessLine(".x myRooPdfs/RooAsymGauss.cxx+") ;

  TString input_filename="";
  TString input_column=Form("_col_%i",column_number);
  input_filename=input_basepath+input_basefilename+"_data_histos";
  TFile *f_input = new TFile(input_filename+".root"); //FOR PD data
 
  //TFile *f_input = new TFile(Form("KEK_data_histos_slot_8_col_%i_00049_allchs.root",column_number));
  Int_t pixelID=column_number+64*(row_number-1);

  
  double my_low_x=-0.7;//-0.5;//1;
  double my_up_x=2.2;//0.8;
  double time_shift=0.0;
  if(raw_spectra){
    TString slot_ID= "";
    TString slot_ID_undsc= "";
    TString slot_ID_bs= "";
    int shift_position=0; //when looking at KEK data to shift the position of the PMT while keeping the PD column numbering scheme
    if(slotID>0) {
      slot_ID = Form("slot_%i",slotID);
      slot_ID_undsc = "_";
      slot_ID_bs= "/";
      shift_position=1;
    }
    TString histo_name_base_slot=slot_ID+slot_ID_bs;
    TString histo_name_base=Form("column_%i/histos_%i_%i/",
				 column_number+shift_position*(pmt_pos-1)*4,
				 column_number+shift_position*(pmt_pos-1)*4,
				 row_number);
    TString histo_name_prefix_time = "time_";
    TString histo_name_prefix_amp = "amp_";
    TString histo_name_slot=slot_ID+slot_ID_undsc;
    TString histo_name_colrow= Form("col_%i_row_%i",
				    column_number+shift_position*(pmt_pos-1)*4,
				    row_number);
    TString histoname_time =histo_name_base_slot+histo_name_base+histo_name_prefix_time+histo_name_slot+histo_name_colrow;
    TString histoname_amp  =histo_name_base_slot+histo_name_base+histo_name_prefix_amp+histo_name_slot+histo_name_colrow;
    cout<<histoname_time<<endl;
    TH1D* h_time = f_input->Get(histoname_time);
    
    Float_t max_bin = h_time->GetMaximumBin();
    TAxis *xaxis = h_time->GetXaxis(); 
    Double_t max_pos = xaxis->GetBinCenter(max_bin);
    time_shift=max_pos;
    my_low_x=max_pos-0.7;//-0.5;//1;
    my_up_x=max_pos+1.5;//70;//0.8;
    
  }
  else{
    TH1D* h_time = f_input->Get(Form("column_%i/histos_%i_%i/h_time",column_number,column_number,row_number));
  }
  
  int column_number_MC=column_number+(pmt_pos-1)*4;
  /*
  TH1D* h_MC_tot = f_input->Get(Form("column_%i/histos_%i_%i/h_MC_tot",column_number_MC,column_number_MC,row_number));

  TH1D* h_MC_tot_ring = f_input->Get(Form("column_%i/histos_%i_%i/h_MC_ring_tot",column_number_MC,column_number_MC,row_number));
  
    TCanvas *can0 = f_input->Get(Form("column_%i/histos_%i_%i/c",column_number,column_number,row_number));
    can0->Draw();
  */

  
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
  bool add_background_component=true;
  int  amplitude_cut = -40;
  
  bool suppress_negligible_first_peak=false;
  bool do_simultaneous_fit=false;
  bool add_third_signal_pos0=false;
  //if(row_number>4){add_third_signal_pos0=true;}
  bool simulate_CB_tail=false;
  
  
  if(!print_prefit_info){MN_output_print_level_prefit=-1;}else{MN_output_print_level_prefit=MN_output_print_level;}
  bool draw_results;
  if(_draw_results=="draw"){draw_results=true;}else if(_draw_results=="blind"){draw_results=false;}else{draw_results=false;}


  
 
  /*
    TAxis *xaxis_MC = h_MC_tot->GetXaxis();
    xaxis_MC->SetRange(0,(TMath::Abs(my_low_x)/(my_up_x-my_low_x))*h_MC_tot->GetNbinsX()-1);
    Float_t max_bin_MC = h_MC_tot->GetMaximumBin();
    Double_t max_first_pos_MC = xaxis_MC->GetBinCenter(max_bin_MC);
    //  cout<<"ciao "<<max_first_pos_MC<<endl;
  */
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
   up_alpha_0=0.8;
  
  low_beta_0=0.5;
   up_beta_0=1.0;

  
   starting_delta_H_0=0.5;
   starting_delta_T_0=0.3;
   starting_sigma_L_0=0.08100;
   starting_sigma_H_0=0.08100;
   starting_sigma_T_0=0.108100;


   
   starting_alpha_0=0.7125;
   starting_beta_0=0.9;
   
   starting_mean_H_0=time_shift;//starting_mean_L_0+starting_delta_H_0;
   low_mean_H_0=starting_mean_H_0-0.175;//315;
   up_mean_H_0=starting_mean_H_0+0.175;//0.315;

   
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
    gROOT->ProcessLine(".L Fit_testbench.cpp"); loop_fit_column("","file","",3,4); > a.log
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
    //basic_fit_data(TString input_basepath, TString input_basefilename, int fit_model_d, int fit_model_r, int slotID, int pmt_pos, int column_number, int row_number)
    my_Results = basic_fit_data("",input_basefilename,2,fit_model_ID,-1,3,column_number,h);   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
