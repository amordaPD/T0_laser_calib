#include <iostream>
using namespace std;

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>

void make_flat_ntuple(int pos, double amplitude_cut=10, float xmin=10., float xmax=30, TString input_path,TString rfile){
  TString out_path="flat_ntuples/";
  double Amplitude (99);
  double Time (-999);
  double weight(-9);
  int Channel (-9);
  int NTimes (-9);
  
  TFile *f_data = new TFile(out_path+rfile+"_out.root","recreate");
  TTree *tree_input = new TTree("tree_input","tree_input");
  tree_input->Branch("Amplitude",&Amplitude,"Amplitude/D");
  tree_input->Branch("Time",&Time,"Time/D");
  tree_input->Branch("Channel",&Channel,"Channel/I");
  tree_input->Branch("NTimes",&NTimes,"NTimes/I");
  tree_input->Branch("weight",&weight,"weight/D");
  


  
  const int maxchans=16;
  const int maxtimes=5;
  int nchannels=0;
  int channels[maxchans];
  int ntimes[maxchans];
  float times[maxchans][maxtimes];
  float amplitudes[maxchans][maxtimes];

  TFile* ifile = new TFile(input_path+"/"+rfile+".root");
  TTree* tree;
  ifile->GetObject("times", tree);
  tree->SetBranchAddress("nch", &(nchannels));
  tree->SetBranchAddress("channels", &(channels[0]));
  tree->SetBranchAddress("ntimes", &(ntimes[0]));
  tree->SetBranchAddress("times", &(times[0][0]));
  tree->SetBranchAddress("amplitudes", &(amplitudes[0][0]));

  TH1F* htimes_0 = new TH1F(Form("fiber0-%d-0",pos), Form("fiber0-%d-0",pos), 200, xmin, xmax);
  TH1F* htimes_1 = new TH1F(Form("fiber0-%d-1",pos), Form("fiber0-%d-1",pos), 200, xmin, xmax);
  TH1F* htimes_2 = new TH1F(Form("fiber0-%d-2",pos), Form("fiber0-%d-2",pos), 200, xmin, xmax);
  TH1F* htimes_3 = new TH1F(Form("fiber0-%d-3",pos), Form("fiber0-%d-3",pos), 200, xmin, xmax);
  TH1F* htimes_4 = new TH1F(Form("fiber0-%d-4",pos), Form("fiber0-%d-4",pos), 200, xmin, xmax);
  TH1F* htimes_5 = new TH1F(Form("fiber0-%d-5",pos), Form("fiber0-%d-5",pos), 200, xmin, xmax);
  TH1F* htimes_6 = new TH1F(Form("fiber0-%d-6",pos), Form("fiber0-%d-6",pos), 200, xmin, xmax);
  TH1F* htimes_7 = new TH1F(Form("fiber0-%d-7",pos), Form("fiber0-%d-7",pos), 200, xmin, xmax);
  TH1F* htimes_8 = new TH1F(Form("fiber0-%d-8",pos), Form("fiber0-%d-8",pos), 200, xmin, xmax);
  TH1F* htimes_9 = new TH1F(Form("fiber0-%d-9",pos), Form("fiber0-%d-9",pos), 200, xmin, xmax);
  TH1F* htimes_10 = new TH1F(Form("fiber0-%d-10",pos), Form("fiber0-%d-10",pos), 200, xmin, xmax);
  TH1F* htimes_11 = new TH1F(Form("fiber0-%d-11",pos), Form("fiber0-%d-11",pos), 200, xmin, xmax);
  TH1F* htimes_12 = new TH1F(Form("fiber0-%d-12",pos), Form("fiber0-%d-12",pos), 200, xmin, xmax);
  TH1F* htimes_13 = new TH1F(Form("fiber0-%d-13",pos), Form("fiber0-%d-13",pos), 200, xmin, xmax);
  TH1F* htimes_14 = new TH1F(Form("fiber0-%d-14",pos), Form("fiber0-%d-14",pos), 200, xmin, xmax);
  TH1F* htimes_15 = new TH1F(Form("fiber0-%d-15",pos), Form("fiber0-%d-15",pos), 200, xmin, xmax);

  for (int k=0; k<tree->GetEntries(); ++k) {
    tree->GetEntry(k);
    for (int l=0; l<nchannels; ++l) {
      for (int j=0; j<ntimes[l]; ++j) {
	if(amplitudes[l][j]<-amplitude_cut){
	  if(channels[l]==0) htimes_0->Fill(times[l][j]);
	  if(channels[l]==1) htimes_1->Fill(times[l][j]);
	  if(channels[l]==2) htimes_2->Fill(times[l][j]);
	  if(channels[l]==3) htimes_3->Fill(times[l][j]);
	  if(channels[l]==4) htimes_4->Fill(times[l][j]);
	  if(channels[l]==5) htimes_5->Fill(times[l][j]);
	  if(channels[l]==6) htimes_6->Fill(times[l][j]);
	  if(channels[l]==7) htimes_7->Fill(times[l][j]);
	  if(channels[l]==8) htimes_8->Fill(times[l][j]);
	  if(channels[l]==9) htimes_9->Fill(times[l][j]);
	  if(channels[l]==10) htimes_10->Fill(times[l][j]);
	  if(channels[l]==11) htimes_11->Fill(times[l][j]);
	  if(channels[l]==12) htimes_12->Fill(times[l][j]);
	  if(channels[l]==13) htimes_13->Fill(times[l][j]);
	  if(channels[l]==14) htimes_14->Fill(times[l][j]);
	  if(channels[l]==15) htimes_15->Fill(times[l][j]);
	}
	Amplitude=amplitudes[l][j];
	Time=times[l][j];
	Channel=channels[l];
	NTimes=ntimes[l];
	if(TMath::Abs(Time-23.2)<0.3) {weight=1.0;}else{weight=0.3/TMath::Abs(Time-23.2);}
	tree_input->Fill();
      }
    }
  }
  

  f_data->cd();
  tree_input->Write();
  htimes_0->Write();
  htimes_1->Write();
  htimes_2->Write();
  htimes_3->Write();
  htimes_4->Write();
  htimes_5->Write();
  htimes_6->Write();
  htimes_7->Write();
  htimes_8->Write();
  htimes_9->Write();
  htimes_10->Write();
  htimes_11->Write();
  htimes_12->Write();
  htimes_13->Write();
  htimes_14->Write();
  htimes_15->Write();
  f_data->Close();
  delete f_data;


  
}



TH1F* plottimes(const char* rfile, int ch, int pos, double amplitude_cut=10, float xmin=10., float xmax=30) {
  const int maxchans=16;
  const int maxtimes=5;
  int nchannels=0;
  int channels[maxchans];
  int ntimes[maxchans];
  float times[maxchans][maxtimes];
  float amplitudes[maxchans][maxtimes];

  TFile* ifile = new TFile(rfile);
  TTree* tree;
  ifile->GetObject("times", tree);
  tree->SetBranchAddress("nch", &(nchannels));
  tree->SetBranchAddress("channels", &(channels[0]));
  tree->SetBranchAddress("ntimes", &(ntimes[0]));
  tree->SetBranchAddress("times", &(times[0][0]));
  tree->SetBranchAddress("amplitudes", &(amplitudes[0][0]));

  TH1F* htimes = new TH1F(Form("fiber0-%d-%d",pos,ch), Form("fiber0-%d-%d",pos,ch), 200, xmin, xmax);

  for (int k=0; k<tree->GetEntries(); ++k) {
    tree->GetEntry(k);
    for (int l=0; l<nchannels; ++l) {
      if (channels[l]==ch) {
        for (int j=0; j<ntimes[l]; ++j) {
          if(amplitudes[l][j]<-amplitude_cut) htimes->Fill(times[l][j]/*-14*/);
        }
      }
    }
  }

  TCanvas* c = new TCanvas("a", "a");
  htimes->Draw();
  return htimes;
}


TH1F* plotamp(const char* rfile, int ch, int pos, double amplitude_cut=10, float xmin=-200, float xmax=0) {
  const int maxchans=16;
  const int maxtimes=5;
  int nchannels=0;
  int channels[maxchans];
  int ntimes[maxchans];
  float times[maxchans][maxtimes];
  float amplitudes[maxchans][maxtimes];

  TFile* ifile = new TFile(rfile);
  TTree* tree;
  ifile->GetObject("times", tree);
  tree->SetBranchAddress("nch", &(nchannels));
  tree->SetBranchAddress("channels", &(channels[0]));
  tree->SetBranchAddress("ntimes", &(ntimes[0]));
  tree->SetBranchAddress("times", &(times[0][0]));
  tree->SetBranchAddress("amplitudes", &(amplitudes[0][0]));

  TH1F* htimes = new TH1F(Form("amp-fiber0-%d-%d",pos,ch), Form("amp-fiber0-%d-%d",pos,ch), 200, xmin, xmax);

  for (int k=0; k<tree->GetEntries(); ++k) {
    tree->GetEntry(k);
    for (int l=0; l<nchannels; ++l) {
      if (channels[l]==ch) {
        for (int j=0; j<ntimes[l]; ++j) {
          if(amplitudes[l][j]<-amplitude_cut) htimes->Fill(amplitudes[l][j]);
        }
      }
    }
  }

  TCanvas* c = new TCanvas("a", "a");
  htimes->Draw();
  return htimes;
}
