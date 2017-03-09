#include <iostream>
using namespace std;

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>

void make_flat_ntuple(TString input_path,TString rfile){
  TString out_path="flat_ntuples/";
  double Amplitude (99);
  double Time (-999);
  int Channel (-9);
  TFile *f_data = new TFile(out_path+rfile+"_out.root","recreate");
  TTree *tree_input = new TTree("tree_input","tree_input");
  tree_input->Branch("Amplitude",&Amplitude,"Amplitude/D");
  tree_input->Branch("Time",&Time,"Time/D");
  tree_input->Branch("Channel",&Channel,"Channel/I");
  


  
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

  //TH1F* htimes = new TH1F(Form("fiber0-%d-%d",pos,ch), Form("fiber0-%d-%d",pos,ch), 200, xmin, xmax);

  for (int k=0; k<tree->GetEntries(); ++k) {
    tree->GetEntry(k);
    for (int l=0; l<nchannels; ++l) {
      for (int j=0; j<ntimes[l]; ++j) {
	//if(amplitudes[l][j]<-amplitude_cut) htimes->Fill(times[l][j]/*-14*/);
	Amplitude=amplitudes[l][j];
	Time=times[l][j];
	Channel=channels[l];
	tree_input->Fill();
      }
    }
  }
  

  f_data->cd();
  tree_input->Write();
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
