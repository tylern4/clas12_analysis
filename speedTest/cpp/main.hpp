/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef CSV_MAKER_H_GUARD
#define CSV_MAKER_H_GUARD
#include <stdlib.h>
#include <stdio.h>
#include <TFileCollection.h>
#include <TFile.h>
#include "TMath.h"
#include <TGraph.h>
#include "THnSparse.h"
#include "TTree.h"
#include "TROOT.h"
#include <TLorentzVector.h>
#include <string.h>
#include <string>
#include <cstring>
#include "time.h"
#include <string>
#include <cstring>
#include "TTree.h"
#include "TROOT.h"
#include "TH2.h"
#include <TFile.h>
#include "TStyle.h"
#include "TCanvas.h"
#include <fstream>
#include "TF1.h"
#include "TLeaf.h"
#include "TChain.h"
#include "TSystem.h"
#include <vector>
#include <fstream>
static const double MASS_E = 0.000511;
static const double MASS_P = 0.93827203;

std::vector<float> *REC_Particle_px;
std::vector<float> *REC_Particle_py;
std::vector<float> *REC_Particle_pz;
std::vector<float> *REC_Particle_beta;

TH1D *momentum = new TH1D("mom", "mom", 500, 0, 10);
TH2D *mom_vs_beta = new TH2D("mom_vs_beta", "mom_vs_beta", 500, 0, 5, 500, -2.5, 2.5);

void test(char *fin, char *fout) {
  TFile *out = new TFile(fout, "RECREATE");
  double P;
  bool electron_cuts;
  // Load chain from branch h10
  TChain chain("clas12");
  chain.Add(fin);

  chain.SetBranchAddress("REC_Particle_px", &REC_Particle_px);
  chain.SetBranchAddress("REC_Particle_py", &REC_Particle_py);
  chain.SetBranchAddress("REC_Particle_pz", &REC_Particle_pz);
  chain.SetBranchAddress("REC_Particle_beta", &REC_Particle_beta);
  int num_of_events = (int)chain.GetEntries();
  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain.GetEntry(current_event);
    if (REC_Particle_px->size() == 0) continue;

    for (int i = 1; i < REC_Particle_px->size(); i++) {
      double px = REC_Particle_px->at(i) * REC_Particle_px->at(i);
      double py = REC_Particle_py->at(i) * REC_Particle_py->at(i);
      double pz = REC_Particle_pz->at(i) * REC_Particle_pz->at(i);

      P = TMath::Sqrt(px + py + pz);
      momentum->Fill(P);
      if (REC_Particle_beta->at(i) != 0) mom_vs_beta->Fill(P, REC_Particle_beta->at(i));
    }
  }
  out->cd();
  momentum->Write();
  mom_vs_beta->Write();
  out->Close();
  chain.Reset();
}
#endif
