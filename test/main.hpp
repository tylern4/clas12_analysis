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
  std::vector<int>     *REC_Particle_pid;
  std::vector<float>   *REC_Particle_px;
  std::vector<float>   *REC_Particle_py;
  std::vector<float>   *REC_Particle_pz;
  std::vector<float>   *REC_Particle_vx;
  std::vector<float>   *REC_Particle_vy;
  std::vector<float>   *REC_Particle_vz;
  std::vector<int>     *REC_Particle_charge;
  std::vector<float>   *REC_Particle_beta;
  std::vector<float>   *REC_Particle_chi2pid;
  std::vector<int>     *REC_Particle_status;
  TH1D *momentum = new TH1D("mom","mom",500,0,10);
  TH1D *W_hist = new TH1D("W","W",250,0,5);
  TH1D *Q2_hist = new TH1D("Q2","Q2",250,0,10);
  TH2D *W_vs_q2 = new TH2D("W_vs_q2","W_vs_q2",250,0,5,250,0,10);
  TH2D *mom_vs_beta = new TH2D("mom_vs_beta","mom_vs_beta",500,0,5,500,-2.5,2.5);


// Calcuating Q^2
// q^mu^2 = (e^mu - e^mu')^2 = -Q^2
double Q2_calc(TLorentzVector e_mu, TLorentzVector e_mu_prime) {
  TLorentzVector q_mu = (e_mu - e_mu_prime);
  return -q_mu.Mag2();
}
//	Calcualting W
//	Gotten from s channel [(gamma - P)^2 == s == w^2]
//	Sqrt√[M_p^2 - Q^2 + 2 M_p gamma]
double W_calc(TLorentzVector e_mu, TLorentzVector e_mu_prime) {
  TLorentzVector q_mu = (e_mu - e_mu_prime);
  TVector3 p_mu_3(0, 0, 0);
  TLorentzVector p_mu;
  p_mu.SetVectM(p_mu_3, MASS_P);
  return (p_mu + q_mu).Mag();
}

void test(char *fin) {
  TFile *out = new TFile("out.root", "RECREATE");
  double P;
  bool electron_cuts;
  // Load chain from branch h10
  TChain chain("clas12");
  chain.Add(fin);

  chain.SetBranchAddress("REC_Particle_pid", &REC_Particle_pid);
  chain.SetBranchAddress("REC_Particle_px", &REC_Particle_px);
  chain.SetBranchAddress("REC_Particle_py", &REC_Particle_py);
  chain.SetBranchAddress("REC_Particle_pz", &REC_Particle_pz);
  chain.SetBranchAddress("REC_Particle_vx", &REC_Particle_vx);
  chain.SetBranchAddress("REC_Particle_vy", &REC_Particle_vy);
  chain.SetBranchAddress("REC_Particle_vz", &REC_Particle_vz);
  chain.SetBranchAddress("REC_Particle_charge", &REC_Particle_charge);
  chain.SetBranchAddress("REC_Particle_beta", &REC_Particle_beta);
  chain.SetBranchAddress("REC_Particle_chi2pid", &REC_Particle_chi2pid);
  chain.SetBranchAddress("REC_Particle_status", &REC_Particle_status);
  int num_of_events = (int)chain.GetEntries();

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain.GetEntry(current_event);
    if (REC_Particle_pid->size() == 0 ) continue;
    std::cerr << REC_Particle_pid->at(0) << std::endl;
    //if (TMath::Abs(REC_Particle_pid->at(0)) != 11 ) continue;
    // Setup scattered electron 4 vector
    TVector3 e_mu_prime_3;
    TLorentzVector e_mu_prime;
    TLorentzVector e_mu(0.0, 0.0, 10.7, 7);
    e_mu_prime_3.SetXYZ(REC_Particle_px->at(0), REC_Particle_py->at(0), REC_Particle_pz->at(0));
    e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
    double W = W_calc(e_mu, e_mu_prime);
    double Q2 = Q2_calc(e_mu, e_mu_prime);

    W_hist->Fill(W);
    Q2_hist->Fill(Q2);
    W_vs_q2->Fill(W,Q2);

    for (int i = 0; i< REC_Particle_pid->size(); i++){
      double px = REC_Particle_px->at(i) * REC_Particle_px->at(i);
      double py = REC_Particle_py->at(i) * REC_Particle_py->at(i);
      double pz = REC_Particle_pz->at(i) * REC_Particle_pz->at(i);

      P = TMath::Sqrt(px + py + pz);
      momentum->Fill(P);
      mom_vs_beta->Fill(P,REC_Particle_beta->at(i));

    }
  }
  out->cd();
  momentum->Write();
  W_hist->Write();
  Q2_hist->Write();
  W_vs_q2->Write();
  mom_vs_beta->Write();
  out->Close();
  chain.Reset();
}
#endif
