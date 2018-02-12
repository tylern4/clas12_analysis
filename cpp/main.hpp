/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD
#include <TFile.h>
#include <TLorentzVector.h>
#include "TH2.h"
#include <TFile.h>
#include <fstream>
#include "TF1.h"
#include "TChain.h"
#include <vector>

// PDG particle masses in GeV/c2
static const double MASS_P = 0.93827203;
static const double MASS_N = 0.93956556;
static const double MASS_E = 0.000511;
static const double MASS_PIP = 0.13957018;
static const double MASS_PIM = 0.13957018;
static const double MASS_PI0 = 0.1349766;
static const double MASS_KP = 0.493677;
static const double MASS_KM = 0.493677;
static const double MASS_G = 0.0;
static const double MASS_OMEGA = 0.78265;
static const double CLAS12_E = 10.7;
const double c_special_units = 29.9792458;
std::vector<int> *REC_Particle_pid;
std::vector<float> *REC_Particle_px;
std::vector<float> *REC_Particle_py;
std::vector<float> *REC_Particle_pz;
std::vector<float> *REC_Particle_vx;
std::vector<float> *REC_Particle_vy;
std::vector<float> *REC_Particle_vz;
std::vector<int> *REC_Particle_charge;
std::vector<float> *REC_Particle_beta;
std::vector<float> *REC_Particle_chi2pid;
std::vector<int> *REC_Particle_status;

std::vector<float> *FTOF_hits_pathLength;
std::vector<float> *FTOF_hits_time;

TH1D *momentum = new TH1D("mom", "mom", 500, 0, 10);
TH1D *W_hist = new TH1D("W", "W", 500, 0, 5);
TH1D *Q2_hist = new TH1D("Q2", "Q2", 500, 0, 10);
TH2D *W_vs_q2 = new TH2D("W_vs_q2", "W_vs_q2", 500, 0, 5, 500, 0, 8);
TH2D *mom_vs_beta = new TH2D("mom_vs_beta", "mom_vs_beta", 500, 0, 5, 500, 0.0, 1.2);
TH2D *mom_vs_beta_0 = new TH2D("mom_vs_beta_0", "mom_vs_beta_0", 500, 0, 5, 500, 0.0, 1.2);
TH2D *mom_vs_beta_pos = new TH2D("mom_vs_beta_pos", "mom_vs_beta_pos", 500, 0, 5, 500, 0.0, 1.2);
TH2D *mom_vs_beta_neg = new TH2D("mom_vs_beta_neg", "mom_vs_beta_neg", 500, 0, 5, 500, 0.0, 1.2);
TH2D *mom_vs_beta_proton =
    new TH2D("mom_vs_beta_proton", "mom_vs_beta_proton", 500, 0, 5, 500, 0.0, 1.2);
TH2D *mom_vs_beta_pion = new TH2D("mom_vs_beta_pion", "mom_vs_beta_pion", 500, 0, 5, 500, 0.0, 1.2);
TH2D *mom_vs_beta_electron =
    new TH2D("mom_vs_beta_electron", "mom_vs_beta_electron", 500, 0, 5, 500, 0.0, 1.2);

TH2D *deltat_proton =
    new TH2D("deltat_proton", "#Deltat assuming mass of proton", 500, -1.0, 10.0, 500, -10.0, 10.0);
TH2D *deltat_pion =
    new TH2D("deltat_pion", "#Deltat assuming mass of proton", 500, -1.0, 10.0, 500, -10.0, 10.0);

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

double vertex_time(double sc_time, double sc_pathlength, double relatavistic_beta) {
  return sc_time - sc_pathlength / (relatavistic_beta * c_special_units);
}

double delta_t(double electron_vertex_time, double mass, double momentum, double sc_t,
               double sc_r) {
  double relatavistic_beta = 1.0 / sqrt(1.0 + (mass / momentum) * (mass / momentum));
  return electron_vertex_time - vertex_time(sc_t, sc_r, relatavistic_beta);
}

void test(char *fin, char *fout) {
  double energy = CLAS12_E;
  if (getenv("CLAS12_E") != NULL)
    energy = atof(getenv("CLAS12_E"));

  TFile *out = new TFile(fout, "RECREATE");
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
  chain.SetBranchAddress("FTOF_hits_pathLength", &FTOF_hits_pathLength);
  chain.SetBranchAddress("FTOF_hits_time", &FTOF_hits_time);
  int num_of_events = (int)chain.GetEntries();
  int total = 0;
  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain.GetEntry(current_event);
    std::cout << "REC\t" << REC_Particle_pid->size() << std::endl;
    std::cout << "FTOF\t" << FTOF_hits_time->size() << std::endl;
    if (REC_Particle_pid->size() == 0 || FTOF_hits_time->size() == 0 ||
        FTOF_hits_time->size() < REC_Particle_pid->size())
      continue;

    double electron_vertex = vertex_time(FTOF_hits_time->at(0), FTOF_hits_pathLength->at(0), 1.0);
    double per = ((double)current_event / (double)num_of_events);
    //    std::cerr << "\t\t" << std::floor((100 * (double)current_event / (double)num_of_events))
    //    << "%\r\r" << std::flush;
    for (int i = 1; i < REC_Particle_pid->size(); i++) {
      double px = REC_Particle_px->at(i) * REC_Particle_px->at(i);
      double py = REC_Particle_py->at(i) * REC_Particle_py->at(i);
      double pz = REC_Particle_pz->at(i) * REC_Particle_pz->at(i);

      P = TMath::Sqrt(px + py + pz);
      if (REC_Particle_pid->at(i) == 0)
        mom_vs_beta_0->Fill(P, REC_Particle_beta->at(i));
      if (REC_Particle_beta->at(i) != 0) {
        momentum->Fill(P);
        mom_vs_beta->Fill(P, REC_Particle_beta->at(i));
        double dt_pion = delta_t(electron_vertex, MASS_PIP, P, FTOF_hits_time->at(i),
                                 FTOF_hits_pathLength->at(i));
        double dt_proton =
            delta_t(electron_vertex, MASS_P, P, FTOF_hits_time->at(i), FTOF_hits_pathLength->at(i));
        deltat_pion->Fill(P, dt_pion);
        deltat_proton->Fill(P, dt_proton);
        if (REC_Particle_charge->at(i) > 0) {
          mom_vs_beta_pos->Fill(P, REC_Particle_beta->at(i));
        } else if (REC_Particle_charge->at(i) < 0) {
          mom_vs_beta_neg->Fill(P, REC_Particle_beta->at(i));
        }

        if (REC_Particle_pid->at(i) == 2212)
          mom_vs_beta_proton->Fill(P, REC_Particle_beta->at(i));
        if (abs(REC_Particle_pid->at(i)) == 211)
          mom_vs_beta_pion->Fill(P, REC_Particle_beta->at(i));
        if (abs(REC_Particle_pid->at(i)) == 11)
          mom_vs_beta_electron->Fill(P, REC_Particle_beta->at(i));
      }

      total++;
      if (REC_Particle_pid->at(0) != 11)
        continue;
      // Setup scattered electron 4 vector
      TVector3 e_mu_prime_3;
      TLorentzVector e_mu_prime;
      TLorentzVector e_mu(0.0, 0.0, energy, energy);
      e_mu_prime_3.SetXYZ(REC_Particle_px->at(0), REC_Particle_py->at(0), REC_Particle_pz->at(0));
      e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
      double W = W_calc(e_mu, e_mu_prime);
      double Q2 = Q2_calc(e_mu, e_mu_prime);

      W_hist->Fill(W);
      Q2_hist->Fill(Q2);
      W_vs_q2->Fill(W, Q2);
    }
  }
  out->cd();
  momentum->Write();
  W_hist->Write();
  Q2_hist->Write();
  W_vs_q2->Write();
  mom_vs_beta->Write();
  mom_vs_beta_pos->Write();
  mom_vs_beta_neg->Write();
  mom_vs_beta_proton->Write();
  mom_vs_beta_pion->Write();
  mom_vs_beta_electron->Write();
  mom_vs_beta_0->Write();
  deltat_pion->Write();
  deltat_proton->Write();
  out->Close();
  chain.Reset();
  std::cerr << "\n" << total << std::endl;
}
#endif
