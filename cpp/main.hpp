/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD
#include <TFile.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <fstream>
#include "TChain.h"
#include <vector>
#include "colors.hpp"
#include "filehandeler.hpp"
#include "constants.hpp"
#include "physics.hpp"
#include "histogram.hpp"
#include "deltat.hpp"

void datahandeler(char *fin, char *fout) {
  double energy = CLAS12_E;
  if (getenv("CLAS12_E") != NULL) energy = atof(getenv("CLAS12_E"));
  TLorentzVector e_mu(0.0, 0.0, energy, energy);

  TFile *out = new TFile(fout, "RECREATE");
  double P;
  bool electron_cuts;
  // Load chain from branch h10
  TChain *chain = filehandeler::addFiles(fin);
  filehandeler::getBranches(chain);

  int num_of_events = (int)chain->GetEntries();
  int total = 0;

  Histogram *hist = new Histogram();

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    if (pid->size() == 0) continue;

    double per = ((double)current_event / (double)num_of_events);
    std::cerr << "\t\t" << std::floor(100 * per) << "%\r\r" << std::flush;

    for (int i = 0; i < pid->size(); i++) {
      if (pid->size() == 0) continue;
      total++;

      double P_x = px->at(i) * px->at(i);
      double P_y = py->at(i) * py->at(i);
      double P_z = pz->at(i) * pz->at(i);

      P = TMath::Sqrt(P_x + P_y + P_z);

      if (pid->at(0) != 11) continue;
      // Setup scattered electron 4 vector
      TVector3 e_mu_prime_3;
      TLorentzVector e_mu_prime;
      e_mu_prime_3.SetXYZ(px->at(0), py->at(0), pz->at(0));
      e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
      double W = physics::W_calc(e_mu, e_mu_prime);
      double Q2 = physics::Q2_calc(e_mu, e_mu_prime);
      hist->Fill_WvsQ2(W, Q2);
    }

    int vertex_id;
    for (int j = 0; j < sc_time->size(); j++) {
      int temp = pindex->at(j);
      if (temp == 0) {
        vertex_id = temp;
        continue;
      }
    }

    for (int j = 0; j < sc_time->size(); j++) {
      if (sc_time->size() == 0) continue;

      Delta_T *dt = new Delta_T(sc_time->at(vertex_id), sc_r->at(vertex_id));
      int index = pindex->at(j);
      if (beta->at(index) <= 0.05) continue;

      double P_x = px->at(index) * px->at(index);
      double P_y = py->at(index) * py->at(index);
      double P_z = pz->at(index) * pz->at(index);
      P = TMath::Sqrt(P_x + P_y + P_z);

      dt->deltat(P, sc_time->at(j), sc_r->at(j));

      if (index == 0) {
        hist->Fill_MomVsBeta_vertex(pid->at(index), charge->at(index), P, beta->at(index));
        hist->Fill_deltat_vertex(pid->at(index), charge->at(index), P, dt);
      } else {
        if (beta->at(index) >= 0.05) {
          hist->Fill_MomVsBeta(pid->at(index), charge->at(index), P, beta->at(index));
          hist->Fill_deltat(pid->at(index), charge->at(index), P, dt);
        }
      }
      delete dt;
    }
  }

  out->cd();
  TDirectory *wvsq2 = out->mkdir("wvsq2");
  wvsq2->cd();
  hist->Write_WvsQ2();

  TDirectory *mom_vs_beta = out->mkdir("mom_vs_beta");
  mom_vs_beta->cd();
  hist->Write_MomVsBeta();

  TDirectory *deltat_ftof = out->mkdir("deltat_ftof");
  deltat_ftof->cd();
  hist->Write_deltat();

  out->Close();
  chain->Reset();
  std::cerr << "\n" << total << "\t" << std::endl;
}

void SinglePi(char *fin, char *fout) {
  double energy = CLAS12_E;
  if (getenv("CLAS12_E") != NULL) energy = atof(getenv("CLAS12_E"));
  TLorentzVector e_mu(0.0, 0.0, energy, energy);
  TFile *out = new TFile(fout, "RECREATE");
  double P;
  bool electron_cuts;
  // Load chain from branch h10
  TChain *chain = filehandeler::addFiles(fin);
  filehandeler::getBranches(chain);

  int num_of_events = (int)chain->GetEntries();
  int total = 0;
  double W = 0;
  double Q2 = 0;
  double P_x = 0;
  double P_y = 0;
  double P_z = 0;
  double per = 0;
  TVector3 e_mu_prime_3;
  TLorentzVector e_mu_prime;

  Histogram *hist = new Histogram();

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    if (pid->size() != 2) continue;

    per = ((double)current_event / (double)num_of_events);
    std::cerr << "\t\t" << std::floor(100 * per) << "%\r\r" << std::flush;
    if (pid->size() == 2) {
      total++;
      if (pid->at(0) != ELECTRON) continue;
      if (pid->at(1) == PIP || pid->at(1) == KP) {
        // Setup scattered electron 4 vector
        e_mu_prime_3.SetXYZ(px->at(0), py->at(0), pz->at(0));
        e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
        W = physics::W_calc(e_mu, e_mu_prime);
        Q2 = physics::Q2_calc(e_mu, e_mu_prime);
        hist->Fill_WvsQ2(W, Q2);
      }
    }
    for (int index = 0; index < pid->size(); index++) {
      P_x = px->at(index) * px->at(index);
      P_y = py->at(index) * py->at(index);
      P_z = pz->at(index) * pz->at(index);
      P = TMath::Sqrt(P_x + P_y + P_z);
      if (index == 0) {
        hist->Fill_MomVsBeta_vertex(pid->at(index), charge->at(index), P, beta->at(index));
      } else {
        if (beta->at(index) >= 0.05) {
          hist->Fill_MomVsBeta(pid->at(index), charge->at(index), P, beta->at(index));
        }
      }
    }
  }

  out->cd();
  TDirectory *wvsq2 = out->mkdir("wvsq2");
  wvsq2->cd();
  hist->Write_WvsQ2();

  TDirectory *mom_vs_beta = out->mkdir("mom_vs_beta");
  mom_vs_beta->cd();
  hist->Write_MomVsBeta();

  out->Close();
  chain->Reset();
  std::cerr << "\n" << total << "\t" << std::endl;
}

void TwoPi(char *fin, char *fout) {
  double energy = CLAS12_E;
  if (getenv("CLAS12_E") != NULL) energy = atof(getenv("CLAS12_E"));
  TLorentzVector e_mu(0.0, 0.0, energy, energy);
  TFile *out = new TFile(fout, "RECREATE");
  double P;
  bool electron_cuts;
  // Load chain from branch h10
  TChain *chain = filehandeler::addFiles(fin);
  filehandeler::getBranches(chain);

  int num_of_events = (int)chain->GetEntries();
  int total = 0;
  double W = 0;
  double Q2 = 0;
  double P_x = 0;
  double P_y = 0;
  double P_z = 0;
  double per = 0;

  Histogram *hist = new Histogram();

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    if (pid->size() != 2) continue;

    per = ((double)current_event / (double)num_of_events);
    std::cerr << "\t\t" << std::floor(100 * per) << "%\r\r" << std::flush;

    for (int i = 0; i < pid->size(); i++) {
      if (pid->size() != 2) continue;
      total++;

      P_x = px->at(i) * px->at(i);
      P_y = py->at(i) * py->at(i);
      P_z = pz->at(i) * pz->at(i);

      P = TMath::Sqrt(P_x + P_y + P_z);

      if (pid->at(0) != ELECTRON) continue;
      if (pid->at(1) != PIP) continue;
      // Setup scattered electron 4 vector
      TVector3 e_mu_prime_3;
      TLorentzVector e_mu_prime;
      e_mu_prime_3.SetXYZ(px->at(0), py->at(0), pz->at(0));
      e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
      W = physics::W_calc(e_mu, e_mu_prime);
      Q2 = physics::Q2_calc(e_mu, e_mu_prime);
      hist->Fill_WvsQ2(W, Q2);
    }
  }

  out->cd();
  TDirectory *wvsq2 = out->mkdir("wvsq2");
  wvsq2->cd();
  hist->Write_WvsQ2();

  TDirectory *mom_vs_beta = out->mkdir("mom_vs_beta");
  mom_vs_beta->cd();
  hist->Write_MomVsBeta();

  out->Close();
  chain->Reset();
  std::cerr << "\n" << total << "\t" << std::endl;
}

#endif
