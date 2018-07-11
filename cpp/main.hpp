/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD
#include <TFile.h>
#include <TLorentzVector.h>
#include <fstream>
#include <vector>
#include "TChain.h"
#include "colors.hpp"
#include "constants.hpp"
#include "deltat.hpp"
#include "filehandeler.hpp"
#include "histogram.hpp"
#include "physics.hpp"

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
  double tot_energy_ec = 0;
  int sc_d = 0;
  double W = 0;
  double Q2 = 0;
  double sf = 0;
  double P_x = 0;
  double P_y = 0;
  double P_z = 0;
  double per = 0;
  int index = 0;
  int num_pip = 0;
  TVector3 e_mu_prime_3;
  TLorentzVector e_mu_prime;
  bool good_e = true;

  Histogram *hist = new Histogram();

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);

    if (pid->size() == 0 || sc_time->size() == 0 || ec_pindex->size() == 0) continue;

    per = ((double)current_event / (double)num_of_events);
    std::cerr << "\t\t" << std::floor(100 * per) << "%\r\r" << std::flush;

    num_pip = 0;
    tot_energy_ec = 0;
    good_e = true;
    for (int j = 0; j < ec_pindex->size(); j++) {
      if (ec_pindex->size() == 0) continue;
      try {
        index = ec_pindex->at(j);
        if (pid->at(index) == ELECTRON && index == 0) {
          e_mu_prime_3.SetXYZ(px->at(index), py->at(index), pz->at(index));
          P = e_mu_prime_3.Mag();
          e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
          tot_energy_ec += ec_energy->at(j);
          good_e = true;
        }
      } catch (std::exception &e) {
        total++;
      }
    }
    sf = tot_energy_ec / e_mu_prime.P();
    if (sf != 0) hist->Fill_EC(sf, e_mu_prime.P());
    // if (!good_e) continue;
    good_e = true;
    for (int j = 0; j < sc_time->size(); j++) {
      if (sc_time->size() == 0) continue;
      try {
        Delta_T *dt = new Delta_T(sc_time->at(0), sc_path->at(0));
        index = sc_pindex->at(j);
        sc_d = sc_detector->at(j);
        // I think 12 is FTOF
        if (sc_d == 12) {
          P_x = px->at(index) * px->at(index);
          P_y = py->at(index) * py->at(index);
          P_z = pz->at(index) * pz->at(index);
          P = TMath::Sqrt(P_x + P_y + P_z);

          dt->deltat(P, sc_time->at(j), sc_path->at(j));

          if (index == 0) {
            hist->Fill_MomVsBeta_vertex(pid->at(index), charge->at(index), P, beta->at(index));
            hist->Fill_deltat_vertex(pid->at(index), charge->at(index), P, dt);
          } else {
            hist->Fill_MomVsBeta(pid->at(index), charge->at(index), P, beta->at(index));
            hist->Fill_deltat(pid->at(index), charge->at(index), P, dt);
          }
        }
        if (pid->at(sc_pindex->at(j)) == PIP && abs(dt->Get_dt_Pi()) < 0.5) num_pip++;
        if (pid->at(sc_pindex->at(j)) == ELECTRON && sc_detector->at(sc_pindex->at(j)) == 12) good_e = true;
        delete dt;
      } catch (std::exception &e) {
        total++;
      }
    }
    if (!good_e) continue;
    // && sf >= 0.07 && sf <= 0.26
    if (px->size() > 0) {
      e_mu_prime_3.SetXYZ(px->at(0), py->at(0), pz->at(0));
      e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
      if (e_mu_prime.P() > 1.5) {
        W = physics::W_calc(e_mu, e_mu_prime);
        Q2 = physics::Q2_calc(e_mu, e_mu_prime);
        hist->Fill_WvsQ2(W, Q2);
        if (num_pip == 1 && pid->size() == 2) hist->Fill_WvsQ2_singlePi(W, Q2);
      }
    }
  }

  out->cd();
  hist->Write_EC();
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
  std::cerr << "\nErrors: " << total << "\t" << std::endl;
}

void datahandeler2(char *fin) {
  double energy = CLAS12_E;
  if (getenv("CLAS12_E") != NULL) energy = atof(getenv("CLAS12_E"));
  TLorentzVector e_mu(0.0, 0.0, energy, energy);

  TFile *out = new TFile("out_reindex.root", "RECREATE");
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
    if (pid->size() == 0 || sc_time->size() == 0) continue;

    double per = ((double)current_event / (double)num_of_events);
    std::cerr << "\t\t" << std::floor(100 * per) << "%\r\r" << std::flush;

    int vertex_id = 0;
    for (int j = 0; j < sc_time->size(); j++) {
      int temp = sc_pindex->at(j);
      if (temp == 0) {
        vertex_id = j;
        continue;
      }
    }

    double tot_energy_ec = 0;
    for (int j = 0; j < ec_pindex->size(); j++) {
      if (ec_pindex->size() == 0) continue;
      int index = ec_pindex->at(j);
      if (pid->at(ec_pindex->at(j)) == ELECTRON && ec_pindex->at(j) == 0) {
        tot_energy_ec += ec_energy->at(j);
        double P_x = px->at(index) * px->at(index);
        double P_y = py->at(index) * py->at(index);
        double P_z = pz->at(index) * pz->at(index);
        P = TMath::Sqrt(P_x + P_y + P_z);
      }
    }
    hist->Fill_EC(tot_energy_ec, P);

    try {
      if (pid->at(vertex_id) == ELECTRON) {
        TVector3 e_mu_prime_3;
        TLorentzVector e_mu_prime;
        e_mu_prime_3.SetXYZ(px->at(vertex_id), py->at(vertex_id), pz->at(vertex_id));

        e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
        double W = physics::W_calc(e_mu, e_mu_prime);
        double Q2 = physics::Q2_calc(e_mu, e_mu_prime);
        hist->Fill_WvsQ2(W, Q2);
      }
    } catch (std::exception &e) {
      total++;
    }

    for (int j = 0; j < sc_time->size(); j++) {
      if (sc_time->size() == 0) continue;
      try {
        Delta_T *dt = new Delta_T(sc_time->at(vertex_id), sc_path->at(vertex_id));
        int index = sc_pindex->at(j);

        double P_x = px->at(index) * px->at(index);
        double P_y = py->at(index) * py->at(index);
        double P_z = pz->at(index) * pz->at(index);
        P = TMath::Sqrt(P_x + P_y + P_z);

        dt->deltat(P, sc_time->at(j), sc_path->at(j));

        if (index == 0) {
          hist->Fill_MomVsBeta_vertex(pid->at(index), charge->at(index), P, beta->at(index));
          hist->Fill_deltat_vertex(pid->at(index), charge->at(index), P, dt);
        } else {
          hist->Fill_MomVsBeta(pid->at(index), charge->at(index), P, beta->at(index));
          hist->Fill_deltat(pid->at(index), charge->at(index), P, dt);
        }
        delete dt;
      } catch (std::exception &e) {
        continue;

        total++;
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
      int e_index = 0;
      int pip_index = 1;
      for (int j = 0; j < sc_time->size(); j++) {
        int temp = sc_pindex->at(j);
        if (temp == 0) {
          e_index = j;
          continue;
        }
      }
      try {
        pip_index = (e_index == 0) ? 1 : 0;
        if (pid->at(e_index) == ELECTRON && (pid->at(pip_index) == PIP || pid->at(pip_index) == KP)) {
          TVector3 e_mu_prime_3;
          TLorentzVector e_mu_prime;
          e_mu_prime_3.SetXYZ(px->at(e_index), py->at(e_index), pz->at(e_index));
          e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
          double W = physics::W_calc(e_mu, e_mu_prime);
          double Q2 = physics::Q2_calc(e_mu, e_mu_prime);
          hist->Fill_WvsQ2(W, Q2);

          for (int index = 0; index < pid->size(); index++) {
            P_x = px->at(index) * px->at(index);
            P_y = py->at(index) * py->at(index);
            P_z = pz->at(index) * pz->at(index);
            P = TMath::Sqrt(P_x + P_y + P_z);
            if (index == 0) {
              hist->Fill_MomVsBeta_vertex(pid->at(index), charge->at(index), P, beta->at(index));
            } else {
              hist->Fill_MomVsBeta(pid->at(index), charge->at(index), P, beta->at(index));
            }
          }
        }
      } catch (std::exception &e) {
        continue;

        total++;
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
  std::cerr << "\nErrors:" << total << "\t" << std::endl;
}

#endif
