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
  TLorentzVector *e_mu = new TLorentzVector(0.0, 0.0, energy, energy);

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
  double per = 0;
  int index = 0;
  int num_pip = 0;
  TLorentzVector e_mu_prime;
  bool good_e = true;

  Histogram *hist = new Histogram();

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    if (pid->size() == 0) continue;

    per = ((double)current_event / (double)num_of_events);
    if (current_event % 1000 == 0) std::cerr << "\t\t" << std::floor(100 * per) << "%\r\r" << std::flush;
    e_mu_prime.SetXYZM(px->at(0), py->at(0), pz->at(0), MASS_E);

    if (e_mu_prime.P() != 0) hist->Fill_EC(ec_tot_energy->at(0) / e_mu_prime.P(), e_mu_prime.P());

    Delta_T *ftof_dt = new Delta_T(sc_ftof_time->at(0), sc_ftof_path->at(0));
    Delta_T *ctof_dt = new Delta_T(sc_ctof_time->at(0), sc_ctof_path->at(0));

    for (int part = 1; part < pid->size(); part++) {
      double P = TMath::Sqrt((px->at(part) * px->at(part) + py->at(part) * py->at(part) + pz->at(part) * pz->at(part)));
      if (beta->at(part) < 0.02 || P < 0.02) continue;
      ftof_dt->deltat(P, sc_ftof_time->at(part), sc_ftof_path->at(part));
      ctof_dt->deltat(P, sc_ctof_time->at(part), sc_ctof_path->at(part));
      hist->Fill_MomVsBeta(pid->at(part), charge->at(part), P, beta->at(part));
      hist->Fill_deltat(pid->at(part), charge->at(part), P, ftof_dt);
      hist->Fill_deltat(pid->at(part), charge->at(part), P, ctof_dt);
    }

    delete ftof_dt;
    delete ctof_dt;
    W = physics::W_calc(*e_mu, e_mu_prime);
    Q2 = physics::Q2_calc(*e_mu, e_mu_prime);
    hist->Fill_WvsQ2(W, Q2);
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
  // delete hist;
}
#endif
