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
#include "filehandeler.hpp"
#include "constants.hpp"
#include "physics.hpp"
#include "histogram.hpp"

void test(char *fin, char *fout) {
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

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    if (pid->size() == 0) continue;

    double per = ((double)current_event / (double)num_of_events);
    std::cerr << "\t\t" << std::floor(100 * per) << "%\r\r" << std::flush;

    for (int i = 0; i < pid->size(); i++) {
      total++;
      double P_x = px->at(i) * px->at(i);
      double P_y = py->at(i) * py->at(i);
      double P_z = pz->at(i) * pz->at(i);

      P = TMath::Sqrt(P_x + P_y + P_z);
      if (i == 0 && beta->at(i) != 0) {
        mom_vs_beta_0th->Fill(P, beta->at(i));
        continue;
      }

      if (beta->at(i) != 0) {
        momentum->Fill(P);
        mom_vs_beta->Fill(P, beta->at(i));

        if (charge->at(i) > 0) {
          mom_vs_beta_pos->Fill(P, beta->at(i));
        } else if (charge->at(i) < 0) {
          mom_vs_beta_neg->Fill(P, beta->at(i));
        }

        if (pid->at(i) == 2212) {
          mom_vs_beta_proton->Fill(P, beta->at(i));
        } else if (abs(pid->at(i)) == 211) {
          mom_vs_beta_pion->Fill(P, beta->at(i));
        } else if (pid->at(i) == 11) {
          mom_vs_beta_electron->Fill(P, beta->at(i));
        }
      }

      if (pid->at(0) != 11) continue;
      // Setup scattered electron 4 vector
      TVector3 e_mu_prime_3;
      TLorentzVector e_mu_prime;
      e_mu_prime_3.SetXYZ(px->at(0), py->at(0), pz->at(0));
      e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
      double W = physics::W_calc(e_mu, e_mu_prime);
      double Q2 = physics::Q2_calc(e_mu, e_mu_prime);

      W_hist->Fill(W);
      Q2_hist->Fill(Q2);
      W_vs_q2->Fill(W, Q2);
    }

    double electron_vertex = 0.0;

    for (int j = 0; j < sc_time->size(); j++) {
      if (sc_time->size() == 0) continue;
      int index = pindex->at(j);

      if (pindex->at(j) == 0) {
        electron_vertex =
            physics::vertex_time(sc_time->at(index), sc_r->at(index), 1.0);
        continue;
      }
    }

    for (int j = 0; j < sc_time->size(); j++) {
      if (sc_time->size() == 0) continue;

      int index = pindex->at(j);

      double P_x = px->at(index) * px->at(index);
      double P_y = py->at(index) * py->at(index);
      double P_z = pz->at(index) * pz->at(index);
      P = TMath::Sqrt(P_x + P_y + P_z);

      double dt_electron = physics::deltat(electron_vertex, MASS_E, P,
                                           sc_time->at(j), sc_r->at(j));
      double dt_pion = physics::deltat(electron_vertex, MASS_PIP, P,
                                       sc_time->at(j), sc_r->at(j));
      double dt_proton = physics::deltat(electron_vertex, MASS_P, P,
                                         sc_time->at(j), sc_r->at(j));

      if (index == 0) {
        deltat_electron_0th->Fill(P, dt_electron);
      }
      if (index == 0 && pid->at(index) == 11) {
        deltat_electron_0th_ID->Fill(P, dt_electron);
      }

      if (index == 0) continue;

      deltat_pion->Fill(P, dt_pion);
      deltat_proton->Fill(P, dt_proton);
      deltat_electron->Fill(P, dt_electron);

      if (pid->at(index) != 2212) deltat_proton_antiID->Fill(P, dt_proton);

      if (pid->at(index) == 2212) {
        deltat_proton_withID->Fill(P, dt_proton);
      } else if (pid->at(index) == 211) {
        deltat_pion_withID->Fill(P, dt_pion);
      } else if (pid->at(index) == 11) {
        deltat_electron_withID->Fill(P, dt_electron);
      }
    }
  }

  out->cd();
  TDirectory *wvsq2 = out->mkdir("wvsq2");
  wvsq2->cd();
  momentum->Write();
  W_hist->Write();
  Q2_hist->Write();
  W_vs_q2->Write();

  TDirectory *mom_vs_beta = out->mkdir("mom_vs_beta");
  mom_vs_beta->cd();
  mom_vs_beta->Write();
  mom_vs_beta_pos->Write();
  mom_vs_beta_neg->Write();
  mom_vs_beta_proton->Write();
  mom_vs_beta_pion->Write();
  mom_vs_beta_electron->Write();
  mom_vs_beta_0th->Write();

  TDirectory *deltat_ftof = out->mkdir("deltat_ftof");
  deltat_ftof->cd();
  deltat_pion->Write();
  deltat_proton->Write();
  deltat_electron->Write();
  deltat_pion_withID->Write();
  deltat_proton_withID->Write();
  deltat_proton_antiID->Write();
  deltat_electron_withID->Write();
  deltat_electron_0th_ID->Write();
  deltat_electron_0th->Write();

  out->Close();
  chain->Reset();
  std::cerr << "\n" << total << "\t" << std::endl;
}
#endif
