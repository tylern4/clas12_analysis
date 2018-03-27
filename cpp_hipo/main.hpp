/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD
#include "libcpp/reader.h"
#include "libcpp/node.h"
#include <TFile.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <fstream>
#include "TChain.h"
#include <vector>
#include "colors.hpp"
#include "constants.hpp"
#include "physics.hpp"
#include "histogram.hpp"
#include "deltat.hpp"

void datahandeler(char *fin, char *fout) {
  double energy = CLAS12_E;
  if (getenv("CLAS12_E") != NULL)
    energy = atof(getenv("CLAS12_E"));
  TLorentzVector e_mu(0.0, 0.0, energy, energy);

  TFile *out = new TFile(fout, "RECREATE");
  double P;
  bool electron_cuts;
  // Load chain from branch h10
  hipo::reader reader;
  reader.open(fin);

  hipo::node<float> *beta = reader.getBranch<float>("REC::Particle", "beta");
  hipo::node<int8_t> *charge =
      reader.getBranch<int8_t>("REC::Particle", "charge");
  hipo::node<float> *chi2pid =
      reader.getBranch<float>("REC::Particle", "chi2pid");
  hipo::node<int32_t> *pid = reader.getBranch<int32_t>("REC::Particle", "pid");
  hipo::node<float> *px = reader.getBranch<float>("REC::Particle", "px");
  hipo::node<float> *py = reader.getBranch<float>("REC::Particle", "py");
  hipo::node<float> *pz = reader.getBranch<float>("REC::Particle", "pz");
  hipo::node<int16_t> *status =
      reader.getBranch<int16_t>("REC::Particle", "status");
  hipo::node<float> *vx = reader.getBranch<float>("REC::Particle", "vx");
  hipo::node<float> *vy = reader.getBranch<float>("REC::Particle", "vy");
  hipo::node<float> *vz = reader.getBranch<float>("REC::Particle", "vz");
  hipo::node<float> *REC_Event_BCG =
      reader.getBranch<float>("REC::Event", "BCG");
  hipo::node<float> *REC_Event_EVNTime =
      reader.getBranch<float>("REC::Event", "EVNTime");
  hipo::node<int16_t> *REC_Event_EvCAT =
      reader.getBranch<int16_t>("REC::Event", "EvCAT");
  hipo::node<int8_t> *REC_Event_Helic =
      reader.getBranch<int8_t>("REC::Event", "Helic");
  hipo::node<double> *REC_Event_LT =
      reader.getBranch<double>("REC::Event", "LT");
  hipo::node<int32_t> *REC_Event_NEVENT =
      reader.getBranch<int32_t>("REC::Event", "NEVENT");
  hipo::node<int16_t> *REC_Event_NPGP =
      reader.getBranch<int16_t>("REC::Event", "NPGP");
  hipo::node<int32_t> *REC_Event_NRUN =
      reader.getBranch<int32_t>("REC::Event", "NRUN");
  hipo::node<float> *REC_Event_PTIME =
      reader.getBranch<float>("REC::Event", "PTIME");
  hipo::node<float> *REC_Event_RFTime =
      reader.getBranch<float>("REC::Event", "RFTime");
  hipo::node<float> *REC_Event_STTime =
      reader.getBranch<float>("REC::Event", "STTime");
  hipo::node<int64_t> *REC_Event_TRG =
      reader.getBranch<int64_t>("REC::Event", "TRG");
  hipo::node<int8_t> *REC_Event_TYPE =
      reader.getBranch<int8_t>("REC::Event", "TYPE");
  hipo::node<float> *REC_Calorimeter_chi2 =
      reader.getBranch<float>("REC::Calorimeter", "chi2");
  hipo::node<int8_t> *REC_Calorimeter_detector =
      reader.getBranch<int8_t>("REC::Calorimeter", "detector");
  hipo::node<float> *REC_Calorimeter_du =
      reader.getBranch<float>("REC::Calorimeter", "du");
  hipo::node<float> *REC_Calorimeter_dv =
      reader.getBranch<float>("REC::Calorimeter", "dv");
  hipo::node<float> *REC_Calorimeter_dw =
      reader.getBranch<float>("REC::Calorimeter", "dw");
  hipo::node<float> *REC_Calorimeter_energy =
      reader.getBranch<float>("REC::Calorimeter", "energy");
  hipo::node<float> *REC_Calorimeter_hx =
      reader.getBranch<float>("REC::Calorimeter", "hx");
  hipo::node<float> *REC_Calorimeter_hy =
      reader.getBranch<float>("REC::Calorimeter", "hy");
  hipo::node<float> *REC_Calorimeter_hz =
      reader.getBranch<float>("REC::Calorimeter", "hz");
  hipo::node<int16_t> *REC_Calorimeter_index =
      reader.getBranch<int16_t>("REC::Calorimeter", "index");
  hipo::node<int8_t> *REC_Calorimeter_layer =
      reader.getBranch<int8_t>("REC::Calorimeter", "layer");
  hipo::node<float> *REC_Calorimeter_lu =
      reader.getBranch<float>("REC::Calorimeter", "lu");
  hipo::node<float> *REC_Calorimeter_lv =
      reader.getBranch<float>("REC::Calorimeter", "lv");
  hipo::node<float> *REC_Calorimeter_lw =
      reader.getBranch<float>("REC::Calorimeter", "lw");
  hipo::node<float> *REC_Calorimeter_m2u =
      reader.getBranch<float>("REC::Calorimeter", "m2u");
  hipo::node<float> *REC_Calorimeter_m2v =
      reader.getBranch<float>("REC::Calorimeter", "m2v");
  hipo::node<float> *REC_Calorimeter_m2w =
      reader.getBranch<float>("REC::Calorimeter", "m2w");
  hipo::node<float> *REC_Calorimeter_path =
      reader.getBranch<float>("REC::Calorimeter", "path");
  hipo::node<int16_t> *REC_Calorimeter_pindex =
      reader.getBranch<int16_t>("REC::Calorimeter", "pindex");
  hipo::node<int8_t> *REC_Calorimeter_sector =
      reader.getBranch<int8_t>("REC::Calorimeter", "sector");
  hipo::node<int16_t> *REC_Calorimeter_status =
      reader.getBranch<int16_t>("REC::Calorimeter", "status");
  hipo::node<float> *REC_Calorimeter_time =
      reader.getBranch<float>("REC::Calorimeter", "time");
  hipo::node<float> *REC_Calorimeter_x =
      reader.getBranch<float>("REC::Calorimeter", "x");
  hipo::node<float> *REC_Calorimeter_y =
      reader.getBranch<float>("REC::Calorimeter", "y");
  hipo::node<float> *REC_Calorimeter_z =
      reader.getBranch<float>("REC::Calorimeter", "z");
  hipo::node<float> *REC_Scintillator_chi2 =
      reader.getBranch<float>("REC::Scintillator", "chi2");
  hipo::node<int16_t> *REC_Scintillator_component =
      reader.getBranch<int16_t>("REC::Scintillator", "component");
  hipo::node<int8_t> *REC_Scintillator_detector =
      reader.getBranch<int8_t>("REC::Scintillator", "detector");
  hipo::node<float> *REC_Scintillator_energy =
      reader.getBranch<float>("REC::Scintillator", "energy");
  hipo::node<float> *REC_Scintillator_hx =
      reader.getBranch<float>("REC::Scintillator", "hx");
  hipo::node<float> *REC_Scintillator_hy =
      reader.getBranch<float>("REC::Scintillator", "hy");
  hipo::node<float> *REC_Scintillator_hz =
      reader.getBranch<float>("REC::Scintillator", "hz");
  hipo::node<int16_t> *REC_Scintillator_index =
      reader.getBranch<int16_t>("REC::Scintillator", "index");
  hipo::node<int8_t> *REC_Scintillator_layer =
      reader.getBranch<int8_t>("REC::Scintillator", "layer");
  hipo::node<float> *sc_r =
      reader.getBranch<float>("REC::Scintillator", "path");
  hipo::node<int16_t> *pindex =
      reader.getBranch<int16_t>("REC::Scintillator", "pindex");
  hipo::node<int8_t> *REC_Scintillator_sector =
      reader.getBranch<int8_t>("REC::Scintillator", "sector");
  hipo::node<int16_t> *REC_Scintillator_status =
      reader.getBranch<int16_t>("REC::Scintillator", "status");
  hipo::node<float> *sc_time =
      reader.getBranch<float>("REC::Scintillator", "time");
  hipo::node<float> *REC_Scintillator_x =
      reader.getBranch<float>("REC::Scintillator", "x");
  hipo::node<float> *REC_Scintillator_y =
      reader.getBranch<float>("REC::Scintillator", "y");
  hipo::node<float> *REC_Scintillator_z =
      reader.getBranch<float>("REC::Scintillator", "z");

  int total = 0;
  Histogram *hist = new Histogram();
  int nrecords = reader.getRecordCount();

  int current_event = 0;
  while (reader.next() == true) {
    current_event++;
    if (pid->getLength() == 0)
      continue;

    if (!std::floor(current_event % 1000))
      std::cerr << "\t\t" << std::floor(current_event) << "\r\r" << std::flush;

    for (int i = 0; i < pid->getLength(); i++) {
      if (pid->getLength() == 0)
        continue;
      total++;
      double P_x = px->getValue(i) * px->getValue(i);
      double P_y = py->getValue(i) * py->getValue(i);
      double P_z = pz->getValue(i) * pz->getValue(i);

      P = TMath::Sqrt(P_x + P_y + P_z);

      if (pid->getValue(0) != 11)
        continue;
      // Setup scattered electron 4 vector
      TVector3 e_mu_prime_3;
      TLorentzVector e_mu_prime;
      e_mu_prime_3.SetXYZ(px->getValue(0), py->getValue(0), pz->getValue(0));
      e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
      double W = physics::W_calc(e_mu, e_mu_prime);
      double Q2 = physics::Q2_calc(e_mu, e_mu_prime);
      hist->Fill_WvsQ2(W, Q2);
    }

    int vertex_id;
    for (int j = 0; j < sc_time->getLength(); j++) {
      int temp = pindex->getValue(j);
      if (temp == 0) {
        vertex_id = temp;
        continue;
      }
    }

    for (int j = 0; j < sc_time->getLength(); j++) {
      if (sc_time->getLength() == 0)
        continue;
      Delta_T *dt =
          new Delta_T(sc_time->getValue(vertex_id), sc_r->getValue(vertex_id));
      int index = pindex->getValue(j);

      double P_x = px->getValue(index) * px->getValue(index);
      double P_y = py->getValue(index) * py->getValue(index);
      double P_z = pz->getValue(index) * pz->getValue(index);
      P = TMath::Sqrt(P_x + P_y + P_z);

      dt->deltat(P, sc_time->getValue(j), sc_r->getValue(j));

      if (index == 0) {
        hist->Fill_MomVsBeta_vertex(pid->getValue(index),
                                    charge->getValue(index), P,
                                    beta->getValue(index));
        hist->Fill_deltat_vertex(pid->getValue(index), charge->getValue(index),
                                 P, dt);
      } else {
        hist->Fill_MomVsBeta(pid->getValue(index), charge->getValue(index), P,
                             beta->getValue(index));
        hist->Fill_deltat(pid->getValue(index), charge->getValue(index), P, dt);
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
  std::cerr << "\n" << total << "\t" << std::endl;
}
#endif
