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
#include "TH1.h"
#include "colors.hpp"
#include "constants.hpp"
#include "node.h"
#include "reader.h"

void datahandeler(char *fin, char *fout) {
  auto start_full = std::chrono::high_resolution_clock::now();
  TFile *out = new TFile(fout, "RECREATE");
  TH1D *hist = new TH1D("momentum", "momentum", 500, 0, 10);

  hipo::reader reader;
  reader.open(fin);

  hipo::node<float> *beta = reader.getBranch<float>("REC::Particle", "beta");
  hipo::node<int8_t> *charge = reader.getBranch<int8_t>("REC::Particle", "charge");
  hipo::node<float> *chi2pid = reader.getBranch<float>("REC::Particle", "chi2pid");
  hipo::node<int32_t> *pid = reader.getBranch<int32_t>("REC::Particle", "pid");
  hipo::node<float> *px = reader.getBranch<float>("REC::Particle", "px");
  hipo::node<float> *py = reader.getBranch<float>("REC::Particle", "py");
  hipo::node<float> *pz = reader.getBranch<float>("REC::Particle", "pz");
  hipo::node<int16_t> *status = reader.getBranch<int16_t>("REC::Particle", "status");
  hipo::node<float> *vx = reader.getBranch<float>("REC::Particle", "vx");
  hipo::node<float> *vy = reader.getBranch<float>("REC::Particle", "vy");
  hipo::node<float> *vz = reader.getBranch<float>("REC::Particle", "vz");
  hipo::node<float> *evnt_BCG = reader.getBranch<float>("REC::Event", "BCG");
  hipo::node<float> *evnt_EVNTime = reader.getBranch<float>("REC::Event", "EVNTime");
  hipo::node<int16_t> *evnt_EvCAT = reader.getBranch<int16_t>("REC::Event", "EvCAT");
  hipo::node<int8_t> *evnt_Helic = reader.getBranch<int8_t>("REC::Event", "Helic");
  hipo::node<double> *evnt_LT = reader.getBranch<double>("REC::Event", "LT");
  hipo::node<int32_t> *evnt_NEVENT = reader.getBranch<int32_t>("REC::Event", "NEVENT");
  hipo::node<int16_t> *evnt_NPGP = reader.getBranch<int16_t>("REC::Event", "NPGP");
  hipo::node<int32_t> *evnt_NRUN = reader.getBranch<int32_t>("REC::Event", "NRUN");
  hipo::node<float> *evnt_PTIME = reader.getBranch<float>("REC::Event", "PTIME");
  hipo::node<float> *evnt_RFTime = reader.getBranch<float>("REC::Event", "RFTime");
  hipo::node<float> *evnt_STTime = reader.getBranch<float>("REC::Event", "STTime");
  hipo::node<int64_t> *evnt_TRG = reader.getBranch<int64_t>("REC::Event", "TRG");
  hipo::node<int8_t> *evnt_TYPE = reader.getBranch<int8_t>("REC::Event", "TYPE");
  hipo::node<float> *ec_chi2 = reader.getBranch<float>("REC::Calorimeter", "chi2");
  hipo::node<int8_t> *ec_detector = reader.getBranch<int8_t>("REC::Calorimeter", "detector");
  hipo::node<float> *ec_du = reader.getBranch<float>("REC::Calorimeter", "du");
  hipo::node<float> *ec_dv = reader.getBranch<float>("REC::Calorimeter", "dv");
  hipo::node<float> *ec_dw = reader.getBranch<float>("REC::Calorimeter", "dw");
  hipo::node<float> *ec_energy = reader.getBranch<float>("REC::Calorimeter", "energy");
  hipo::node<float> *ec_hx = reader.getBranch<float>("REC::Calorimeter", "hx");
  hipo::node<float> *ec_hy = reader.getBranch<float>("REC::Calorimeter", "hy");
  hipo::node<float> *ec_hz = reader.getBranch<float>("REC::Calorimeter", "hz");
  hipo::node<int16_t> *ec_index = reader.getBranch<int16_t>("REC::Calorimeter", "index");
  hipo::node<int8_t> *ec_layer = reader.getBranch<int8_t>("REC::Calorimeter", "layer");
  hipo::node<float> *ec_lu = reader.getBranch<float>("REC::Calorimeter", "lu");
  hipo::node<float> *ec_lv = reader.getBranch<float>("REC::Calorimeter", "lv");
  hipo::node<float> *ec_lw = reader.getBranch<float>("REC::Calorimeter", "lw");
  hipo::node<float> *ec_m2u = reader.getBranch<float>("REC::Calorimeter", "m2u");
  hipo::node<float> *ec_m2v = reader.getBranch<float>("REC::Calorimeter", "m2v");
  hipo::node<float> *ec_m2w = reader.getBranch<float>("REC::Calorimeter", "m2w");
  hipo::node<float> *ec_path = reader.getBranch<float>("REC::Calorimeter", "path");
  hipo::node<int16_t> *ec_pindex = reader.getBranch<int16_t>("REC::Calorimeter", "pindex");
  hipo::node<int8_t> *ec_sector = reader.getBranch<int8_t>("REC::Calorimeter", "sector");
  hipo::node<int16_t> *ec_status = reader.getBranch<int16_t>("REC::Calorimeter", "status");
  hipo::node<float> *ec_time = reader.getBranch<float>("REC::Calorimeter", "time");
  hipo::node<float> *ec_x = reader.getBranch<float>("REC::Calorimeter", "x");
  hipo::node<float> *ec_y = reader.getBranch<float>("REC::Calorimeter", "y");
  hipo::node<float> *ec_z = reader.getBranch<float>("REC::Calorimeter", "z");
  hipo::node<float> *sc_chi2 = reader.getBranch<float>("REC::Scintillator", "chi2");
  hipo::node<int16_t> *sc_component = reader.getBranch<int16_t>("REC::Scintillator", "component");
  hipo::node<int8_t> *sc_detector = reader.getBranch<int8_t>("REC::Scintillator", "detector");
  hipo::node<float> *sc_energy = reader.getBranch<float>("REC::Scintillator", "energy");
  hipo::node<float> *sc_hx = reader.getBranch<float>("REC::Scintillator", "hx");
  hipo::node<float> *sc_hy = reader.getBranch<float>("REC::Scintillator", "hy");
  hipo::node<float> *sc_hz = reader.getBranch<float>("REC::Scintillator", "hz");
  hipo::node<int16_t> *sc_index = reader.getBranch<int16_t>("REC::Scintillator", "index");
  hipo::node<int8_t> *sc_layer = reader.getBranch<int8_t>("REC::Scintillator", "layer");
  hipo::node<float> *sc_path = reader.getBranch<float>("REC::Scintillator", "path");
  hipo::node<int16_t> *sc_pindex = reader.getBranch<int16_t>("REC::Scintillator", "pindex");
  hipo::node<int8_t> *sc_sector = reader.getBranch<int8_t>("REC::Scintillator", "sector");
  hipo::node<int16_t> *sc_status = reader.getBranch<int16_t>("REC::Scintillator", "status");
  hipo::node<float> *sc_time = reader.getBranch<float>("REC::Scintillator", "time");
  hipo::node<float> *sc_x = reader.getBranch<float>("REC::Scintillator", "x");
  hipo::node<float> *sc_y = reader.getBranch<float>("REC::Scintillator", "y");
  hipo::node<float> *sc_z = reader.getBranch<float>("REC::Scintillator", "z");
  hipo::node<float> *cc_chi2 = reader.getBranch<float>("REC::Cherenkov", "chi2");
  hipo::node<int8_t> *cc_detector = reader.getBranch<int8_t>("REC::Cherenkov", "detector");
  hipo::node<float> *cc_dphi = reader.getBranch<float>("REC::Cherenkov", "dphi");
  hipo::node<float> *cc_dtheta = reader.getBranch<float>("REC::Cherenkov", "dtheta");
  hipo::node<int16_t> *cc_index = reader.getBranch<int16_t>("REC::Cherenkov", "index");
  hipo::node<float> *cc_nphe = reader.getBranch<float>("REC::Cherenkov", "nphe");
  hipo::node<float> *cc_path = reader.getBranch<float>("REC::Cherenkov", "path");
  hipo::node<float> *cc_phi = reader.getBranch<float>("REC::Cherenkov", "phi");
  hipo::node<int16_t> *cc_pindex = reader.getBranch<int16_t>("REC::Cherenkov", "pindex");
  hipo::node<int8_t> *cc_sector = reader.getBranch<int8_t>("REC::Cherenkov", "sector");
  hipo::node<int16_t> *cc_status = reader.getBranch<int16_t>("REC::Cherenkov", "status");
  hipo::node<float> *cc_theta = reader.getBranch<float>("REC::Cherenkov", "theta");
  hipo::node<float> *cc_time = reader.getBranch<float>("REC::Cherenkov", "time");
  hipo::node<float> *cc_x = reader.getBranch<float>("REC::Cherenkov", "x");
  hipo::node<float> *cc_y = reader.getBranch<float>("REC::Cherenkov", "y");
  hipo::node<float> *cc_z = reader.getBranch<float>("REC::Cherenkov", "z");

  int num_of_events = 0;
  while (reader.next() == true) {
    num_of_events++;
    if (px->getLength() == 0) continue;

    for (int i = 0; i < px->getLength(); i++) {
      hist->Fill(TMath::Sqrt(px->getValue(i) * px->getValue(i) + py->getValue(0) * py->getValue(0) +
                             pz->getValue(0) * pz->getValue(0)));
    }
  }

  out->cd();
  hist->Write();
  out->Close();

  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start_full);
  std::cout << "Elapsed time: " << elapsed_full.count() << " s" << std::endl;
  std::cout << "Sec/Event: " << elapsed_full.count() / num_of_events << std::endl;
  std::cout << "Events/Sec: " << num_of_events / elapsed_full.count() << " Hz" << std::endl;
}

#endif
