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
#include "reaction.hpp"

void datahandeler(std::string fin, std::string fout) {
  auto out = std::make_unique<TFile>(fout.c_str(), "RECREATE");

  bool electron_cuts;
  // Load chain from branch h10
  auto chain = filehandeler::addFiles(fin);
  filehandeler::getBranches(chain);

  int num_of_events = (int)chain->GetEntries();
  double per = 0;

  auto hist = std::make_unique<Histogram>();

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    if (pid->size() == 0) continue;

    per = ((double)current_event / (double)num_of_events);
    if (current_event % 1000 == 0) std::cerr << "\t\t" << std::floor(100 * per) << "%\r\r" << std::flush;

    auto event = std::make_unique<Reaction>();
    event->SetElec(px->at(0), py->at(0), pz->at(0));

    if (p->at(0) != 0) hist->Fill_EC(ec_tot_energy->at(0) / p->at(0), p->at(0));

    auto dt = std::make_unique<Delta_T>(sc_ftof_1b_time->at(0), sc_ftof_1b_path->at(0), sc_ftof_1a_time->at(0),
                                        sc_ftof_1a_path->at(0), sc_ftof_2_time->at(0), sc_ftof_2_path->at(0));

    for (int part = 0; part < pid->size(); part++) {
      dt->dt_calc(p->at(part), sc_ftof_1b_time->at(part), sc_ftof_1b_path->at(part), sc_ftof_1a_time->at(part),
                  sc_ftof_1a_path->at(part), sc_ftof_2_time->at(part), sc_ftof_2_path->at(part), sc_ctof_time->at(part),
                  sc_ctof_path->at(part));

      hist->Fill_MomVsBeta(pid->at(part), charge->at(part), p->at(part), beta->at(part));
      hist->Fill_deltat_pi(pid->at(part), charge->at(part), dt->dt_Pi(), p->at(part));

      if (charge->at(part) == 1 && abs(dt->dt_Pi()) < 2.0)
        event->SetPip(px->at(part), py->at(part), pz->at(part));
      else if (charge->at(part) == 1 && abs(dt->dt_P()) < 2.0)
        event->SetProton(px->at(part), py->at(part), pz->at(part));
      else if (charge->at(part) == -1 && abs(dt->dt_Pi()) < 2.0)
        event->SetPim(px->at(part), py->at(part), pz->at(part));
      else
        event->SetOther(px->at(part), py->at(part), pz->at(part), pid->at(part));
    }

    hist->Fill_WvsQ2(event->W(), event->Q2(), ec_pcal_sec->at(0));
    if (event->SinglePip()) hist->Fill_WvsQ2_singlePi(event->W(), event->Q2(), event->MM(), ec_pcal_sec->at(0));
  }

  out->cd();
  hist->Write_EC();
  auto wvsq2 = out->mkdir("wvsq2");
  wvsq2->cd();
  hist->Write_WvsQ2();

  auto mom_vs_beta = out->mkdir("mom_vs_beta");
  mom_vs_beta->cd();
  hist->Write_MomVsBeta();

  auto deltat_ftof = out->mkdir("deltat_ftof");
  deltat_ftof->cd();
  hist->Write_deltat();

  out->Close();
  chain->Reset();
}
#endif
