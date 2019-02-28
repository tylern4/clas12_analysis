/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD
#include <TFile.h>
#include <TLorentzVector.h>
#include <chrono>
#include <fstream>
#include <vector>
#include "TChain.h"
#include "clas12branches.hpp"
#include "colors.hpp"
#include "constants.hpp"
#include "deltat.hpp"
#include "histogram.hpp"
#include "physics.hpp"
#include "reaction.hpp"

void datahandeler(std::string fin, std::string fout) {
  auto out = std::make_unique<TFile>(fout.c_str(), "RECREATE");

  bool electron_cuts;
  // Load chain from branch h10
  TChain *clas12 = new TChain("clas12", "clas12");
  clas12->Add(fin.c_str());

  auto data = std::make_shared<Clas12Branches>(clas12);
  auto hist = std::make_shared<Histogram>();

  size_t num_of_events = data->GetEntries();
  auto start_full = std::chrono::high_resolution_clock::now();
  for (size_t current_event = 0; current_event < num_of_events; current_event++) {
    clas12->GetEntry(current_event);
    if (data->gpart() == 0) continue;
    if (data->charge(0) != -1) continue;

    if (current_event % 1000 == 0)
      std::cerr << "\t\t" << std::floor(100 * ((double)current_event / (double)num_of_events)) << "%\r\r" << std::flush;

    auto event = std::make_unique<Reaction>();
    event->SetElec(data->px(0), data->py(0), data->pz(0));

    if (data->p(0) != 0) hist->Fill_EC(data->ec_tot_energy(0) / data->p(0), data->p(0));

    auto dt = std::make_unique<Delta_T>(data);

    for (int part = 1; part < data->gpart(); part++) {
      hist->Fill_MomVsBeta(data->pid(part), data->charge(part), data->p(part), data->beta(part));
      hist->Fill_deltat_pi(data->pid(part), data->charge(part), dt->dt_Pi(part), data->p(part));
      hist->Fill_deltat_prot(data->pid(part), data->charge(part), dt->dt_P(part), data->p(part));

      if (data->charge(part) == 1 && abs(dt->dt_Pi(part)) < 0.5)
        event->SetPip(data->px(part), data->py(part), data->pz(part));
      else if (data->charge(part) == 1 && abs(dt->dt_P(part)) < 0.5)
        event->SetProton(data->px(part), data->py(part), data->pz(part));
      else if (data->charge(part) == -1 && abs(dt->dt_Pi(part)) < 0.5)
        event->SetPim(data->px(part), data->py(part), data->pz(part));
      else
        event->SetOther(data->px(part), data->py(part), data->pz(part), data->pid(part));
    }

    hist->Fill_WvsQ2(event->W(), event->Q2(), data->ec_pcal_sec(0));
    if (event->SinglePip()) hist->Fill_WvsQ2_singlePi(event->W(), event->Q2(), event->MM(), data->ec_pcal_sec(0));
    if (event->SinglePip() && event->MM() > 0.85 && event->MM() < 1.1)
      hist->Fill_WvsQ2_Npip(event->W(), event->Q2(), event->MM(), data->ec_pcal_sec(0));
  }

  out->cd();
  hist->Write_EC();
  auto wvsq2 = out->mkdir("wvsq2");
  wvsq2->cd();
  hist->Write_WvsQ2(out.get());

  auto mom_vs_beta = out->mkdir("mom_vs_beta");
  mom_vs_beta->cd();
  hist->Write_MomVsBeta();

  auto deltat_ftof = out->mkdir("deltat_ftof");
  deltat_ftof->cd();
  hist->Write_deltat();

  out->Close();
  clas12->Reset();

  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start_full);
  std::cout << "Elapsed time for " << num_of_events << " events: " << elapsed_full.count() << " s" << std::endl;
  std::cout << "Events/Sec: " << (num_of_events / elapsed_full.count()) << " Hz" << std::endl;
}
#endif
