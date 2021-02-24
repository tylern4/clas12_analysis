/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD

#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "branches.hpp"
#include "colors.hpp"
#include "cuts.hpp"
#include "histogram.hpp"
#include "reaction.hpp"

template <class CutType>
size_t run(std::shared_ptr<TChain> _chain, const std::shared_ptr<Histogram>& _hists, int thread_id) {
  // Get the number of events in this thread
  size_t num_of_events = (int)_chain->GetEntries();

  float beam_energy = rga_E0;
  if (std::is_same<CutType, rga_Cuts>::value) {
    beam_energy = rga_E0;
  } else if (std::is_same<CutType, uconn_Cuts>::value) {
    beam_energy = rga_E0;
  } else if (std::is_same<CutType, rgf_Cuts>::value) {
    beam_energy = rgf_E0;
  }
  // else if (std::is_same<CutType, rgk_Cuts>::value) {
  //   beam_energy = rgk_E0;
  // }

  // Print some information for each thread
  std::cout << "=============== " << RED << "Thread " << thread_id << DEF << " =============== " << BLUE
            << num_of_events << " Events " << DEF << "===============\n";

  // Make a data object which all the branches can be accessed from
  auto data = std::make_shared<Branches12>(_chain);

  // Total number of events "Processed"
  size_t total = 0;
  // For each event
  for (size_t current_event = 0; current_event < num_of_events; current_event++) {
    // Get current event
    _chain->GetEntry(current_event);
    // If we are the 0th thread print the progress of the thread every 1000 events
    if (thread_id == 0 && current_event % 1000 == 0)
      std::cout << "\t" << (100 * current_event / num_of_events) << " %\r" << std::flush;

    auto cuts = std::make_unique<CutType>(data);
    if (!cuts->ElectronCuts()) continue;

    _hists->Fill_EC(data);
    // If we pass electron cuts the event is processed
    total++;

    // Make a reaction class from the data given
    auto event = std::make_shared<Reaction>(data, beam_energy);
    auto dt = std::make_shared<Delta_T>(data);

    // For each particle in the event
    for (int part = 1; part < data->gpart(); part++) {
      dt->dt_calc(part);
      _hists->Fill_MomVsBeta(data, part);
      _hists->Fill_deltat_pi(data, dt, part);
      _hists->Fill_deltat_prot(data, dt, part);

      // Check particle ID's and fill the reaction class
      if (cuts->IsPip(part)) {
        event->SetPip(part);
      } else if (cuts->IsProton(part)) {
        event->SetProton(part);
      } else if (cuts->IsPim(part)) {
        event->SetPim(part);
      } else {
        event->SetOther(part);
      }
    }
    // Check the reaction class what kind of even it is and fill the appropriate histograms
    _hists->Fill_WvsQ2(event);
    if (event->SinglePip()) _hists->Fill_WvsQ2_singlePip(event);
    if (event->NeutronPip()) _hists->Fill_WvsQ2_Npip(event);
  }
  std::cout << "Percent = " << 100.0 * total / num_of_events << std::endl;
  // Return the total number of events
  return num_of_events;
}
#endif
