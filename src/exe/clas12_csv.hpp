/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef CSV_H_GUARD
#define CSV_H_GUARD

#include <iostream>
#include "TFile.h"
#include "branches.hpp"
#include "colors.hpp"
#include "cuts.hpp"
#include "deltat.hpp"
#include "reaction.hpp"

std::string run(std::shared_ptr<TChain> _chain, int thread_id);
std::string run_files(std::vector<std::string> inputs, int thread_id);

std::string run_files(std::vector<std::string> inputs, int thread_id) {
  // Called once for each thread
  // Make a new chain to process for this thread
  auto chain = std::make_shared<TChain>("clas12");
  // Add every file to the chain
  for (auto in : inputs) chain->Add(in.c_str());

  // Run the function over each thread
  return run(chain, thread_id);
}

std::string run(std::shared_ptr<TChain> _chain, int thread_id) {
  std::string csv_out;
  // Get the number of events in this thread
  size_t num_of_events = (int)_chain->GetEntries();
  float beam_energy = NAN;
  if (getenv("BEAM_E") != NULL) beam_energy = atof(getenv("BEAM_E"));

  // Print some information for each thread
  std::cout << "=============== " << RED << "Thread " << thread_id << DEF << " =============== " << BLUE
            << num_of_events << " Events " << DEF << "===============\n";

  // Make a data object which all the branches can be accessed from
  auto data = std::make_shared<Branches12>(_chain, true);

  // Total number of events "Processed"
  size_t total = 0;
  // For each event
  for (size_t current_event = 0; current_event < num_of_events; current_event++) {
    // Get current event
    data->GetEntry(current_event);

    // If we are the 0th thread print the progress of the thread every 1000 events
    if (thread_id == 0 && current_event % 1000 == 0)
      std::cout << BLUE << "\t" << (100 * current_event / num_of_events) << " %\r" << DEF << std::flush;

    total++;
    if (data->mc_npart() > 1) {
      // Make a reaction class from the data given
      auto mc_event = std::make_shared<MCReaction>(data, beam_energy);
      for (int part = 1; part < data->mc_npart(); part++) {
        // Check particle ID's and fill the reaction class
        if (data->mc_pid(part) == PIP) {
          mc_event->SetPip(part);
        } else if (data->mc_pid(part) == PROTON) {
          mc_event->SetProton(part);
        } else if (data->mc_pid(part) == PIM) {
          mc_event->SetPim(part);
        } else {
          mc_event->SetOther(part);
        }
      }
    }

    auto cuts = std::make_shared<Cuts>(data);
    if (!cuts->ElectronCuts()) continue;

    auto event = std::make_shared<Reaction>(data, beam_energy);
    auto dt = std::make_shared<Delta_T>(data);
    // For each particle in the event
    // For each particle in the event
    for (int part = 1; part < data->gpart(); part++) {
      dt->dt_calc(part);
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
    // if (event->NeutronPip())
  }

  // Return the total number of events
  return csv_out;
}
#endif
