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

    auto mc_event = std::make_shared<MCReaction>(data, beam_energy);
    csv_out += mc_event->ReacToCsv();
  }

  // Return the total number of events
  return csv_out;
}
#endif
