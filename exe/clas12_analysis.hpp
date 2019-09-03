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
#include "histogram.hpp"
#include "reaction.hpp"

size_t run(std::shared_ptr<TChain> _chain, std::shared_ptr<Histogram> _hists, int thread_id);
size_t run_files(std::vector<std::string> inputs, std::shared_ptr<Histogram> hists, int thread_id);

size_t run_files(std::vector<std::string> inputs, std::shared_ptr<Histogram> hists, int thread_id) {
  // Called once for each thread
  // Make a new chain to process for this thread
  auto chain = std::make_shared<TChain>("clas12");
  // Add every file to the chain
  for (auto in : inputs) chain->Add(in.c_str());

  // Run the function over each thread
  return run(chain, hists, thread_id);
}

size_t run(std::shared_ptr<TChain> _chain, std::shared_ptr<Histogram> _hists, int thread_id) {
  // Get the number of events in this thread
  size_t num_of_events = (int)_chain->GetEntries();

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

    if (data->gpart() < 1) continue;
    bool elec = true;
    elec &= (data->charge(0) == NEGATIVE);
    elec &= (data->pid(0) == 11);
    if (!elec) continue;

    // If we pass electron cuts the event is processed
    total++;

    // Make a reaction class from the data given
    auto event = std::make_shared<Reaction>(data);
    auto dt = std::make_unique<Delta_T>(data);

    //// Fill_EC(data->, data->p(0));

    // For each particle in the event
    for (int part = 1; part < data->gpart(); part++) {
      dt->dt_calc(part);
      _hists->Fill_MomVsBeta(data->pid(part), data->charge(part), data->p(part), data->beta(part));
      _hists->Fill_deltat_pi(data->pid(part), data->charge(part), dt->dt_Pi(), dt->momentum(), dt->ctof());
      _hists->Fill_deltat_prot(data->pid(part), data->charge(part), dt->dt_P(), dt->momentum(), dt->ctof());

      // Check particle ID's and fill the reaction class
      if (abs(dt->dt_Pi()) < 0.5 && data->charge(part) == POSITIVE) {
        event->SetPip(part);
      } else if (abs(dt->dt_P()) < 0.5 && data->charge(part) == POSITIVE) {
        event->SetProton(part);
      } else if (abs(dt->dt_Pi()) < 0.5 && data->charge(part) == NEGATIVE) {
        event->SetPim(part);
      } else {
        event->SetOther(part);
      }
    }
    // Check the reaction class what kind of even it is and fill the appropriate histograms
    _hists->Fill_WvsQ2(event);
    if (event->SinglePip()) _hists->Fill_WvsQ2_singlePi(event);
    if (event->NeutronPip()) _hists->Fill_WvsQ2_Npip(event);
  }

  // Return the total number of events
  return total;
}
#endif
