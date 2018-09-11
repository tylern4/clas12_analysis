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
#include <iostream>
#include <vector>
#include "TChain.h"
#include "TH1.h"
#include "colors.hpp"
#include "constants.hpp"
#include "filehandeler.hpp"

void datahandeler(char *fin, char *fout) {
  auto start_full = std::chrono::high_resolution_clock::now();
  TFile *out = new TFile(fout, "RECREATE");
  TH1D *hist = new TH1D("momentum", "momentum", 500, 0, 10);

  TChain *chain = filehandeler::addFiles(fin);
  filehandeler::getBranches(chain);

  int num_of_events = (int)chain->GetEntries();

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);

    if (px->size() == 0) continue;

    for (int i = 0; i < px->size(); i++) {
      hist->Fill(TMath::Sqrt(px->at(i) * px->at(i) + py->at(0) * py->at(0) + pz->at(0) * pz->at(0)));
    }
  }

  out->cd();
  hist->Write();
  out->Close();
  chain->Reset();
  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start_full);
  std::cout << "Elapsed time: " << elapsed_full.count() << " s" << std::endl;
  std::cout << "Sec/Event: " << elapsed_full.count() / num_of_events << std::endl;
  std::cout << "Events/Sec: " << num_of_events / elapsed_full.count() << " Hz" << std::endl;
}

#endif
