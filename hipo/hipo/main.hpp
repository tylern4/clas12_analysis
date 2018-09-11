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
#include "../libcpp/node.h"
#include "../libcpp/reader.h"
#include "TChain.h"
#include "TH1.h"
#include "colors.hpp"
#include "constants.hpp"

void datahandeler(char *fin, char *fout) {
  auto start_full = std::chrono::high_resolution_clock::now();
  TFile *out = new TFile(fout, "RECREATE");
  TH1D *hist = new TH1D("momentum", "momentum", 500, 0, 10);

  hipo::reader reader;
  reader.open(fin);

  hipo::node<float> *px = reader.getBranch<float>("REC::Particle", "px");
  hipo::node<float> *py = reader.getBranch<float>("REC::Particle", "py");
  hipo::node<float> *pz = reader.getBranch<float>("REC::Particle", "pz");
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
