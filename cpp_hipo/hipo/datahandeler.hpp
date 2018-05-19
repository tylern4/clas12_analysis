/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef DATAHANDELER_H_GUARD
#define DATAHANDELER_H_GUARD
#include "colors.hpp"
#include "constants.hpp"
#include "deltat.hpp"
#include "histogram.hpp"
#include "main.hpp"
#include "physics.hpp"

class DataHandeler {
 private:
  Histogram *hist;
  TFile *out;
  TLorentzVector *e_mu;
  double energy = CLAS12_E;
  int total = 0;
  std::vector<std::string> input_files;

 public:
  DataHandeler(std::vector<std::string> fin, std::string fout);
  ~DataHandeler();
  void file_handeler(std::string fin);
  void run();
};

#endif
