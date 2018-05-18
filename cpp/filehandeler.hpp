/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef FILEHANDELER_H_GUARD
#define FILEHANDELER_H_GUARD
#include <vector>
#include "TChain.h"

std::vector<int> *pid;
std::vector<float> *px;
std::vector<float> *py;
std::vector<float> *pz;
std::vector<float> *vx;
std::vector<float> *vy;
std::vector<float> *vz;
std::vector<int> *charge;
std::vector<float> *beta;
std::vector<float> *chi2pid;
std::vector<int> *status;

std::vector<int> *sc_pindex;
std::vector<int> *sc_detector;
std::vector<float> *sc_time;
std::vector<float> *sc_r;

std::vector<int> *ec_pindex;
std::vector<float> *etot;

namespace filehandeler {
void getBranches(TTree *myTree) {
  myTree->SetBranchAddress("REC_Particle_pid", &pid);
  myTree->SetBranchAddress("REC_Particle_px", &px);
  myTree->SetBranchAddress("REC_Particle_py", &py);
  myTree->SetBranchAddress("REC_Particle_pz", &pz);
  myTree->SetBranchAddress("REC_Particle_vx", &vx);
  myTree->SetBranchAddress("REC_Particle_vy", &vy);
  myTree->SetBranchAddress("REC_Particle_vz", &vz);
  myTree->SetBranchAddress("REC_Particle_charge", &charge);
  myTree->SetBranchAddress("REC_Particle_beta", &beta);
  myTree->SetBranchAddress("REC_Particle_chi2pid", &chi2pid);
  myTree->SetBranchAddress("REC_Particle_status", &status);

  myTree->SetBranchAddress("REC_Scintillator_pindex", &sc_pindex);
  myTree->SetBranchAddress("REC_Scintillator_detector", &sc_detector);
  myTree->SetBranchAddress("REC_Scintillator_time", &sc_time);
  myTree->SetBranchAddress("REC_Scintillator_path", &sc_r);

  myTree->SetBranchAddress("REC_Calorimeter_pindex", &ec_pindex);
  myTree->SetBranchAddress("REC_Calorimeter_energy", &etot);

  myTree->SetBranchStatus("*", 1);
}

TChain *addFiles(char *fin) {
  TChain *clas12 = new TChain("clas12", "clas12");
  clas12->Add(fin);
  return clas12;
}
}  // namespace filehandeler

#endif
