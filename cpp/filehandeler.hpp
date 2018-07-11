/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef FILEHANDELER_H_GUARD
#define FILEHANDELER_H_GUARD
#include <vector>
#include "TChain.h"

std::vector<int> *REC_Event_NRUN;
std::vector<int> *REC_Event_NEVENT;
std::vector<float> *REC_Event_EVNTime;
std::vector<int> *REC_Event_TYPE;
std::vector<int> *REC_Event_TRG;
std::vector<float> *REC_Event_BCG;
std::vector<float> *REC_Event_STTime;
std::vector<float> *REC_Event_RFTime;
std::vector<int> *REC_Event_Helic;
std::vector<int> *pid;
std::vector<float> *px;
std::vector<float> *py;
std::vector<float> *pz;
std::vector<float> *vx;
std::vector<float> *vy;
std::vector<float> *vz;
std::vector<int> *charge;
std::vector<float> *beta;
std::vector<int> *status;
std::vector<int> *ec_index;
std::vector<int> *ec_pindex;
std::vector<int> *ec_detector;
std::vector<int> *ec_sector;
std::vector<int> *ec_layer;
std::vector<float> *ec_energy;
std::vector<float> *ec_time;
std::vector<float> *ec_path;
std::vector<float> *ec_chi2;
std::vector<float> *ec_x;
std::vector<float> *ec_y;
std::vector<float> *ec_z;
std::vector<float> *ec_lu;
std::vector<float> *ec_lv;
std::vector<float> *ec_lw;
std::vector<int> *ec_status;
std::vector<int> *cc_pindex;
std::vector<int> *cc_detector;
std::vector<int> *cc_sector;
std::vector<float> *cc_nphe;
std::vector<float> *cc_time;
std::vector<float> *cc_path;
std::vector<float> *cc_x;
std::vector<float> *cc_y;
std::vector<float> *cc_z;
std::vector<float> *cc_theta;
std::vector<float> *cc_phi;
std::vector<int> *ft_index;
std::vector<int> *ft_pindex;
std::vector<int> *ft_detector;
std::vector<float> *ft_energy;
std::vector<float> *ft_time;
std::vector<float> *ft_path;
std::vector<float> *ft_chi2;
std::vector<float> *ft_x;
std::vector<float> *ft_y;
std::vector<float> *ft_z;
std::vector<float> *ft_dx;
std::vector<float> *ft_dy;
std::vector<float> *ft_radius;
std::vector<int> *ft_size;
std::vector<int> *ft_status;
std::vector<int> *sc_pindex;
std::vector<int> *sc_detector;
std::vector<int> *sc_sector;
std::vector<int> *sc_layer;
std::vector<int> *sc_component;
std::vector<float> *sc_energy;
std::vector<float> *sc_time;
std::vector<float> *sc_path;

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
  myTree->SetBranchAddress("REC_Particle_status", &status);

  myTree->SetBranchAddress("REC_Event_NRUN", &REC_Event_NRUN);
  myTree->SetBranchAddress("REC_Event_NEVENT", &REC_Event_NEVENT);
  myTree->SetBranchAddress("REC_Event_EVNTime", &REC_Event_EVNTime);
  myTree->SetBranchAddress("REC_Event_TYPE", &REC_Event_TYPE);
  myTree->SetBranchAddress("REC_Event_TRG", &REC_Event_TRG);
  myTree->SetBranchAddress("REC_Event_BCG", &REC_Event_BCG);
  myTree->SetBranchAddress("REC_Event_STTime", &REC_Event_STTime);
  myTree->SetBranchAddress("REC_Event_RFTime", &REC_Event_RFTime);
  myTree->SetBranchAddress("REC_Event_Helic", &REC_Event_Helic);
  myTree->SetBranchAddress("REC_Calorimeter_pindex", &ec_pindex);
  myTree->SetBranchAddress("REC_Calorimeter_detector", &ec_detector);
  myTree->SetBranchAddress("REC_Calorimeter_sector", &ec_sector);
  myTree->SetBranchAddress("REC_Calorimeter_layer", &ec_layer);
  myTree->SetBranchAddress("REC_Calorimeter_energy", &ec_energy);
  myTree->SetBranchAddress("REC_Calorimeter_time", &ec_time);
  myTree->SetBranchAddress("REC_Calorimeter_path", &ec_path);
  myTree->SetBranchAddress("REC_Calorimeter_x", &ec_x);
  myTree->SetBranchAddress("REC_Calorimeter_y", &ec_y);
  myTree->SetBranchAddress("REC_Calorimeter_z", &ec_z);
  myTree->SetBranchAddress("REC_Calorimeter_lu", &ec_lu);
  myTree->SetBranchAddress("REC_Calorimeter_lv", &ec_lv);
  myTree->SetBranchAddress("REC_Calorimeter_lw", &ec_lw);
  myTree->SetBranchAddress("REC_Cherenkov_pindex", &cc_pindex);
  myTree->SetBranchAddress("REC_Cherenkov_detector", &cc_detector);
  myTree->SetBranchAddress("REC_Cherenkov_sector", &cc_sector);
  myTree->SetBranchAddress("REC_Cherenkov_nphe", &cc_nphe);
  myTree->SetBranchAddress("REC_Cherenkov_time", &cc_time);
  myTree->SetBranchAddress("REC_Cherenkov_path", &cc_path);
  myTree->SetBranchAddress("REC_Cherenkov_x", &cc_x);
  myTree->SetBranchAddress("REC_Cherenkov_y", &cc_y);
  myTree->SetBranchAddress("REC_Cherenkov_z", &cc_z);
  myTree->SetBranchAddress("REC_Cherenkov_theta", &cc_theta);
  myTree->SetBranchAddress("REC_Cherenkov_phi", &cc_phi);
  myTree->SetBranchAddress("REC_ForwardTagger_pindex", &ft_pindex);
  myTree->SetBranchAddress("REC_ForwardTagger_detector", &ft_detector);
  myTree->SetBranchAddress("REC_ForwardTagger_energy", &ft_energy);
  myTree->SetBranchAddress("REC_ForwardTagger_time", &ft_time);
  myTree->SetBranchAddress("REC_ForwardTagger_path", &ft_path);
  myTree->SetBranchAddress("REC_ForwardTagger_chi2", &ft_chi2);
  myTree->SetBranchAddress("REC_ForwardTagger_x", &ft_x);
  myTree->SetBranchAddress("REC_ForwardTagger_y", &ft_y);
  myTree->SetBranchAddress("REC_ForwardTagger_z", &ft_z);
  myTree->SetBranchAddress("REC_ForwardTagger_dx", &ft_dx);
  myTree->SetBranchAddress("REC_ForwardTagger_dy", &ft_dy);
  myTree->SetBranchAddress("REC_ForwardTagger_radius", &ft_radius);
  myTree->SetBranchAddress("REC_ForwardTagger_size", &ft_size);
  myTree->SetBranchAddress("REC_ForwardTagger_status", &ft_status);
  myTree->SetBranchAddress("REC_Scintillator_pindex", &sc_pindex);
  myTree->SetBranchAddress("REC_Scintillator_detector", &sc_detector);
  myTree->SetBranchAddress("REC_Scintillator_sector", &sc_sector);
  myTree->SetBranchAddress("REC_Scintillator_layer", &sc_layer);
  myTree->SetBranchAddress("REC_Scintillator_component", &sc_component);
  myTree->SetBranchAddress("REC_Scintillator_energy", &sc_energy);
  myTree->SetBranchAddress("REC_Scintillator_time", &sc_time);
  myTree->SetBranchAddress("REC_Scintillator_path", &sc_path);

  myTree->SetBranchStatus("*", 1);
}

TChain *addFiles(char *fin) {
  TChain *clas12 = new TChain("clas12", "clas12");
  clas12->Add(fin);
  return clas12;
}
}  // namespace filehandeler

#endif
