/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef FILEHANDELER_H_GUARD
#define FILEHANDELER_H_GUARD
#include <vector>
#include "TChain.h"

std::vector<int> *run;
std::vector<int> *event;
std::vector<float> *torus;
std::vector<float> *solenoid;
std::vector<int> *crate;
std::vector<int> *slot;
std::vector<int> *channel;
std::vector<int> *helicity;
std::vector<int> *quartet;
std::vector<int> *value;
std::vector<int> *NRUN;
std::vector<int> *NEVENT;
std::vector<float> *EVNTime;
std::vector<int> *TYPE;
std::vector<int> *TRG;
std::vector<float> *BCG;
std::vector<float> *STTime;
std::vector<float> *RFTime;
std::vector<int> *Helic;

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

std::vector<int> *cal_pindex;
std::vector<int> *cal_detector;
std::vector<int> *cal_sector;
std::vector<int> *cal_layer;
std::vector<float> *cal_energy;
std::vector<float> *cal_time;
std::vector<float> *cal_path;
std::vector<float> *cal_x;
std::vector<float> *cal_y;
std::vector<float> *cal_z;
std::vector<float> *cal_lu;
std::vector<float> *cal_lv;
std::vector<float> *cal_lw;

std::vector<int> *chern_pindex;
std::vector<int> *chern_detector;
std::vector<int> *chern_sector;
std::vector<float> *chern_nphe;
std::vector<float> *chern_time;
std::vector<float> *chern_path;
std::vector<float> *chern_theta;
std::vector<float> *chern_phi;

std::vector<int> *fortag_pindex;
std::vector<int> *fortag_detector;
std::vector<float> *fortag_energy;
std::vector<float> *fortag_time;
std::vector<float> *fortag_path;
std::vector<float> *fortag_x;
std::vector<float> *fortag_y;
std::vector<float> *fortag_z;
std::vector<float> *fortag_dx;
std::vector<float> *fortag_dy;
std::vector<float> *fortag_radius;
std::vector<int> *fortag_size;

std::vector<int> *scint_pindex;
std::vector<int> *scint_detector;
std::vector<int> *scint_sector;
std::vector<int> *scint_layer;
std::vector<int> *scint_component;
std::vector<float> *scint_energy;
std::vector<float> *scint_time;
std::vector<float> *scint_path;

std::vector<int> *dc_sector;
std::vector<float> *dc_px;
std::vector<float> *dc_py;
std::vector<float> *dc_pz;
std::vector<float> *dc_vx;
std::vector<float> *dc_vy;
std::vector<float> *dc_vz;

std::vector<float> *cvt_px;
std::vector<float> *cvt_py;
std::vector<float> *cvt_pz;
std::vector<float> *cvt_vx;
std::vector<float> *cvt_vy;
std::vector<float> *cvt_vz;

std::vector<float> *ec_tot_energy;
std::vector<float> *ec_pcal_energy;
std::vector<int> *ec_pcal_sec;
std::vector<float> *ec_pcal_time;
std::vector<float> *ec_pcal_path;
std::vector<float> *ec_pcal_x;
std::vector<float> *ec_pcal_y;
std::vector<float> *ec_pcal_z;
std::vector<float> *ec_pcal_lu;
std::vector<float> *ec_pcal_lv;
std::vector<float> *ec_pcal_lw;

std::vector<float> *ec_ecin_energy;
std::vector<int> *ec_ecin_sec;
std::vector<float> *ec_ecin_time;
std::vector<float> *ec_ecin_path;
std::vector<float> *ec_ecin_x;
std::vector<float> *ec_ecin_y;
std::vector<float> *ec_ecin_z;
std::vector<float> *ec_ecin_lu;
std::vector<float> *ec_ecin_lv;
std::vector<float> *ec_ecin_lw;

std::vector<float> *ec_ecout_energy;
std::vector<int> *ec_ecout_sec;
std::vector<float> *ec_ecout_time;
std::vector<float> *ec_ecout_path;
std::vector<float> *ec_ecout_x;
std::vector<float> *ec_ecout_y;
std::vector<float> *ec_ecout_z;
std::vector<float> *ec_ecout_lu;
std::vector<float> *ec_ecout_lv;
std::vector<float> *ec_ecout_lw;

std::vector<float> *cc_nphe_tot;
std::vector<int> *cc_ltcc_sec;
std::vector<float> *cc_ltcc_nphe;
std::vector<float> *cc_ltcc_time;
std::vector<float> *cc_ltcc_path;
std::vector<float> *cc_ltcc_theta;
std::vector<float> *cc_ltcc_phi;
std::vector<int> *cc_htcc_sec;
std::vector<float> *cc_htcc_nphe;
std::vector<float> *cc_htcc_time;
std::vector<float> *cc_htcc_path;
std::vector<float> *cc_htcc_theta;
std::vector<float> *cc_htcc_phi;

std::vector<int> *sc_ftof_sec;
std::vector<float> *sc_ftof_time;
std::vector<float> *sc_ftof_path;
std::vector<float> *sc_ftof_layer;
std::vector<float> *sc_ftof_energy;
std::vector<float> *sc_ftof_x;
std::vector<float> *sc_ftof_y;
std::vector<float> *sc_ftof_z;
std::vector<float> *sc_ftof_hx;
std::vector<float> *sc_ftof_hy;
std::vector<float> *sc_ftof_hz;

std::vector<float> *sc_ctof_time;
std::vector<float> *sc_ctof_path;
std::vector<float> *sc_ctof_energy;
std::vector<float> *sc_ctof_x;
std::vector<float> *sc_ctof_y;
std::vector<float> *sc_ctof_z;
std::vector<float> *sc_ctof_hx;
std::vector<float> *sc_ctof_hy;
std::vector<float> *sc_ctof_hz;

std::vector<float> *ft_cal_energy;
std::vector<float> *ft_cal_time;
std::vector<float> *ft_cal_path;
std::vector<float> *ft_cal_x;
std::vector<float> *ft_cal_y;
std::vector<float> *ft_cal_z;
std::vector<float> *ft_cal_dx;
std::vector<float> *ft_cal_dy;
std::vector<float> *ft_cal_radius;

std::vector<float> *ft_hodo_energy;
std::vector<float> *ft_hodo_time;
std::vector<float> *ft_hodo_path;
std::vector<float> *ft_hodo_x;
std::vector<float> *ft_hodo_y;
std::vector<float> *ft_hodo_z;
std::vector<float> *ft_hodo_dx;
std::vector<float> *ft_hodo_dy;
std::vector<float> *ft_hodo_radius;

namespace filehandeler {
void getBranches(TTree *myTree) {
  myTree->SetBranchStatus("*", 0);

  myTree->SetBranchAddress("run", &run);
  myTree->SetBranchAddress("event", &event);
  myTree->SetBranchAddress("torus", &torus);
  myTree->SetBranchAddress("solenoid", &solenoid);
  myTree->SetBranchAddress("crate", &crate);
  myTree->SetBranchAddress("slot", &slot);
  myTree->SetBranchAddress("channel", &channel);
  myTree->SetBranchAddress("helicity", &helicity);
  myTree->SetBranchAddress("quartet", &quartet);
  myTree->SetBranchAddress("value", &value);
  myTree->SetBranchAddress("STTime", &STTime);
  myTree->SetBranchAddress("RFTime", &RFTime);
  myTree->SetBranchAddress("pid", &pid);
  myTree->SetBranchAddress("px", &px);
  myTree->SetBranchAddress("py", &py);
  myTree->SetBranchAddress("pz", &pz);
  myTree->SetBranchAddress("vx", &vx);
  myTree->SetBranchAddress("vy", &vy);
  myTree->SetBranchAddress("vz", &vz);
  myTree->SetBranchAddress("charge", &charge);
  myTree->SetBranchAddress("beta", &beta);
  myTree->SetBranchAddress("chi2pid", &chi2pid);
  myTree->SetBranchAddress("status", &status);
  myTree->SetBranchAddress("ec_tot_energy", &ec_tot_energy);
  myTree->SetBranchAddress("ec_pcal_energy", &ec_pcal_energy);
  myTree->SetBranchAddress("ec_pcal_sec", &ec_pcal_sec);
  myTree->SetBranchAddress("ec_pcal_time", &ec_pcal_time);
  myTree->SetBranchAddress("ec_pcal_path", &ec_pcal_path);
  myTree->SetBranchAddress("ec_pcal_x", &ec_pcal_x);
  myTree->SetBranchAddress("ec_pcal_y", &ec_pcal_y);
  myTree->SetBranchAddress("ec_pcal_z", &ec_pcal_z);
  myTree->SetBranchAddress("ec_pcal_lu", &ec_pcal_lu);
  myTree->SetBranchAddress("ec_pcal_lv", &ec_pcal_lv);
  myTree->SetBranchAddress("ec_pcal_lw", &ec_pcal_lw);
  myTree->SetBranchAddress("ec_ecin_energy", &ec_ecin_energy);
  myTree->SetBranchAddress("ec_ecin_sec", &ec_ecin_sec);
  myTree->SetBranchAddress("ec_ecin_time", &ec_ecin_time);
  myTree->SetBranchAddress("ec_ecin_path", &ec_ecin_path);
  myTree->SetBranchAddress("ec_ecin_x", &ec_ecin_x);
  myTree->SetBranchAddress("ec_ecin_y", &ec_ecin_y);
  myTree->SetBranchAddress("ec_ecin_z", &ec_ecin_z);
  myTree->SetBranchAddress("ec_ecin_lu", &ec_ecin_lu);
  myTree->SetBranchAddress("ec_ecin_lv", &ec_ecin_lv);
  myTree->SetBranchAddress("ec_ecin_lw", &ec_ecin_lw);
  myTree->SetBranchAddress("ec_ecout_energy", &ec_ecout_energy);
  myTree->SetBranchAddress("ec_ecout_sec", &ec_ecout_sec);
  myTree->SetBranchAddress("ec_ecout_time", &ec_ecout_time);
  myTree->SetBranchAddress("ec_ecout_path", &ec_ecout_path);
  myTree->SetBranchAddress("ec_ecout_x", &ec_ecout_x);
  myTree->SetBranchAddress("ec_ecout_y", &ec_ecout_y);
  myTree->SetBranchAddress("ec_ecout_z", &ec_ecout_z);
  myTree->SetBranchAddress("ec_ecout_lu", &ec_ecout_lu);
  myTree->SetBranchAddress("ec_ecout_lv", &ec_ecout_lv);
  myTree->SetBranchAddress("ec_ecout_lw", &ec_ecout_lw);
  myTree->SetBranchAddress("dc_sector", &dc_sector);
  myTree->SetBranchAddress("dc_px", &dc_px);
  myTree->SetBranchAddress("dc_py", &dc_py);
  myTree->SetBranchAddress("dc_pz", &dc_pz);
  myTree->SetBranchAddress("dc_vx", &dc_vx);
  myTree->SetBranchAddress("dc_vy", &dc_vy);
  myTree->SetBranchAddress("dc_vz", &dc_vz);
  myTree->SetBranchAddress("cvt_px", &cvt_px);
  myTree->SetBranchAddress("cvt_py", &cvt_py);
  myTree->SetBranchAddress("cvt_pz", &cvt_pz);
  myTree->SetBranchAddress("cvt_vx", &cvt_vx);
  myTree->SetBranchAddress("cvt_vy", &cvt_vy);
  myTree->SetBranchAddress("cvt_vz", &cvt_vz);
  myTree->SetBranchAddress("cc_nphe_tot", &cc_nphe_tot);
  myTree->SetBranchAddress("cc_ltcc_sec", &cc_ltcc_sec);
  myTree->SetBranchAddress("cc_ltcc_nphe", &cc_ltcc_nphe);
  myTree->SetBranchAddress("cc_ltcc_time", &cc_ltcc_time);
  myTree->SetBranchAddress("cc_ltcc_path", &cc_ltcc_path);
  myTree->SetBranchAddress("cc_ltcc_theta", &cc_ltcc_theta);
  myTree->SetBranchAddress("cc_ltcc_phi", &cc_ltcc_phi);
  myTree->SetBranchAddress("cc_htcc_sec", &cc_htcc_sec);
  myTree->SetBranchAddress("cc_htcc_nphe", &cc_htcc_nphe);
  myTree->SetBranchAddress("cc_htcc_time", &cc_htcc_time);
  myTree->SetBranchAddress("cc_htcc_path", &cc_htcc_path);
  myTree->SetBranchAddress("cc_htcc_theta", &cc_htcc_theta);
  myTree->SetBranchAddress("cc_htcc_phi", &cc_htcc_phi);
  myTree->SetBranchAddress("sc_ftof_sec", &sc_ftof_sec);
  myTree->SetBranchAddress("sc_ftof_time", &sc_ftof_time);
  myTree->SetBranchAddress("sc_ftof_path", &sc_ftof_path);
  myTree->SetBranchAddress("sc_ftof_layer", &sc_ftof_layer);
  myTree->SetBranchAddress("sc_ftof_energy", &sc_ftof_energy);
  myTree->SetBranchAddress("sc_ctof_time", &sc_ctof_time);
  myTree->SetBranchAddress("sc_ctof_path", &sc_ctof_path);
  myTree->SetBranchAddress("sc_ctof_energy", &sc_ctof_energy);
  myTree->SetBranchAddress("ft_cal_energy", &ft_cal_energy);
  myTree->SetBranchAddress("ft_cal_time", &ft_cal_time);
  myTree->SetBranchAddress("ft_cal_path", &ft_cal_path);
  myTree->SetBranchAddress("ft_cal_x", &ft_cal_x);
  myTree->SetBranchAddress("ft_cal_y", &ft_cal_y);
  myTree->SetBranchAddress("ft_cal_z", &ft_cal_z);
  myTree->SetBranchAddress("ft_cal_dx", &ft_cal_dx);
  myTree->SetBranchAddress("ft_cal_dy", &ft_cal_dy);
  myTree->SetBranchAddress("ft_cal_radius", &ft_cal_radius);
  myTree->SetBranchAddress("ft_hodo_energy", &ft_hodo_energy);
  myTree->SetBranchAddress("ft_hodo_time", &ft_hodo_time);
  myTree->SetBranchAddress("ft_hodo_path", &ft_hodo_path);
  myTree->SetBranchAddress("ft_hodo_x", &ft_hodo_x);
  myTree->SetBranchAddress("ft_hodo_y", &ft_hodo_y);
  myTree->SetBranchAddress("ft_hodo_z", &ft_hodo_z);
  myTree->SetBranchAddress("ft_hodo_dx", &ft_hodo_dx);
  myTree->SetBranchAddress("ft_hodo_dy", &ft_hodo_dy);
  myTree->SetBranchAddress("ft_hodo_radius", &ft_hodo_radius);
}

TChain *addFiles(char *fin) {
  TChain *clas12 = new TChain("clas12", "clas12");
  clas12->Add(fin);
  return clas12;
}
}  // namespace filehandeler

#endif
