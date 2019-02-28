/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#include "clas12branches.hpp"

Clas12Branches::Clas12Branches(TChain *tree) {
  _clas12Chain = tree;
  this->getBranches();
}

Clas12Branches::Clas12Branches(std::unique_ptr<TChain> tree) {
  _clas12Chain = tree.get();
  getBranches();
}

Clas12Branches::Clas12Branches(TChain *tree, bool MC) {
  _MC = MC;
  _clas12Chain = tree;
  getBranches();
}

Clas12Branches::~Clas12Branches() {}

void Clas12Branches::getBranches() {
  _clas12Chain->SetBranchAddress("NRUN", &_NRUN);
  _clas12Chain->SetBranchAddress("NEVENT", &_NEVENT);
  _clas12Chain->SetBranchAddress("EVNTime", &_EVNTime);
  _clas12Chain->SetBranchAddress("TYPE", &_TYPE);
  _clas12Chain->SetBranchAddress("TRG", &_TRG);
  _clas12Chain->SetBranchAddress("BCG", &_BCG);
  _clas12Chain->SetBranchAddress("STTime", &_STTime);
  _clas12Chain->SetBranchAddress("RFTime", &_RFTime);
  _clas12Chain->SetBranchAddress("Helic", &_Helic);
  _clas12Chain->SetBranchAddress("EvCAT", &_EvCAT);
  _clas12Chain->SetBranchAddress("NPGP", &_NPGP);
  _clas12Chain->SetBranchAddress("LT", &_LT);
  _clas12Chain->SetBranchAddress("PTIME", &_PTIME);

  _clas12Chain->SetBranchAddress("pid", &_pid);
  _clas12Chain->SetBranchAddress("p", &_p);
  _clas12Chain->SetBranchAddress("p2", &_p2);
  _clas12Chain->SetBranchAddress("px", &_px);
  _clas12Chain->SetBranchAddress("py", &_py);
  _clas12Chain->SetBranchAddress("pz", &_pz);
  _clas12Chain->SetBranchAddress("vx", &_vx);
  _clas12Chain->SetBranchAddress("vy", &_vy);
  _clas12Chain->SetBranchAddress("vz", &_vz);
  _clas12Chain->SetBranchAddress("charge", &_charge);
  _clas12Chain->SetBranchAddress("beta", &_beta);
  _clas12Chain->SetBranchAddress("chi2pid", &_chi2pid);
  _clas12Chain->SetBranchAddress("status", &_status);
  _clas12Chain->SetBranchAddress("ec_tot_energy", &_ec_tot_energy);
  _clas12Chain->SetBranchAddress("ec_pcal_energy", &_ec_pcal_energy);
  _clas12Chain->SetBranchAddress("ec_pcal_sec", &_ec_pcal_sec);
  _clas12Chain->SetBranchAddress("ec_pcal_time", &_ec_pcal_time);
  _clas12Chain->SetBranchAddress("ec_pcal_path", &_ec_pcal_path);
  _clas12Chain->SetBranchAddress("ec_pcal_x", &_ec_pcal_x);
  _clas12Chain->SetBranchAddress("ec_pcal_y", &_ec_pcal_y);
  _clas12Chain->SetBranchAddress("ec_pcal_z", &_ec_pcal_z);
  _clas12Chain->SetBranchAddress("ec_pcal_lu", &_ec_pcal_lu);
  _clas12Chain->SetBranchAddress("ec_pcal_lv", &_ec_pcal_lv);
  _clas12Chain->SetBranchAddress("ec_pcal_lw", &_ec_pcal_lw);
  _clas12Chain->SetBranchAddress("ec_ecin_energy", &_ec_ecin_energy);
  _clas12Chain->SetBranchAddress("ec_ecin_sec", &_ec_ecin_sec);
  _clas12Chain->SetBranchAddress("ec_ecin_time", &_ec_ecin_time);
  _clas12Chain->SetBranchAddress("ec_ecin_path", &_ec_ecin_path);
  _clas12Chain->SetBranchAddress("ec_ecin_x", &_ec_ecin_x);
  _clas12Chain->SetBranchAddress("ec_ecin_y", &_ec_ecin_y);
  _clas12Chain->SetBranchAddress("ec_ecin_z", &_ec_ecin_z);
  _clas12Chain->SetBranchAddress("ec_ecin_lu", &_ec_ecin_lu);
  _clas12Chain->SetBranchAddress("ec_ecin_lv", &_ec_ecin_lv);
  _clas12Chain->SetBranchAddress("ec_ecin_lw", &_ec_ecin_lw);
  _clas12Chain->SetBranchAddress("ec_ecout_energy", &_ec_ecout_energy);
  _clas12Chain->SetBranchAddress("ec_ecout_sec", &_ec_ecout_sec);
  _clas12Chain->SetBranchAddress("ec_ecout_time", &_ec_ecout_time);
  _clas12Chain->SetBranchAddress("ec_ecout_path", &_ec_ecout_path);
  _clas12Chain->SetBranchAddress("ec_ecout_x", &_ec_ecout_x);
  _clas12Chain->SetBranchAddress("ec_ecout_y", &_ec_ecout_y);
  _clas12Chain->SetBranchAddress("ec_ecout_z", &_ec_ecout_z);
  _clas12Chain->SetBranchAddress("ec_ecout_lu", &_ec_ecout_lu);
  _clas12Chain->SetBranchAddress("ec_ecout_lv", &_ec_ecout_lv);
  _clas12Chain->SetBranchAddress("ec_ecout_lw", &_ec_ecout_lw);
  _clas12Chain->SetBranchAddress("dc_sec", &_dc_sec);
  _clas12Chain->SetBranchAddress("dc_px", &_dc_px);
  _clas12Chain->SetBranchAddress("dc_py", &_dc_py);
  _clas12Chain->SetBranchAddress("dc_pz", &_dc_pz);
  _clas12Chain->SetBranchAddress("dc_vx", &_dc_vx);
  _clas12Chain->SetBranchAddress("dc_vy", &_dc_vy);
  _clas12Chain->SetBranchAddress("dc_vz", &_dc_vz);
  _clas12Chain->SetBranchAddress("cvt_px", &_cvt_px);
  _clas12Chain->SetBranchAddress("cvt_py", &_cvt_py);
  _clas12Chain->SetBranchAddress("cvt_pz", &_cvt_pz);
  _clas12Chain->SetBranchAddress("cvt_vx", &_cvt_vx);
  _clas12Chain->SetBranchAddress("cvt_vy", &_cvt_vy);
  _clas12Chain->SetBranchAddress("cvt_vz", &_cvt_vz);
  _clas12Chain->SetBranchAddress("cc_nphe_tot", &_cc_nphe_tot);
  _clas12Chain->SetBranchAddress("cc_ltcc_sec", &_cc_ltcc_sec);
  _clas12Chain->SetBranchAddress("cc_ltcc_nphe", &_cc_ltcc_nphe);
  _clas12Chain->SetBranchAddress("cc_ltcc_time", &_cc_ltcc_time);
  _clas12Chain->SetBranchAddress("cc_ltcc_path", &_cc_ltcc_path);
  _clas12Chain->SetBranchAddress("cc_ltcc_theta", &_cc_ltcc_theta);
  _clas12Chain->SetBranchAddress("cc_ltcc_phi", &_cc_ltcc_phi);
  _clas12Chain->SetBranchAddress("cc_htcc_sec", &_cc_htcc_sec);
  _clas12Chain->SetBranchAddress("cc_htcc_nphe", &_cc_htcc_nphe);
  _clas12Chain->SetBranchAddress("cc_htcc_time", &_cc_htcc_time);
  _clas12Chain->SetBranchAddress("cc_htcc_path", &_cc_htcc_path);
  _clas12Chain->SetBranchAddress("cc_htcc_theta", &_cc_htcc_theta);
  _clas12Chain->SetBranchAddress("cc_htcc_phi", &_cc_htcc_phi);
  _clas12Chain->SetBranchAddress("cc_rich_sec", &_cc_rich_sec);
  _clas12Chain->SetBranchAddress("cc_rich_nphe", &_cc_rich_nphe);
  _clas12Chain->SetBranchAddress("cc_rich_time", &_cc_rich_time);
  _clas12Chain->SetBranchAddress("cc_rich_path", &_cc_rich_path);
  _clas12Chain->SetBranchAddress("cc_rich_theta", &_cc_rich_theta);
  _clas12Chain->SetBranchAddress("cc_rich_phi", &_cc_rich_phi);
  _clas12Chain->SetBranchAddress("sc_ftof_1a_sec", &_sc_ftof_1a_sec);
  _clas12Chain->SetBranchAddress("sc_ftof_1a_time", &_sc_ftof_1a_time);
  _clas12Chain->SetBranchAddress("sc_ftof_1a_path", &_sc_ftof_1a_path);
  _clas12Chain->SetBranchAddress("sc_ftof_1a_energy", &_sc_ftof_1a_energy);
  _clas12Chain->SetBranchAddress("sc_ftof_1a_component", &_sc_ftof_1a_component);
  _clas12Chain->SetBranchAddress("sc_ftof_1a_x", &_sc_ftof_1a_x);
  _clas12Chain->SetBranchAddress("sc_ftof_1a_y", &_sc_ftof_1a_y);
  _clas12Chain->SetBranchAddress("sc_ftof_1a_z", &_sc_ftof_1a_z);
  _clas12Chain->SetBranchAddress("sc_ftof_1a_hx", &_sc_ftof_1a_hx);
  _clas12Chain->SetBranchAddress("sc_ftof_1a_hy", &_sc_ftof_1a_hy);
  _clas12Chain->SetBranchAddress("sc_ftof_1a_hz", &_sc_ftof_1a_hz);
  _clas12Chain->SetBranchAddress("sc_ftof_1b_sec", &_sc_ftof_1b_sec);
  _clas12Chain->SetBranchAddress("sc_ftof_1b_time", &_sc_ftof_1b_time);
  _clas12Chain->SetBranchAddress("sc_ftof_1b_path", &_sc_ftof_1b_path);
  _clas12Chain->SetBranchAddress("sc_ftof_1b_energy", &_sc_ftof_1b_energy);
  _clas12Chain->SetBranchAddress("sc_ftof_1b_component", &_sc_ftof_1b_component);
  _clas12Chain->SetBranchAddress("sc_ftof_1b_x", &_sc_ftof_1b_x);
  _clas12Chain->SetBranchAddress("sc_ftof_1b_y", &_sc_ftof_1b_y);
  _clas12Chain->SetBranchAddress("sc_ftof_1b_z", &_sc_ftof_1b_z);
  _clas12Chain->SetBranchAddress("sc_ftof_1b_hx", &_sc_ftof_1b_hx);
  _clas12Chain->SetBranchAddress("sc_ftof_1b_hy", &_sc_ftof_1b_hy);
  _clas12Chain->SetBranchAddress("sc_ftof_1b_hz", &_sc_ftof_1b_hz);
  _clas12Chain->SetBranchAddress("sc_ftof_2_sec", &_sc_ftof_2_sec);
  _clas12Chain->SetBranchAddress("sc_ftof_2_time", &_sc_ftof_2_time);
  _clas12Chain->SetBranchAddress("sc_ftof_2_path", &_sc_ftof_2_path);
  _clas12Chain->SetBranchAddress("sc_ftof_2_energy", &_sc_ftof_2_energy);
  _clas12Chain->SetBranchAddress("sc_ftof_2_component", &_sc_ftof_2_component);
  _clas12Chain->SetBranchAddress("sc_ftof_2_x", &_sc_ftof_2_x);
  _clas12Chain->SetBranchAddress("sc_ftof_2_y", &_sc_ftof_2_y);
  _clas12Chain->SetBranchAddress("sc_ftof_2_z", &_sc_ftof_2_z);
  _clas12Chain->SetBranchAddress("sc_ftof_2_hx", &_sc_ftof_2_hx);
  _clas12Chain->SetBranchAddress("sc_ftof_2_hy", &_sc_ftof_2_hy);
  _clas12Chain->SetBranchAddress("sc_ftof_2_hz", &_sc_ftof_2_hz);
  _clas12Chain->SetBranchAddress("sc_ctof_time", &_sc_ctof_time);
  _clas12Chain->SetBranchAddress("sc_ctof_path", &_sc_ctof_path);
  _clas12Chain->SetBranchAddress("sc_ctof_energy", &_sc_ctof_energy);
  _clas12Chain->SetBranchAddress("sc_ctof_component", &_sc_ctof_component);
  _clas12Chain->SetBranchAddress("sc_ctof_x", &_sc_ctof_x);
  _clas12Chain->SetBranchAddress("sc_ctof_y", &_sc_ctof_y);
  _clas12Chain->SetBranchAddress("sc_ctof_z", &_sc_ctof_z);
  _clas12Chain->SetBranchAddress("sc_ctof_hx", &_sc_ctof_hx);
  _clas12Chain->SetBranchAddress("sc_ctof_hy", &_sc_ctof_hy);
  _clas12Chain->SetBranchAddress("sc_ctof_hz", &_sc_ctof_hz);
  _clas12Chain->SetBranchAddress("sc_cnd_time", &_sc_cnd_time);
  _clas12Chain->SetBranchAddress("sc_cnd_path", &_sc_cnd_path);
  _clas12Chain->SetBranchAddress("sc_cnd_energy", &_sc_cnd_energy);
  _clas12Chain->SetBranchAddress("sc_cnd_component", &_sc_cnd_component);
  _clas12Chain->SetBranchAddress("sc_cnd_x", &_sc_cnd_x);
  _clas12Chain->SetBranchAddress("sc_cnd_y", &_sc_cnd_y);
  _clas12Chain->SetBranchAddress("sc_cnd_z", &_sc_cnd_z);
  _clas12Chain->SetBranchAddress("sc_cnd_hx", &_sc_cnd_hx);
  _clas12Chain->SetBranchAddress("sc_cnd_hy", &_sc_cnd_hy);
  _clas12Chain->SetBranchAddress("sc_cnd_hz", &_sc_cnd_hz);
  _clas12Chain->SetBranchAddress("ft_cal_energy", &_ft_cal_energy);
  _clas12Chain->SetBranchAddress("ft_cal_time", &_ft_cal_time);
  _clas12Chain->SetBranchAddress("ft_cal_path", &_ft_cal_path);
  _clas12Chain->SetBranchAddress("ft_cal_x", &_ft_cal_x);
  _clas12Chain->SetBranchAddress("ft_cal_y", &_ft_cal_y);
  _clas12Chain->SetBranchAddress("ft_cal_z", &_ft_cal_z);
  _clas12Chain->SetBranchAddress("ft_cal_dx", &_ft_cal_dx);
  _clas12Chain->SetBranchAddress("ft_cal_dy", &_ft_cal_dy);
  _clas12Chain->SetBranchAddress("ft_cal_radius", &_ft_cal_radius);
  _clas12Chain->SetBranchAddress("ft_hodo_energy", &_ft_hodo_energy);
  _clas12Chain->SetBranchAddress("ft_hodo_time", &_ft_hodo_time);
  _clas12Chain->SetBranchAddress("ft_hodo_path", &_ft_hodo_path);
  _clas12Chain->SetBranchAddress("ft_hodo_x", &_ft_hodo_x);
  _clas12Chain->SetBranchAddress("ft_hodo_y", &_ft_hodo_y);
  _clas12Chain->SetBranchAddress("ft_hodo_z", &_ft_hodo_z);
  _clas12Chain->SetBranchAddress("ft_hodo_dx", &_ft_hodo_dx);
  _clas12Chain->SetBranchAddress("ft_hodo_dy", &_ft_hodo_dy);
  _clas12Chain->SetBranchAddress("ft_hodo_radius", &_ft_hodo_radius);

  auto cachesize = 128000000U;                // 128 MBytes
  _clas12Chain->SetCacheSize(cachesize);      //<<<
  _clas12Chain->AddBranchToCache("*", true);  //<<< add all branches to the cache
}

int Clas12Branches::gpart() { return _pid->size(); }

int Clas12Branches::pid(int i) {
  if (i > _pid->size()) return -9999;
  return _pid->at(i);
}

int Clas12Branches::charge(int i) {
  if (i > _pid->size()) return -9999;
  return _charge->at(i);
}

int Clas12Branches::status(int i) {
  if (i > _pid->size()) return -9999;
  return _status->at(i);
}

int Clas12Branches::ec_pcal_sec(int i) {
  if (i > _pid->size()) return -9999;
  return _ec_pcal_sec->at(i);
}

int Clas12Branches::ec_ecin_sec(int i) {
  if (i > _pid->size()) return -9999;
  return _ec_ecin_sec->at(i);
}

int Clas12Branches::ec_ecout_sec(int i) {
  if (i > _pid->size()) return -9999;
  return _ec_ecout_sec->at(i);
}

int Clas12Branches::dc_sec(int i) {
  if (i > _pid->size()) return -9999;
  return _dc_sec->at(i);
}

int Clas12Branches::cc_ltcc_sec(int i) {
  if (i > _pid->size()) return -9999;
  return _cc_ltcc_sec->at(i);
}

int Clas12Branches::cc_htcc_sec(int i) {
  if (i > _pid->size()) return -9999;
  return _cc_htcc_sec->at(i);
}

int Clas12Branches::cc_rich_sec(int i) {
  if (i > _pid->size()) return -9999;
  return _cc_rich_sec->at(i);
}

int Clas12Branches::sc_ftof_1a_sec(int i) {
  if (i > _pid->size()) return -9999;
  return _sc_ftof_1a_sec->at(i);
}

int Clas12Branches::sc_ftof_1a_component(int i) {
  if (i > _pid->size()) return -9999;
  return _sc_ftof_1a_component->at(i);
}

int Clas12Branches::sc_ftof_1b_sec(int i) {
  if (i > _pid->size()) return -9999;
  return _sc_ftof_1b_sec->at(i);
}

int Clas12Branches::sc_ftof_1b_component(int i) {
  if (i > _pid->size()) return -9999;
  return _sc_ftof_1b_component->at(i);
}

int Clas12Branches::sc_ftof_2_sec(int i) {
  if (i > _pid->size()) return -9999;
  return _sc_ftof_2_sec->at(i);
}

int Clas12Branches::sc_ftof_2_component(int i) {
  if (i > _pid->size()) return -9999;
  return _sc_ftof_2_component->at(i);
}

int Clas12Branches::sc_ctof_component(int i) {
  if (i > _pid->size()) return -9999;
  return _sc_ctof_component->at(i);
}

int Clas12Branches::sc_cnd_component(int i) {
  if (i > _pid->size()) return -9999;
  return _sc_cnd_component->at(i);
}

float Clas12Branches::p(int i) {
  if (i > _pid->size()) return NAN;
  return _p->at(i);
}
float Clas12Branches::p2(int i) {
  if (i > _pid->size()) return NAN;
  return _p2->at(i);
}
float Clas12Branches::px(int i) {
  if (i > _pid->size()) return NAN;
  return _px->at(i);
}
float Clas12Branches::py(int i) {
  if (i > _pid->size()) return NAN;
  return _py->at(i);
}
float Clas12Branches::pz(int i) {
  if (i > _pid->size()) return NAN;
  return _pz->at(i);
}
float Clas12Branches::vx(int i) {
  if (i > _pid->size()) return NAN;
  return _vx->at(i);
}
float Clas12Branches::vy(int i) {
  if (i > _pid->size()) return NAN;
  return _vy->at(i);
}
float Clas12Branches::vz(int i) {
  if (i > _pid->size()) return NAN;
  return _vz->at(i);
}

float Clas12Branches::beta(int i) {
  if (i > _pid->size()) return NAN;
  return _beta->at(i);
}
float Clas12Branches::chi2pid(int i) {
  if (i > _pid->size()) return NAN;
  return _chi2pid->at(i);
}

float Clas12Branches::ec_tot_energy(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_tot_energy->at(i);
}
float Clas12Branches::ec_pcal_energy(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_pcal_energy->at(i);
}

float Clas12Branches::ec_pcal_time(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_pcal_time->at(i);
}
float Clas12Branches::ec_pcal_path(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_pcal_path->at(i);
}
float Clas12Branches::ec_pcal_x(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_pcal_x->at(i);
}
float Clas12Branches::ec_pcal_y(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_pcal_y->at(i);
}
float Clas12Branches::ec_pcal_z(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_pcal_z->at(i);
}
float Clas12Branches::ec_pcal_lu(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_pcal_lu->at(i);
}
float Clas12Branches::ec_pcal_lv(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_pcal_lv->at(i);
}
float Clas12Branches::ec_pcal_lw(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_pcal_lw->at(i);
}
float Clas12Branches::ec_ecin_energy(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_ecin_energy->at(i);
}

float Clas12Branches::ec_ecin_time(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_ecin_time->at(i);
}
float Clas12Branches::ec_ecin_path(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_ecin_path->at(i);
}
float Clas12Branches::ec_ecin_x(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_ecin_x->at(i);
}
float Clas12Branches::ec_ecin_y(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_ecin_y->at(i);
}
float Clas12Branches::ec_ecin_z(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_ecin_z->at(i);
}
float Clas12Branches::ec_ecin_lu(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_ecin_lu->at(i);
}
float Clas12Branches::ec_ecin_lv(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_ecin_lv->at(i);
}
float Clas12Branches::ec_ecin_lw(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_ecin_lw->at(i);
}
float Clas12Branches::ec_ecout_energy(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_ecout_energy->at(i);
}

float Clas12Branches::ec_ecout_time(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_ecout_time->at(i);
}
float Clas12Branches::ec_ecout_path(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_ecout_path->at(i);
}
float Clas12Branches::ec_ecout_x(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_ecout_x->at(i);
}
float Clas12Branches::ec_ecout_y(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_ecout_y->at(i);
}
float Clas12Branches::ec_ecout_z(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_ecout_z->at(i);
}
float Clas12Branches::ec_ecout_lu(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_ecout_lu->at(i);
}
float Clas12Branches::ec_ecout_lv(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_ecout_lv->at(i);
}
float Clas12Branches::ec_ecout_lw(int i) {
  if (i > _pid->size()) return NAN;
  return _ec_ecout_lw->at(i);
}

float Clas12Branches::dc_px(int i) {
  if (i > _pid->size()) return NAN;
  return _dc_px->at(i);
}
float Clas12Branches::dc_py(int i) {
  if (i > _pid->size()) return NAN;
  return _dc_py->at(i);
}
float Clas12Branches::dc_pz(int i) {
  if (i > _pid->size()) return NAN;
  return _dc_pz->at(i);
}
float Clas12Branches::dc_vx(int i) {
  if (i > _pid->size()) return NAN;
  return _dc_vx->at(i);
}
float Clas12Branches::dc_vy(int i) {
  if (i > _pid->size()) return NAN;
  return _dc_vy->at(i);
}
float Clas12Branches::dc_vz(int i) {
  if (i > _pid->size()) return NAN;
  return _dc_vz->at(i);
}
float Clas12Branches::cvt_px(int i) {
  if (i > _pid->size()) return NAN;
  return _cvt_px->at(i);
}
float Clas12Branches::cvt_py(int i) {
  if (i > _pid->size()) return NAN;
  return _cvt_py->at(i);
}
float Clas12Branches::cvt_pz(int i) {
  if (i > _pid->size()) return NAN;
  return _cvt_pz->at(i);
}
float Clas12Branches::cvt_vx(int i) {
  if (i > _pid->size()) return NAN;
  return _cvt_vx->at(i);
}
float Clas12Branches::cvt_vy(int i) {
  if (i > _pid->size()) return NAN;
  return _cvt_vy->at(i);
}
float Clas12Branches::cvt_vz(int i) {
  if (i > _pid->size()) return NAN;
  return _cvt_vz->at(i);
}
float Clas12Branches::cc_nphe_tot(int i) {
  if (i > _pid->size()) return NAN;
  return _cc_nphe_tot->at(i);
}

float Clas12Branches::cc_ltcc_nphe(int i) {
  if (i > _pid->size()) return NAN;
  return _cc_ltcc_nphe->at(i);
}
float Clas12Branches::cc_ltcc_time(int i) {
  if (i > _pid->size()) return NAN;
  return _cc_ltcc_time->at(i);
}
float Clas12Branches::cc_ltcc_path(int i) {
  if (i > _pid->size()) return NAN;
  return _cc_ltcc_path->at(i);
}
float Clas12Branches::cc_ltcc_theta(int i) {
  if (i > _pid->size()) return NAN;
  return _cc_ltcc_theta->at(i);
}
float Clas12Branches::cc_ltcc_phi(int i) {
  if (i > _pid->size()) return NAN;
  return _cc_ltcc_phi->at(i);
}

float Clas12Branches::cc_htcc_nphe(int i) {
  if (i > _pid->size()) return NAN;
  return _cc_htcc_nphe->at(i);
}
float Clas12Branches::cc_htcc_time(int i) {
  if (i > _pid->size()) return NAN;
  return _cc_htcc_time->at(i);
}
float Clas12Branches::cc_htcc_path(int i) {
  if (i > _pid->size()) return NAN;
  return _cc_htcc_path->at(i);
}
float Clas12Branches::cc_htcc_theta(int i) {
  if (i > _pid->size()) return NAN;
  return _cc_htcc_theta->at(i);
}
float Clas12Branches::cc_htcc_phi(int i) {
  if (i > _pid->size()) return NAN;
  return _cc_htcc_phi->at(i);
}

float Clas12Branches::cc_rich_nphe(int i) {
  if (i > _pid->size()) return NAN;
  return _cc_rich_nphe->at(i);
}
float Clas12Branches::cc_rich_time(int i) {
  if (i > _pid->size()) return NAN;
  return _cc_rich_time->at(i);
}
float Clas12Branches::cc_rich_path(int i) {
  if (i > _pid->size()) return NAN;
  return _cc_rich_path->at(i);
}
float Clas12Branches::cc_rich_theta(int i) {
  if (i > _pid->size()) return NAN;
  return _cc_rich_theta->at(i);
}
float Clas12Branches::cc_rich_phi(int i) {
  if (i > _pid->size()) return NAN;
  return _cc_rich_phi->at(i);
}

float Clas12Branches::sc_ftof_1a_time(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_1a_time->at(i);
}
float Clas12Branches::sc_ftof_1a_path(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_1a_path->at(i);
}
float Clas12Branches::sc_ftof_1a_energy(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_1a_energy->at(i);
}

float Clas12Branches::sc_ftof_1a_x(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_1a_x->at(i);
}
float Clas12Branches::sc_ftof_1a_y(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_1a_y->at(i);
}
float Clas12Branches::sc_ftof_1a_z(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_1a_z->at(i);
}
float Clas12Branches::sc_ftof_1a_hx(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_1a_hx->at(i);
}
float Clas12Branches::sc_ftof_1a_hy(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_1a_hy->at(i);
}
float Clas12Branches::sc_ftof_1a_hz(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_1a_hz->at(i);
}

float Clas12Branches::sc_ftof_1b_time(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_1b_time->at(i);
}
float Clas12Branches::sc_ftof_1b_path(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_1b_path->at(i);
}
float Clas12Branches::sc_ftof_1b_energy(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_1b_energy->at(i);
}

float Clas12Branches::sc_ftof_1b_x(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_1b_x->at(i);
}
float Clas12Branches::sc_ftof_1b_y(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_1b_y->at(i);
}
float Clas12Branches::sc_ftof_1b_z(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_1b_z->at(i);
}
float Clas12Branches::sc_ftof_1b_hx(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_1b_hx->at(i);
}
float Clas12Branches::sc_ftof_1b_hy(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_1b_hy->at(i);
}
float Clas12Branches::sc_ftof_1b_hz(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_1b_hz->at(i);
}

float Clas12Branches::sc_ftof_2_time(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_2_time->at(i);
}
float Clas12Branches::sc_ftof_2_path(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_2_path->at(i);
}
float Clas12Branches::sc_ftof_2_energy(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_2_energy->at(i);
}

float Clas12Branches::sc_ftof_2_x(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_2_x->at(i);
}
float Clas12Branches::sc_ftof_2_y(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_2_y->at(i);
}
float Clas12Branches::sc_ftof_2_z(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_2_z->at(i);
}
float Clas12Branches::sc_ftof_2_hx(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_2_hx->at(i);
}
float Clas12Branches::sc_ftof_2_hy(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_2_hy->at(i);
}
float Clas12Branches::sc_ftof_2_hz(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ftof_2_hz->at(i);
}
float Clas12Branches::sc_ctof_time(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ctof_time->at(i);
}
float Clas12Branches::sc_ctof_path(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ctof_path->at(i);
}
float Clas12Branches::sc_ctof_energy(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ctof_energy->at(i);
}

float Clas12Branches::sc_ctof_x(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ctof_x->at(i);
}
float Clas12Branches::sc_ctof_y(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ctof_y->at(i);
}
float Clas12Branches::sc_ctof_z(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ctof_z->at(i);
}
float Clas12Branches::sc_ctof_hx(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ctof_hx->at(i);
}
float Clas12Branches::sc_ctof_hy(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ctof_hy->at(i);
}
float Clas12Branches::sc_ctof_hz(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_ctof_hz->at(i);
}
float Clas12Branches::sc_cnd_time(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_cnd_time->at(i);
}
float Clas12Branches::sc_cnd_path(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_cnd_path->at(i);
}
float Clas12Branches::sc_cnd_energy(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_cnd_energy->at(i);
}

float Clas12Branches::sc_cnd_x(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_cnd_x->at(i);
}
float Clas12Branches::sc_cnd_y(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_cnd_y->at(i);
}
float Clas12Branches::sc_cnd_z(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_cnd_z->at(i);
}
float Clas12Branches::sc_cnd_hx(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_cnd_hx->at(i);
}
float Clas12Branches::sc_cnd_hy(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_cnd_hy->at(i);
}
float Clas12Branches::sc_cnd_hz(int i) {
  if (i > _pid->size()) return NAN;
  return _sc_cnd_hz->at(i);
}
float Clas12Branches::ft_cal_energy(int i) {
  if (i > _pid->size()) return NAN;
  return _ft_cal_energy->at(i);
}
float Clas12Branches::ft_cal_time(int i) {
  if (i > _pid->size()) return NAN;
  return _ft_cal_time->at(i);
}
float Clas12Branches::ft_cal_path(int i) {
  if (i > _pid->size()) return NAN;
  return _ft_cal_path->at(i);
}
float Clas12Branches::ft_cal_x(int i) {
  if (i > _pid->size()) return NAN;
  return _ft_cal_x->at(i);
}
float Clas12Branches::ft_cal_y(int i) {
  if (i > _pid->size()) return NAN;
  return _ft_cal_y->at(i);
}
float Clas12Branches::ft_cal_z(int i) {
  if (i > _pid->size()) return NAN;
  return _ft_cal_z->at(i);
}
float Clas12Branches::ft_cal_dx(int i) {
  if (i > _pid->size()) return NAN;
  return _ft_cal_dx->at(i);
}
float Clas12Branches::ft_cal_dy(int i) {
  if (i > _pid->size()) return NAN;
  return _ft_cal_dy->at(i);
}
float Clas12Branches::ft_cal_radius(int i) {
  if (i > _pid->size()) return NAN;
  return _ft_cal_radius->at(i);
}
float Clas12Branches::ft_hodo_energy(int i) {
  if (i > _pid->size()) return NAN;
  return _ft_hodo_energy->at(i);
}
float Clas12Branches::ft_hodo_time(int i) {
  if (i > _pid->size()) return NAN;
  return _ft_hodo_time->at(i);
}
float Clas12Branches::ft_hodo_path(int i) {
  if (i > _pid->size()) return NAN;
  return _ft_hodo_path->at(i);
}
float Clas12Branches::ft_hodo_x(int i) {
  if (i > _pid->size()) return NAN;
  return _ft_hodo_x->at(i);
}
float Clas12Branches::ft_hodo_y(int i) {
  if (i > _pid->size()) return NAN;
  return _ft_hodo_y->at(i);
}
float Clas12Branches::ft_hodo_z(int i) {
  if (i > _pid->size()) return NAN;
  return _ft_hodo_z->at(i);
}
float Clas12Branches::ft_hodo_dx(int i) {
  if (i > _pid->size()) return NAN;
  return _ft_hodo_dx->at(i);
}
float Clas12Branches::ft_hodo_dy(int i) {
  if (i > _pid->size()) return NAN;
  return _ft_hodo_dy->at(i);
}
float Clas12Branches::ft_hodo_radius(int i) {
  if (i > _pid->size()) return NAN;
  return _ft_hodo_radius->at(i);
}
