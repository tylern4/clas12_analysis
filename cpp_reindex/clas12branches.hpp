/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef CLAS12BRANCHES_H_GUARD
#define CLAS12BRANCHES_H_GUARD
#include <iostream>
#include <vector>
#include "TChain.h"

class Clas12Branches {
 private:
  TChain *_clas12Chain;
  bool _MC;

  size_t _entry = 0;

  int _NRUN;
  int _NEVENT;
  float _EVNTime;
  int _TYPE;
  int _TRG;
  float _BCG;
  float _STTime;
  float _RFTime;
  int _Helic;
  int _EvCAT;
  int _NPGP;
  double _LT;
  float _PTIME;

  std::vector<int> *_pid;
  std::vector<float> *_p;
  std::vector<float> *_p2;
  std::vector<float> *_px;
  std::vector<float> *_py;
  std::vector<float> *_pz;
  std::vector<float> *_vx;
  std::vector<float> *_vy;
  std::vector<float> *_vz;
  std::vector<int> *_charge;
  std::vector<float> *_beta;
  std::vector<float> *_chi2pid;
  std::vector<int> *_status;
  std::vector<float> *_ec_tot_energy;
  std::vector<float> *_ec_pcal_energy;
  std::vector<int> *_ec_pcal_sec;
  std::vector<float> *_ec_pcal_time;
  std::vector<float> *_ec_pcal_path;
  std::vector<float> *_ec_pcal_x;
  std::vector<float> *_ec_pcal_y;
  std::vector<float> *_ec_pcal_z;
  std::vector<float> *_ec_pcal_lu;
  std::vector<float> *_ec_pcal_lv;
  std::vector<float> *_ec_pcal_lw;
  std::vector<float> *_ec_ecin_energy;
  std::vector<int> *_ec_ecin_sec;
  std::vector<float> *_ec_ecin_time;
  std::vector<float> *_ec_ecin_path;
  std::vector<float> *_ec_ecin_x;
  std::vector<float> *_ec_ecin_y;
  std::vector<float> *_ec_ecin_z;
  std::vector<float> *_ec_ecin_lu;
  std::vector<float> *_ec_ecin_lv;
  std::vector<float> *_ec_ecin_lw;
  std::vector<float> *_ec_ecout_energy;
  std::vector<int> *_ec_ecout_sec;
  std::vector<float> *_ec_ecout_time;
  std::vector<float> *_ec_ecout_path;
  std::vector<float> *_ec_ecout_x;
  std::vector<float> *_ec_ecout_y;
  std::vector<float> *_ec_ecout_z;
  std::vector<float> *_ec_ecout_lu;
  std::vector<float> *_ec_ecout_lv;
  std::vector<float> *_ec_ecout_lw;
  std::vector<int> *_dc_sec;
  std::vector<float> *_dc_px;
  std::vector<float> *_dc_py;
  std::vector<float> *_dc_pz;
  std::vector<float> *_dc_vx;
  std::vector<float> *_dc_vy;
  std::vector<float> *_dc_vz;
  std::vector<float> *_cvt_px;
  std::vector<float> *_cvt_py;
  std::vector<float> *_cvt_pz;
  std::vector<float> *_cvt_vx;
  std::vector<float> *_cvt_vy;
  std::vector<float> *_cvt_vz;
  std::vector<float> *_cc_nphe_tot;
  std::vector<int> *_cc_ltcc_sec;
  std::vector<float> *_cc_ltcc_nphe;
  std::vector<float> *_cc_ltcc_time;
  std::vector<float> *_cc_ltcc_path;
  std::vector<float> *_cc_ltcc_theta;
  std::vector<float> *_cc_ltcc_phi;
  std::vector<int> *_cc_htcc_sec;
  std::vector<float> *_cc_htcc_nphe;
  std::vector<float> *_cc_htcc_time;
  std::vector<float> *_cc_htcc_path;
  std::vector<float> *_cc_htcc_theta;
  std::vector<float> *_cc_htcc_phi;
  std::vector<int> *_cc_rich_sec;
  std::vector<float> *_cc_rich_nphe;
  std::vector<float> *_cc_rich_time;
  std::vector<float> *_cc_rich_path;
  std::vector<float> *_cc_rich_theta;
  std::vector<float> *_cc_rich_phi;
  std::vector<int> *_sc_ftof_1a_sec;
  std::vector<float> *_sc_ftof_1a_time;
  std::vector<float> *_sc_ftof_1a_path;
  std::vector<float> *_sc_ftof_1a_energy;
  std::vector<int> *_sc_ftof_1a_component;
  std::vector<float> *_sc_ftof_1a_x;
  std::vector<float> *_sc_ftof_1a_y;
  std::vector<float> *_sc_ftof_1a_z;
  std::vector<float> *_sc_ftof_1a_hx;
  std::vector<float> *_sc_ftof_1a_hy;
  std::vector<float> *_sc_ftof_1a_hz;
  std::vector<int> *_sc_ftof_1b_sec;
  std::vector<float> *_sc_ftof_1b_time;
  std::vector<float> *_sc_ftof_1b_path;
  std::vector<float> *_sc_ftof_1b_energy;
  std::vector<int> *_sc_ftof_1b_component;
  std::vector<float> *_sc_ftof_1b_x;
  std::vector<float> *_sc_ftof_1b_y;
  std::vector<float> *_sc_ftof_1b_z;
  std::vector<float> *_sc_ftof_1b_hx;
  std::vector<float> *_sc_ftof_1b_hy;
  std::vector<float> *_sc_ftof_1b_hz;
  std::vector<int> *_sc_ftof_2_sec;
  std::vector<float> *_sc_ftof_2_time;
  std::vector<float> *_sc_ftof_2_path;
  std::vector<float> *_sc_ftof_2_energy;
  std::vector<int> *_sc_ftof_2_component;
  std::vector<float> *_sc_ftof_2_x;
  std::vector<float> *_sc_ftof_2_y;
  std::vector<float> *_sc_ftof_2_z;
  std::vector<float> *_sc_ftof_2_hx;
  std::vector<float> *_sc_ftof_2_hy;
  std::vector<float> *_sc_ftof_2_hz;
  std::vector<float> *_sc_ctof_time;
  std::vector<float> *_sc_ctof_path;
  std::vector<float> *_sc_ctof_energy;
  std::vector<int> *_sc_ctof_component;
  std::vector<float> *_sc_ctof_x;
  std::vector<float> *_sc_ctof_y;
  std::vector<float> *_sc_ctof_z;
  std::vector<float> *_sc_ctof_hx;
  std::vector<float> *_sc_ctof_hy;
  std::vector<float> *_sc_ctof_hz;
  std::vector<float> *_sc_cnd_time;
  std::vector<float> *_sc_cnd_path;
  std::vector<float> *_sc_cnd_energy;
  std::vector<int> *_sc_cnd_component;
  std::vector<float> *_sc_cnd_x;
  std::vector<float> *_sc_cnd_y;
  std::vector<float> *_sc_cnd_z;
  std::vector<float> *_sc_cnd_hx;
  std::vector<float> *_sc_cnd_hy;
  std::vector<float> *_sc_cnd_hz;
  std::vector<float> *_ft_cal_energy;
  std::vector<float> *_ft_cal_time;
  std::vector<float> *_ft_cal_path;
  std::vector<float> *_ft_cal_x;
  std::vector<float> *_ft_cal_y;
  std::vector<float> *_ft_cal_z;
  std::vector<float> *_ft_cal_dx;
  std::vector<float> *_ft_cal_dy;
  std::vector<float> *_ft_cal_radius;
  std::vector<float> *_ft_hodo_energy;
  std::vector<float> *_ft_hodo_time;
  std::vector<float> *_ft_hodo_path;
  std::vector<float> *_ft_hodo_x;
  std::vector<float> *_ft_hodo_y;
  std::vector<float> *_ft_hodo_z;
  std::vector<float> *_ft_hodo_dx;
  std::vector<float> *_ft_hodo_dy;
  std::vector<float> *_ft_hodo_radius;

 public:
  Clas12Branches(TChain *tree);
  Clas12Branches(std::unique_ptr<TChain> tree);
  Clas12Branches(TChain *tree, bool MC);
  ~Clas12Branches();

  void getBranches();
  void addFile(std::string fin);

  bool getNext();
  int GetEntry(long long e);
  inline size_t GetEntries() { return _clas12Chain->GetEntries(); }

  inline int NRUN() { return _NRUN; }
  inline int NEVENT() { return _NEVENT; }
  inline float EVNTime() { return _EVNTime; }
  inline int TYPE() { return _TYPE; }
  inline int TRG() { return _TRG; }
  inline float BCG() { return _BCG; }
  inline float STTime() { return _STTime; }
  inline float RFTime() { return _RFTime; }
  inline int Helic() { return _Helic; }
  inline int EvCAT() { return _EvCAT; }
  inline int NPGP() { return _NPGP; }
  inline double LT() { return _LT; }
  inline float PTIME() { return _PTIME; }

  int gpart();

  int charge(int i);
  int status(int i);
  int ec_pcal_sec(int i);
  int ec_ecin_sec(int i);
  int ec_ecout_sec(int i);
  int dc_sec(int i);
  int cc_ltcc_sec(int i);
  int cc_htcc_sec(int i);
  int cc_rich_sec(int i);
  int sc_ftof_1a_sec(int i);
  int sc_ftof_1a_component(int i);
  int sc_ftof_1b_sec(int i);
  int sc_ftof_1b_component(int i);
  int sc_ftof_2_sec(int i);
  int sc_ftof_2_component(int i);
  int sc_ctof_component(int i);
  int sc_cnd_component(int i);

  int pid(int i);
  float p(int i);
  float p2(int i);
  float px(int i);
  float py(int i);
  float pz(int i);
  float vx(int i);
  float vy(int i);
  float vz(int i);
  float beta(int i);
  float chi2pid(int i);
  float ec_tot_energy(int i);
  float ec_pcal_energy(int i);
  float ec_pcal_time(int i);
  float ec_pcal_path(int i);
  float ec_pcal_x(int i);
  float ec_pcal_y(int i);
  float ec_pcal_z(int i);
  float ec_pcal_lu(int i);
  float ec_pcal_lv(int i);
  float ec_pcal_lw(int i);
  float ec_ecin_energy(int i);
  float ec_ecin_time(int i);
  float ec_ecin_path(int i);
  float ec_ecin_x(int i);
  float ec_ecin_y(int i);
  float ec_ecin_z(int i);
  float ec_ecin_lu(int i);
  float ec_ecin_lv(int i);
  float ec_ecin_lw(int i);
  float ec_ecout_energy(int i);
  float ec_ecout_time(int i);
  float ec_ecout_path(int i);
  float ec_ecout_x(int i);
  float ec_ecout_y(int i);
  float ec_ecout_z(int i);
  float ec_ecout_lu(int i);
  float ec_ecout_lv(int i);
  float ec_ecout_lw(int i);
  float dc_px(int i);
  float dc_py(int i);
  float dc_pz(int i);
  float dc_vx(int i);
  float dc_vy(int i);
  float dc_vz(int i);
  float cvt_px(int i);
  float cvt_py(int i);
  float cvt_pz(int i);
  float cvt_vx(int i);
  float cvt_vy(int i);
  float cvt_vz(int i);
  float cc_nphe_tot(int i);
  float cc_ltcc_nphe(int i);
  float cc_ltcc_time(int i);
  float cc_ltcc_path(int i);
  float cc_ltcc_theta(int i);
  float cc_ltcc_phi(int i);
  float cc_htcc_nphe(int i);
  float cc_htcc_time(int i);
  float cc_htcc_path(int i);
  float cc_htcc_theta(int i);
  float cc_htcc_phi(int i);
  float cc_rich_nphe(int i);
  float cc_rich_time(int i);
  float cc_rich_path(int i);
  float cc_rich_theta(int i);
  float cc_rich_phi(int i);
  float sc_ftof_1a_time(int i);
  float sc_ftof_1a_path(int i);
  float sc_ftof_1a_energy(int i);
  float sc_ftof_1a_x(int i);
  float sc_ftof_1a_y(int i);
  float sc_ftof_1a_z(int i);
  float sc_ftof_1a_hx(int i);
  float sc_ftof_1a_hy(int i);
  float sc_ftof_1a_hz(int i);
  float sc_ftof_1b_time(int i);
  float sc_ftof_1b_path(int i);
  float sc_ftof_1b_energy(int i);
  float sc_ftof_1b_x(int i);
  float sc_ftof_1b_y(int i);
  float sc_ftof_1b_z(int i);
  float sc_ftof_1b_hx(int i);
  float sc_ftof_1b_hy(int i);
  float sc_ftof_1b_hz(int i);
  float sc_ftof_2_time(int i);
  float sc_ftof_2_path(int i);
  float sc_ftof_2_energy(int i);
  float sc_ftof_2_x(int i);
  float sc_ftof_2_y(int i);
  float sc_ftof_2_z(int i);
  float sc_ftof_2_hx(int i);
  float sc_ftof_2_hy(int i);
  float sc_ftof_2_hz(int i);
  float sc_ctof_time(int i);
  float sc_ctof_path(int i);
  float sc_ctof_energy(int i);
  float sc_ctof_x(int i);
  float sc_ctof_y(int i);
  float sc_ctof_z(int i);
  float sc_ctof_hx(int i);
  float sc_ctof_hy(int i);
  float sc_ctof_hz(int i);
  float sc_cnd_time(int i);
  float sc_cnd_path(int i);
  float sc_cnd_energy(int i);
  float sc_cnd_x(int i);
  float sc_cnd_y(int i);
  float sc_cnd_z(int i);
  float sc_cnd_hx(int i);
  float sc_cnd_hy(int i);
  float sc_cnd_hz(int i);
  float ft_cal_energy(int i);
  float ft_cal_time(int i);
  float ft_cal_path(int i);
  float ft_cal_x(int i);
  float ft_cal_y(int i);
  float ft_cal_z(int i);
  float ft_cal_dx(int i);
  float ft_cal_dy(int i);
  float ft_cal_radius(int i);
  float ft_hodo_energy(int i);
  float ft_hodo_time(int i);
  float ft_hodo_path(int i);
  float ft_hodo_x(int i);
  float ft_hodo_y(int i);
  float ft_hodo_z(int i);
  float ft_hodo_dx(int i);
  float ft_hodo_dy(int i);
  float ft_hodo_radius(int i);
};

#endif
