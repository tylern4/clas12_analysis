/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef FILEHANDELER_H_GUARD
#define FILEHANDELER_H_GUARD
#include <vector>
#include "TChain.h"
using v_int = std::vector<int> *;
using v_float = std::vector<float> *;

int NRUN;
int NEVENT;
float EVNTime;
int TYPE;
int TRG;
float BCG;
float STTime;
float RFTime;
int Helic;
int EvCAT;
int NPGP;
double LT;
float PTIME;

v_int pid;
v_float p;
v_float p2;
v_float px;
v_float py;
v_float pz;
v_float vx;
v_float vy;
v_float vz;
v_int charge;
v_float beta;
v_float chi2pid;
v_int status;

v_int dc_sec;
v_float dc_px;
v_float dc_py;
v_float dc_pz;
v_float dc_vx;
v_float dc_vy;
v_float dc_vz;

v_float cvt_px;
v_float cvt_py;
v_float cvt_pz;
v_float cvt_vx;
v_float cvt_vy;
v_float cvt_vz;

v_float ec_tot_energy;
v_float ec_pcal_energy;
v_int ec_pcal_sec;
v_float ec_pcal_time;
v_float ec_pcal_path;
v_float ec_pcal_x;
v_float ec_pcal_y;
v_float ec_pcal_z;
v_float ec_pcal_hx;
v_float ec_pcal_hy;
v_float ec_pcal_hz;
v_float ec_pcal_lu;
v_float ec_pcal_lv;
v_float ec_pcal_lw;
v_float ec_pcal_du;
v_float ec_pcal_dv;
v_float ec_pcal_dw;
v_float ec_pcal_m2u;
v_float ec_pcal_m2v;
v_float ec_pcal_m2w;
v_float ec_pcal_m3u;
v_float ec_pcal_m3v;
v_float ec_pcal_m3w;

v_float ec_ecin_energy;
v_int ec_ecin_sec;
v_float ec_ecin_time;
v_float ec_ecin_path;
v_float ec_ecin_x;
v_float ec_ecin_y;
v_float ec_ecin_z;
v_float ec_ecin_hx;
v_float ec_ecin_hy;
v_float ec_ecin_hz;
v_float ec_ecin_lu;
v_float ec_ecin_lv;
v_float ec_ecin_lw;
v_float ec_ecin_du;
v_float ec_ecin_dv;
v_float ec_ecin_dw;
v_float ec_ecin_m2u;
v_float ec_ecin_m2v;
v_float ec_ecin_m2w;
v_float ec_ecin_m3u;
v_float ec_ecin_m3v;
v_float ec_ecin_m3w;

v_float ec_ecout_energy;
v_int ec_ecout_sec;
v_float ec_ecout_time;
v_float ec_ecout_path;
v_float ec_ecout_x;
v_float ec_ecout_y;
v_float ec_ecout_z;
v_float ec_ecout_hx;
v_float ec_ecout_hy;
v_float ec_ecout_hz;
v_float ec_ecout_lu;
v_float ec_ecout_lv;
v_float ec_ecout_lw;
v_float ec_ecout_du;
v_float ec_ecout_dv;
v_float ec_ecout_dw;
v_float ec_ecout_m2u;
v_float ec_ecout_m2v;
v_float ec_ecout_m2w;
v_float ec_ecout_m3u;
v_float ec_ecout_m3v;
v_float ec_ecout_m3w;

v_float cc_nphe_tot;
v_int cc_ltcc_sec;
v_float cc_ltcc_nphe;
v_float cc_ltcc_time;
v_float cc_ltcc_path;
v_float cc_ltcc_theta;
v_float cc_ltcc_phi;
v_float cc_ltcc_x;
v_float cc_ltcc_y;
v_float cc_ltcc_z;
v_int cc_htcc_sec;
v_float cc_htcc_nphe;
v_float cc_htcc_time;
v_float cc_htcc_path;
v_float cc_htcc_theta;
v_float cc_htcc_phi;
v_float cc_htcc_x;
v_float cc_htcc_y;
v_float cc_htcc_z;
v_int cc_rich_sec;
v_float cc_rich_nphe;
v_float cc_rich_time;
v_float cc_rich_path;
v_float cc_rich_theta;
v_float cc_rich_phi;
v_float cc_rich_x;
v_float cc_rich_y;
v_float cc_rich_z;

v_int sc_ftof_1a_sec;
v_float sc_ftof_1a_time;
v_float sc_ftof_1a_path;
v_float sc_ftof_1a_energy;
v_int sc_ftof_1a_component;
v_float sc_ftof_1a_x;
v_float sc_ftof_1a_y;
v_float sc_ftof_1a_z;
v_float sc_ftof_1a_hx;
v_float sc_ftof_1a_hy;
v_float sc_ftof_1a_hz;

v_int sc_ftof_1b_sec;
v_float sc_ftof_1b_time;
v_float sc_ftof_1b_path;
v_float sc_ftof_1b_energy;
v_int sc_ftof_1b_component;
v_float sc_ftof_1b_x;
v_float sc_ftof_1b_y;
v_float sc_ftof_1b_z;
v_float sc_ftof_1b_hx;
v_float sc_ftof_1b_hy;
v_float sc_ftof_1b_hz;

v_int sc_ftof_2_sec;
v_float sc_ftof_2_time;
v_float sc_ftof_2_path;
v_float sc_ftof_2_energy;
v_int sc_ftof_2_component;
v_float sc_ftof_2_x;
v_float sc_ftof_2_y;
v_float sc_ftof_2_z;
v_float sc_ftof_2_hx;
v_float sc_ftof_2_hy;
v_float sc_ftof_2_hz;

v_float sc_ctof_time;
v_float sc_ctof_path;
v_float sc_ctof_energy;
v_int sc_ctof_component;
v_float sc_ctof_x;
v_float sc_ctof_y;
v_float sc_ctof_z;
v_float sc_ctof_hx;
v_float sc_ctof_hy;
v_float sc_ctof_hz;

v_float sc_cnd_time;
v_float sc_cnd_path;
v_float sc_cnd_energy;
v_int sc_cnd_component;
v_float sc_cnd_x;
v_float sc_cnd_y;
v_float sc_cnd_z;
v_float sc_cnd_hx;
v_float sc_cnd_hy;
v_float sc_cnd_hz;

v_float ft_cal_energy;
v_float ft_cal_time;
v_float ft_cal_path;
v_float ft_cal_x;
v_float ft_cal_y;
v_float ft_cal_z;
v_float ft_cal_dx;
v_float ft_cal_dy;
v_float ft_cal_radius;

v_float ft_hodo_energy;
v_float ft_hodo_time;
v_float ft_hodo_path;
v_float ft_hodo_x;
v_float ft_hodo_y;
v_float ft_hodo_z;
v_float ft_hodo_dx;
v_float ft_hodo_dy;
v_float ft_hodo_radius;

v_int MC_pid;
v_float MC_helicity;
v_float MC_px;
v_float MC_py;
v_float MC_pz;
v_float MC_vx;
v_float MC_vy;
v_float MC_vz;
v_float MC_vt;

v_int Lund_pid;
v_float Lund_px;
v_float Lund_py;
v_float Lund_pz;
v_float Lund_E;
v_float Lund_vx;
v_float Lund_vy;
v_float Lund_vz;
v_float Lund_ltime;

v_float CovMat_11;
v_float CovMat_12;
v_float CovMat_13;
v_float CovMat_14;
v_float CovMat_15;
v_float CovMat_22;
v_float CovMat_23;
v_float CovMat_24;
v_float CovMat_25;
v_float CovMat_33;
v_float CovMat_34;
v_float CovMat_35;
v_float CovMat_44;
v_float CovMat_45;
v_float CovMat_55;

v_int VertDoca_index1_vec;
v_int VertDoca_index2_vec;
v_float VertDoca_x_vec;
v_float VertDoca_y_vec;
v_float VertDoca_z_vec;
v_float VertDoca_x1_vec;
v_float VertDoca_y1_vec;
v_float VertDoca_z1_vec;
v_float VertDoca_cx1_vec;
v_float VertDoca_cy1_vec;
v_float VertDoca_cz1_vec;
v_float VertDoca_x2_vec;
v_float VertDoca_y2_vec;
v_float VertDoca_z2_vec;
v_float VertDoca_cx2_vec;
v_float VertDoca_cy2_vec;
v_float VertDoca_cz2_vec;
v_float VertDoca_r_vec;

v_int traj_pindex_vec;
v_int traj_index_vec;
v_int traj_detId_vec;
v_int traj_q_vec;
v_float traj_x_vec;
v_float traj_y_vec;
v_float traj_z_vec;
v_float traj_cx_vec;
v_float traj_cy_vec;
v_float traj_cz_vec;
v_float traj_pathlength_vec;

namespace filehandeler {
void getBranches(TTree *myTree) {
  myTree->SetBranchStatus("*", 0);

  myTree->SetBranchAddress("NRUN", &NRUN);
  myTree->SetBranchAddress("NEVENT", &NEVENT);
  myTree->SetBranchAddress("EVNTime", &EVNTime);
  myTree->SetBranchAddress("TYPE", &TYPE);
  myTree->SetBranchAddress("TRG", &TRG);
  myTree->SetBranchAddress("BCG", &BCG);
  myTree->SetBranchAddress("STTime", &STTime);
  myTree->SetBranchAddress("RFTime", &RFTime);
  myTree->SetBranchAddress("Helic", &Helic);
  myTree->SetBranchAddress("EvCAT", &EvCAT);
  myTree->SetBranchAddress("NPGP", &NPGP);
  myTree->SetBranchAddress("LT", &LT);
  myTree->SetBranchAddress("PTIME", &PTIME);
  myTree->SetBranchAddress("pid", &pid);
  myTree->SetBranchAddress("p", &p);
  myTree->SetBranchAddress("p2", &p2);
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
  myTree->SetBranchAddress("dc_sec", &dc_sec);
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
  myTree->SetBranchAddress("cc_rich_sec", &cc_rich_sec);
  myTree->SetBranchAddress("cc_rich_nphe", &cc_rich_nphe);
  myTree->SetBranchAddress("cc_rich_time", &cc_rich_time);
  myTree->SetBranchAddress("cc_rich_path", &cc_rich_path);
  myTree->SetBranchAddress("cc_rich_theta", &cc_rich_theta);
  myTree->SetBranchAddress("cc_rich_phi", &cc_rich_phi);
  myTree->SetBranchAddress("sc_ftof_1a_sec", &sc_ftof_1a_sec);
  myTree->SetBranchAddress("sc_ftof_1a_time", &sc_ftof_1a_time);
  myTree->SetBranchAddress("sc_ftof_1a_path", &sc_ftof_1a_path);
  myTree->SetBranchAddress("sc_ftof_1a_energy", &sc_ftof_1a_energy);
  myTree->SetBranchAddress("sc_ftof_1a_component", &sc_ftof_1a_component);
  myTree->SetBranchAddress("sc_ftof_1a_x", &sc_ftof_1a_x);
  myTree->SetBranchAddress("sc_ftof_1a_y", &sc_ftof_1a_y);
  myTree->SetBranchAddress("sc_ftof_1a_z", &sc_ftof_1a_z);
  myTree->SetBranchAddress("sc_ftof_1a_hx", &sc_ftof_1a_hx);
  myTree->SetBranchAddress("sc_ftof_1a_hy", &sc_ftof_1a_hy);
  myTree->SetBranchAddress("sc_ftof_1a_hz", &sc_ftof_1a_hz);
  myTree->SetBranchAddress("sc_ftof_1b_sec", &sc_ftof_1b_sec);
  myTree->SetBranchAddress("sc_ftof_1b_time", &sc_ftof_1b_time);
  myTree->SetBranchAddress("sc_ftof_1b_path", &sc_ftof_1b_path);
  myTree->SetBranchAddress("sc_ftof_1b_energy", &sc_ftof_1b_energy);
  myTree->SetBranchAddress("sc_ftof_1b_component", &sc_ftof_1b_component);
  myTree->SetBranchAddress("sc_ftof_1b_x", &sc_ftof_1b_x);
  myTree->SetBranchAddress("sc_ftof_1b_y", &sc_ftof_1b_y);
  myTree->SetBranchAddress("sc_ftof_1b_z", &sc_ftof_1b_z);
  myTree->SetBranchAddress("sc_ftof_1b_hx", &sc_ftof_1b_hx);
  myTree->SetBranchAddress("sc_ftof_1b_hy", &sc_ftof_1b_hy);
  myTree->SetBranchAddress("sc_ftof_1b_hz", &sc_ftof_1b_hz);
  myTree->SetBranchAddress("sc_ftof_2_sec", &sc_ftof_2_sec);
  myTree->SetBranchAddress("sc_ftof_2_time", &sc_ftof_2_time);
  myTree->SetBranchAddress("sc_ftof_2_path", &sc_ftof_2_path);
  myTree->SetBranchAddress("sc_ftof_2_energy", &sc_ftof_2_energy);
  myTree->SetBranchAddress("sc_ftof_2_component", &sc_ftof_2_component);
  myTree->SetBranchAddress("sc_ftof_2_x", &sc_ftof_2_x);
  myTree->SetBranchAddress("sc_ftof_2_y", &sc_ftof_2_y);
  myTree->SetBranchAddress("sc_ftof_2_z", &sc_ftof_2_z);
  myTree->SetBranchAddress("sc_ftof_2_hx", &sc_ftof_2_hx);
  myTree->SetBranchAddress("sc_ftof_2_hy", &sc_ftof_2_hy);
  myTree->SetBranchAddress("sc_ftof_2_hz", &sc_ftof_2_hz);
  myTree->SetBranchAddress("sc_ctof_time", &sc_ctof_time);
  myTree->SetBranchAddress("sc_ctof_path", &sc_ctof_path);
  myTree->SetBranchAddress("sc_ctof_energy", &sc_ctof_energy);
  myTree->SetBranchAddress("sc_ctof_component", &sc_ctof_component);
  myTree->SetBranchAddress("sc_ctof_x", &sc_ctof_x);
  myTree->SetBranchAddress("sc_ctof_y", &sc_ctof_y);
  myTree->SetBranchAddress("sc_ctof_z", &sc_ctof_z);
  myTree->SetBranchAddress("sc_ctof_hx", &sc_ctof_hx);
  myTree->SetBranchAddress("sc_ctof_hy", &sc_ctof_hy);
  myTree->SetBranchAddress("sc_ctof_hz", &sc_ctof_hz);
  myTree->SetBranchAddress("sc_cnd_time", &sc_cnd_time);
  myTree->SetBranchAddress("sc_cnd_path", &sc_cnd_path);
  myTree->SetBranchAddress("sc_cnd_energy", &sc_cnd_energy);
  myTree->SetBranchAddress("sc_cnd_component", &sc_cnd_component);
  myTree->SetBranchAddress("sc_cnd_x", &sc_cnd_x);
  myTree->SetBranchAddress("sc_cnd_y", &sc_cnd_y);
  myTree->SetBranchAddress("sc_cnd_z", &sc_cnd_z);
  myTree->SetBranchAddress("sc_cnd_hx", &sc_cnd_hx);
  myTree->SetBranchAddress("sc_cnd_hy", &sc_cnd_hy);
  myTree->SetBranchAddress("sc_cnd_hz", &sc_cnd_hz);
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

TChain *addFiles(std::string fin) {
  TChain *clas12 = new TChain("clas12", "clas12");
  clas12->Add(fin.c_str());
  return clas12;
}
}  // namespace filehandeler

#endif
