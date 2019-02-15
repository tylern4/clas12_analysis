/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef HIST_H_GUARD
#define HIST_H_GUARD
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "constants.hpp"
#include "deltat.hpp"

class Histogram {
 private:
  int bins = 500;
  double p_min = 0.0;
  double p_max = 10.0;
  double Dt_max = 10.0;
  double Dt_min = -Dt_max;
  double q2_max = 8.0;
  double w_max = 5.0;

  double zero = 0.0;
  std::string hname;
  std::string htitle;

  static const short particle_num = 4;  // 0-e 1-Pi 2-P 3-K
  std::string particle_name[particle_num] = {"e", "pi", "P", "K"};
  static const short charge_num = 2;  // 0-pos 1-neg
  std::string charge_name[charge_num] = {"positive", "negative"};
  static const short with_id_num = 3;  // 0-without 1-with 2-anti
  std::string id_name[with_id_num] = {"withoutID", "withID", "antiID"};

  // Kinematics
  TH1D *momentum;
  TH1D *W_hist;
  TH1D *Q2_hist;
  TH2D *W_vs_q2;

  static const short num_sectors = 6;
  TH2D *W_vs_q2_sec[num_sectors];
  TH1D *W_sec[num_sectors];

  TH2D *W_vs_q2_singlePi_sec[num_sectors];
  TH1D *W_singlePi_sec[num_sectors];

  TH2D *W_vs_q2_Npip_sec[num_sectors];
  TH2D *W_vs_MM_singlePi[num_sectors];
  TH1D *W_Npip_sec[num_sectors];
  TH1D *MM_Npip_sec[num_sectors];

  TH1D *MM_neutron;
  TH1D *MM_neutron_sec[num_sectors];

  TH1D *W_hist_singlePi;
  TH1D *Q2_hist_singlePi;
  TH2D *W_vs_q2_singlePi;

  // EC Sampling Fraction
  TH2D *EC_sampling_fraction;
  // EC Sampling Fraction

  // Mom vs Beta
  TH2D *momvsbeta_hist[particle_num][charge_num][with_id_num];
  TH2D *momvsbeta_vertex[with_id_num];
  // Mom vs Beta

  // Delta T
  TH2D *delta_t_hist[particle_num][charge_num][with_id_num];
  TH2D *delta_t_vertex[with_id_num];
  // Delta T

 public:
  Histogram();
  ~Histogram();

  // W and Q^2
  void makeHists_sector();
  void Fill_WvsQ2(double W, double Q2, int sec);
  void Fill_WvsQ2_singlePi(double W, double Q2, double mm, int sec);
  void Fill_WvsQ2_Npip(double W, double Q2, double mm, int sec);
  void Write_WvsQ2(TFile *out);

  // P and E
  void makeHists_MomVsBeta();
  void Fill_momentum(double P);
  void Fill_MomVsBeta_vertex(int pid, int charge, double P, double beta);
  void Fill_MomVsBeta(int pid, int charge, double P, double beta);
  void Write_MomVsBeta();

  // Delta T
  void makeHists_deltat();
  void Fill_deltat_vertex(int pid, int charge, float dt, float momentum);
  void Fill_deltat_pi(int pid, int charge, float dt, float momentum);
  void Fill_deltat_prot(int pid, int charge, float dt, float momentum);
  void Write_deltat();

  // EC Sampling Fraction
  void Fill_EC(double etot, double momentum);
  void Write_EC();
};

#endif
