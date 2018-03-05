/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef HIST_H_GUARD
#define HIST_H_GUARD
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"
#include "constants.hpp"
#include "deltat.hpp"

class Histogram {
 private:
  int bins = 500;
  double p_min = 0.0;
  double p_max = 10;
  double Dt_max = 10.0;
  double Dt_min = -Dt_max;

  double zero = 0.0;
  std::string hname;
  std::string htitle;

  static const short particle_num = 4;  // 0-e 1-Pi 2-P 3-K
  std::string particle_name[particle_num] = {"e", "pi", "P", "K"};
  static const short charge_num = 3;  // 0-un 1-pos 2-neg
  std::string charge_name[charge_num] = {"both", "positive", "negative"};
  static const short with_id_num = 3;  // 0-without 1-with 2-anti
  std::string id_name[with_id_num] = {"withoutID", "withID", "antiID"};

  // Kinematics
  TH1D *momentum = new TH1D("mom", "mom", bins, zero, 10);
  TH1D *W_hist = new TH1D("W", "W", bins, zero, 5);
  TH1D *Q2_hist = new TH1D("Q2", "Q2", bins, zero, 10);
  TH2D *W_vs_q2 = new TH2D("W_vs_q2", "W_vs_q2", bins, zero, 5, bins, zero, 8);

  // Mom vs Beta
  TH2D *mom_vs_beta = new TH2D("mom_vs_beta", "mom_vs_beta", bins, zero, 5, bins, zero, 1.2);
  TH2D *mom_vs_beta_0th = new TH2D("mom_vs_beta_0th", "mom_vs_beta_0th", bins, zero, 5, bins, zero, 1.2);
  TH2D *mom_vs_beta_pos = new TH2D("mom_vs_beta_pos", "mom_vs_beta_pos", bins, zero, 5, bins, zero, 1.2);
  TH2D *mom_vs_beta_neg = new TH2D("mom_vs_beta_neg", "mom_vs_beta_neg", bins, zero, 5, bins, zero, 1.2);
  TH2D *mom_vs_beta_proton =
      new TH2D("mom_vs_beta_proton", "mom_vs_beta_proton", bins, zero, 5, bins, zero, 1.2);
  TH2D *mom_vs_beta_pion = new TH2D("mom_vs_beta_pion", "mom_vs_beta_pion", bins, zero, 5, bins, zero, 1.2);
  TH2D *mom_vs_beta_electron =
      new TH2D("mom_vs_beta_electron", "mom_vs_beta_electron", bins, zero, 5, bins, zero, 1.2);
  // Mom vs Beta

  // Delta T
  TH2D *delta_t_hist[particle_num][charge_num][with_id_num];
  // Delta T
 public:
  Histogram();
  ~Histogram();

  // W and Q^2
  void Fill_WvsQ2(double W, double Q2);
  void Write_WvsQ2();

  // P and E
  void Fill_momentum(double P);
  void Fill_MomVsBeta(int pid, int charge, double P, double beta);
  void Write_MomVsBeta();

  // Delta T
  void makeHists_deltat();
  void Fill_deltat(int pid, int charge, double P, Delta_T *dt);
  void Write_deltat();
};

#endif
