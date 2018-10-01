/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef HIST_H_GUARD
#define HIST_H_GUARD
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
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
  static const short charge_num = 2;  // 0-un 1-pos 2-neg
  std::string charge_name[charge_num] = {"positive", "negative"};
  static const short with_id_num = 3;  // 0-without 1-with 2-anti
  std::string id_name[with_id_num] = {"withoutID", "withID", "antiID"};

  // Kinematics
  TH1D *momentum = new TH1D("mom", "mom", bins, p_min, p_max);
  TH1D *W_hist = new TH1D("W", "W", bins, zero, w_max);
  TH1D *Q2_hist = new TH1D("Q2", "Q2", bins, zero, q2_max);
  TH2D *W_vs_q2 = new TH2D("W_vs_q2", "W_vs_q2", bins, zero, w_max, bins, zero, q2_max);

  TH1D *W_hist_lower = new TH1D("W_lower", "W_lower", bins, zero, w_max);
  TH1D *Q2_hist_lower = new TH1D("Q2_lower", "Q2_lower", bins, zero, 0.4);
  TH2D *W_vs_q2_lower = new TH2D("W_vs_q2_lower", "W_vs_q2_lower", bins, zero, w_max, bins, zero, 0.4);

  TH1D *W_hist_upper = new TH1D("W_upper", "W_upper", bins, zero, w_max);
  TH1D *Q2_hist_upper = new TH1D("Q2_upper", "Q2_upper", bins, 0.4, q2_max);
  TH2D *W_vs_q2_upper = new TH2D("W_vs_q2_upper", "W_vs_q2_upper", bins, zero, w_max, bins, 0.4, q2_max);

  TH1D *W_hist_singlePi = new TH1D("W_singlePi", "W_singlePi", bins, zero, w_max);
  TH1D *Q2_hist_singlePi = new TH1D("Q2_singlePi", "Q2_singlePi", bins, zero, q2_max);
  TH2D *W_vs_q2_singlePi = new TH2D("W_vs_q2_singlePi", "W_vs_q2_singlePi", bins, zero, w_max, bins, zero, q2_max);

  TH1D *W_hist_lower_singlePi = new TH1D("W_lower_singlePi", "W_lower_singlePi", bins, zero, w_max);
  TH1D *Q2_hist_lower_singlePi = new TH1D("Q2_lower_singlePi", "Q2_lower_singlePi", bins, zero, 0.4);
  TH2D *W_vs_q2_lower_singlePi =
      new TH2D("W_vs_q2_lower_singlePi", "W_vs_q2_lower_singlePi", bins, zero, w_max, bins, zero, 0.4);

  TH1D *W_hist_upper_singlePi = new TH1D("W_upper_singlePi", "W_upper_singlePi", bins, zero, w_max);
  TH1D *Q2_hist_upper_singlePi = new TH1D("Q2_upper_singlePi", "Q2_upper_singlePi", bins, 0.4, q2_max);
  TH2D *W_vs_q2_upper_singlePi =
      new TH2D("W_vs_q2_upper_singlePi", "W_vs_q2_upper_singlePi", bins, zero, w_max, bins, 0.4, q2_max);

  // Mom vs Beta
  TH2D *momvsbeta_hist[particle_num][charge_num][with_id_num];
  TH2D *momvsbeta_vertex[with_id_num];
  // Mom vs Beta

  // Delta T
  TH2D *delta_t_hist[particle_num][charge_num][with_id_num];
  TH2D *delta_t_vertex[with_id_num];
  // Delta T

  // EC Sampling Fraction
  TH2D *EC_sampling_fraction =
      new TH2D("EC_sampling_fraction", "EC_sampling_fraction", bins, p_min, p_max, bins, zero, 1.0);
  // EC Sampling Fraction

 public:
  Histogram();
  ~Histogram();

  // W and Q^2
  void Fill_WvsQ2(double W, double Q2);
  void Fill_WvsQ2_singlePi(double W, double Q2);
  void Write_WvsQ2();

  // P and E
  void makeHists_MomVsBeta();
  void Fill_momentum(double P);
  void Fill_MomVsBeta_vertex(int pid, int charge, double P, double beta);
  void Fill_MomVsBeta(int pid, int charge, double P, double beta);
  void Write_MomVsBeta();

  // Delta T
  void makeHists_deltat();
  void Fill_deltat_vertex(int pid, int charge, double P, Delta_T *dt);
  void Fill_deltat(int pid, int charge, double P, Delta_T *dt);
  void Write_deltat();

  // EC Sampling Fraction
  void Fill_EC(double etot, double momentum);
  void Write_EC();
};

#endif
