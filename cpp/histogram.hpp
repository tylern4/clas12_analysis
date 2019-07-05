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
#include "colors.hpp"
#include "constants.hpp"
#include "deltat.hpp"

using TH2D_ptr = std::shared_ptr<TH2D>;
using TH1D_ptr = std::shared_ptr<TH1D>;

class Histogram {
 protected:
  std::shared_ptr<TFile> RootOutputFile;
  TCanvas* def;

  int bins = 500;
  double p_min = 0.0;
  double p_max = 10.0;
  double Dt_max = 10.0;
  double Dt_min = -Dt_max;
  double q2_max = 8.0;
  double w_max = 5.0;

  double zero = 0.0;

  static const short particle_num = 4;  // 0-e 1-Pi 2-P 3-K
  std::string particle_name[particle_num] = {"e", "pi", "P", "K"};
  static const short charge_num = 2;  // 0-pos 1-neg
  std::string charge_name[charge_num] = {"positive", "negative"};
  static const short with_id_num = 3;  // 0-without 1-with 2-anti
  std::string id_name[with_id_num] = {"withoutID", "withID", "antiID"};

  // Kinematics
  TH1D_ptr momentum;
  TH1D_ptr W_hist;
  TH1D_ptr Q2_hist;
  TH2D_ptr W_vs_q2;

  static const short num_sectors = 6;
  TH2D_ptr W_vs_q2_sec[num_sectors];
  TH1D_ptr W_sec[num_sectors];
  TH1D_ptr W_det[3];

  TH2D_ptr W_vs_q2_singlePi_sec[num_sectors];
  TH1D_ptr W_singlePi_sec[num_sectors];

  TH2D_ptr W_vs_q2_Npip_sec[num_sectors];
  TH2D_ptr W_vs_MM_singlePi[num_sectors];
  TH1D_ptr W_Npip_sec[num_sectors];
  TH1D_ptr MM_Npip_sec[num_sectors];

  TH1D_ptr MM_neutron;
  TH1D_ptr MM_neutron_sec[num_sectors];

  TH1D_ptr W_hist_singlePi;
  TH1D_ptr Q2_hist_singlePi;
  TH2D_ptr W_vs_q2_singlePi;

  // EC Sampling Fraction
  TH2D_ptr EC_sampling_fraction;
  // EC Sampling Fraction

  // Mom vs Beta
  TH2D_ptr momvsbeta_hist[particle_num][charge_num][with_id_num];
  // Mom vs Beta

  // Delta T
  TH2D_ptr delta_t_hist[particle_num][charge_num][with_id_num][2];
  // Delta T

 public:
  Histogram(const std::string& output_file);
  ~Histogram();

  // W and Q^2
  void makeHists_sector();
  void Fill_WvsQ2(double W, double Q2, int sec);
  void Fill_WvsQ2_det(double W, double Q2, int det);
  void Fill_WvsQ2_singlePi(double W, double Q2, double mm, int sec);
  void Fill_WvsQ2_Npip(double W, double Q2, double mm, int sec);
  void Write_WvsQ2();

  // P and E
  void makeHists_MomVsBeta();
  void Fill_MomVsBeta(int pid, int charge, double P, double beta);
  void Write_MomVsBeta();

  // Delta T
  void makeHists_deltat();
  void Fill_deltat_pi(int pid, int charge, float dt, float momentum, bool fc);
  void Fill_deltat_prot(int pid, int charge, float dt, float momentum, bool fc);
  void Write_deltat();

  // EC Sampling Fraction
  void Fill_EC(double etot, double momentum);
  void Write_EC();

  //
  void Write();
};

#endif
