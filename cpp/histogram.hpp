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
  double zero = 0.0;

  TH1D *momentum = new TH1D("mom", "mom", bins, zero, 10);
  TH1D *W_hist = new TH1D("W", "W", bins, zero, 5);
  TH1D *Q2_hist = new TH1D("Q2", "Q2", bins, zero, 10);
  TH2D *W_vs_q2 = new TH2D("W_vs_q2", "W_vs_q2", bins, zero, 5, bins, zero, 8);
  TH2D *mom_vs_beta =
      new TH2D("mom_vs_beta", "mom_vs_beta", bins, zero, 5, bins, zero, 1.2);
  TH2D *mom_vs_beta_0th = new TH2D("mom_vs_beta_0th", "mom_vs_beta_0th", bins,
                                   zero, 5, bins, zero, 1.2);
  TH2D *mom_vs_beta_pos = new TH2D("mom_vs_beta_pos", "mom_vs_beta_pos", bins,
                                   zero, 5, bins, zero, 1.2);
  TH2D *mom_vs_beta_neg = new TH2D("mom_vs_beta_neg", "mom_vs_beta_neg", bins,
                                   zero, 5, bins, zero, 1.2);
  TH2D *mom_vs_beta_proton =
      new TH2D("mom_vs_beta_proton", "mom_vs_beta_proton", bins, zero, 5, bins,
               zero, 1.2);
  TH2D *mom_vs_beta_pion = new TH2D("mom_vs_beta_pion", "mom_vs_beta_pion",
                                    bins, zero, 5, bins, zero, 1.2);
  TH2D *mom_vs_beta_electron =
      new TH2D("mom_vs_beta_electron", "mom_vs_beta_electron", bins, zero, 5,
               bins, zero, 1.2);

  TH2D *deltat_proton =
      new TH2D("deltat_proton", "#Deltat assuming mass of proton", bins, p_min,
               p_max, bins, -10.0, 10.0);
  TH2D *deltat_pion = new TH2D("deltat_pion", "#Deltat assuming mass of pion",
                               bins, p_min, p_max, bins, -10.0, 10.0);
  TH2D *deltat_electron =
      new TH2D("deltat_electron", "#Deltat assuming mass of electron", bins,
               p_min, p_max, bins, -10.0, 10.0);
  TH2D *deltat_proton_withID = new TH2D(
      "deltat_proton_withID", "#Deltat assuming mass of proton with ID", bins,
      p_min, p_max, bins, -10.0, 10.0);

  TH2D *deltat_proton_antiID = new TH2D(
      "deltat_proton_antiID", "#Deltat assuming mass of proton anti ID", bins,
      p_min, p_max, bins, -10.0, 10.0);

  TH2D *deltat_pion_withID =
      new TH2D("deltat_pion_withID", "#Deltat assuming mass of pion with ID",
               bins, p_min, p_max, bins, -10.0, 10.0);
  TH2D *deltat_electron_withID = new TH2D(
      "deltat_electron_withID", "#Deltat assuming mass of electron with ID",
      bins, p_min, p_max, bins, -10.0, 10.0);
  TH2D *deltat_electron_0th = new TH2D(
      "deltat_electron_0th", "#Deltat assuming mass of electron at 0th", bins,
      p_min, p_max, bins, -10.0, 10.0);
  TH2D *deltat_electron_0th_ID =
      new TH2D("deltat_electron_0th_ID",
               "#Deltat assuming mass of electron with ID at 0th", bins, p_min,
               p_max, bins, -10.0, 10.0);

 public:
  Histogram();
  ~Histogram();

  // W and Q^2
  void Fill_WvsQ2(double W, double Q2);
  void Write_WvsQ2();

  // P and E
  void Fill_momentum(double P);
  void Fill_MomVsBeta(int pid, int charge, double P, double beta);
  void Fill_mom_vs_beta_0th(double P, double beta);
  void Write_MomVsBeta();

  // Delta T
  void Fill_deltat(int pid, double P, double dt_proton, double dt_pion,
                   double dt_electron);
  void Fill_deltat(int pid, double P, Delta_T *dt);
  void Fill_mom_vs_beta_0th(int pid, double P, double dt);
  void Fill_mom_vs_beta_0th(int pid, double P, Delta_T *dt);
  void Write_deltat();
};

#endif
