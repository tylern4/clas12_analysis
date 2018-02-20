/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef HIST_H_GUARD
#define HIST_H_GUARD
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"

TH1D *momentum = new TH1D("mom", "mom", 500, 0, 10);
TH1D *W_hist = new TH1D("W", "W", 500, 0, 5);
TH1D *Q2_hist = new TH1D("Q2", "Q2", 500, 0, 10);
TH2D *W_vs_q2 = new TH2D("W_vs_q2", "W_vs_q2", 500, 0, 5, 500, 0, 8);
TH2D *mom_vs_beta =
    new TH2D("mom_vs_beta", "mom_vs_beta", 500, 0, 5, 500, 0.0, 1.2);
TH2D *mom_vs_beta_0th =
    new TH2D("mom_vs_beta_0th", "mom_vs_beta_0th", 500, 0, 5, 500, 0.0, 1.2);
TH2D *mom_vs_beta_pos =
    new TH2D("mom_vs_beta_pos", "mom_vs_beta_pos", 500, 0, 5, 500, 0.0, 1.2);
TH2D *mom_vs_beta_neg =
    new TH2D("mom_vs_beta_neg", "mom_vs_beta_neg", 500, 0, 5, 500, 0.0, 1.2);
TH2D *mom_vs_beta_proton = new TH2D("mom_vs_beta_proton", "mom_vs_beta_proton",
                                    500, 0, 5, 500, 0.0, 1.2);
TH2D *mom_vs_beta_pion =
    new TH2D("mom_vs_beta_pion", "mom_vs_beta_pion", 500, 0, 5, 500, 0.0, 1.2);
TH2D *mom_vs_beta_electron = new TH2D(
    "mom_vs_beta_electron", "mom_vs_beta_electron", 500, 0, 5, 500, 0.0, 1.2);

TH2D *deltat_proton =
    new TH2D("deltat_proton", "#Deltat assuming mass of proton", 500, 0.0, 10.0,
             500, -10.0, 10.0);
TH2D *deltat_pion = new TH2D("deltat_pion", "#Deltat assuming mass of pion",
                             500, 0.0, 10.0, 500, -10.0, 10.0);
TH2D *deltat_electron =
    new TH2D("deltat_electron", "#Deltat assuming mass of electron", 500, 0.0,
             10.0, 500, -10.0, 10.0);
TH2D *deltat_proton_withID =
    new TH2D("deltat_proton_withID", "#Deltat assuming mass of proton with ID",
             500, 0.0, 10.0, 500, -10.0, 10.0);

TH2D *deltat_proton_antiID =
    new TH2D("deltat_proton_antiID", "#Deltat assuming mass of proton anti ID",
             500, 0.0, 10.0, 500, -10.0, 10.0);

TH2D *deltat_pion_withID =
    new TH2D("deltat_pion_withID", "#Deltat assuming mass of pion with ID", 500,
             0.0, 10.0, 500, -10.0, 10.0);
TH2D *deltat_electron_withID = new TH2D(
    "deltat_electron_withID", "#Deltat assuming mass of electron with ID", 500,
    0.0, 10.0, 500, -10.0, 10.0);
TH2D *deltat_electron_0th =
    new TH2D("deltat_electron_0th", "#Deltat assuming mass of electron at 0th",
             500, 0.0, 10.0, 500, -10.0, 10.0);
TH2D *deltat_electron_0th_ID =
    new TH2D("deltat_electron_0th_ID",
             "#Deltat assuming mass of electron with ID at 0th", 500, 0.0, 10.0,
             500, -10.0, 10.0);

#endif
