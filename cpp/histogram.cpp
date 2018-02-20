/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "histogram.hpp"
#include "constants.hpp"

Histogram::Histogram() {}

Histogram::~Histogram() {
  delete momentum;
  delete W_hist;
  delete Q2_hist;
  delete W_vs_q2;
  delete mom_vs_beta;
  delete mom_vs_beta_0th;
  delete mom_vs_beta_pos;
  delete mom_vs_beta_neg;
  delete mom_vs_beta_proton;
  delete mom_vs_beta_pion;
  delete mom_vs_beta_electron;
  delete deltat_proton;
  delete deltat_pion;
  delete deltat_electron;
  delete deltat_proton_withID;
  delete deltat_proton_antiID;
  delete deltat_pion_withID;
  delete deltat_electron_withID;
  delete deltat_electron_0th;
  delete deltat_electron_0th_ID;
}

// W and Q^2
void Histogram::Fill_WvsQ2(double W, double Q2) {
  W_vs_q2->Fill(W, Q2);
  W_hist->Fill(W);
  Q2_hist->Fill(Q2);
}

void Histogram::Write_WvsQ2() {
  W_vs_q2->SetXTitle("W (GeV)");
  W_vs_q2->SetYTitle("Q^{2} (GeV^{2})");
  W_vs_q2->SetOption("COLZ");
  W_vs_q2->Write();

  W_hist->SetXTitle("W (GeV)");
  W_hist->Write();

  Q2_hist->SetXTitle("Q^{2} (GeV^{2})");
  Q2_hist->Write();
}

void Histogram::Fill_MomVsBeta(int pid, int charge, double P, double beta) {
  if (beta != 0) {
    momentum->Fill(P);
    mom_vs_beta->Fill(P, beta);

    if (charge > 0) {
      mom_vs_beta_pos->Fill(P, beta);
    } else if (charge < 0) {
      mom_vs_beta_neg->Fill(P, beta);
    }
    if (pid == PROTON) {
      mom_vs_beta_proton->Fill(P, beta);
    } else if (pid == PIP) {
      mom_vs_beta_pion->Fill(P, beta);
    } else if (pid == ELECTRON) {
      mom_vs_beta_electron->Fill(P, beta);
    }
  }
}
void Histogram::Fill_mom_vs_beta_0th(double P, double beta) {
  mom_vs_beta_0th->Fill(P, beta);
}

void Histogram::Fill_deltat(int pid, double P, double dt_proton, double dt_pion,
                            double dt_electron) {
  deltat_pion->Fill(P, dt_pion);
  deltat_proton->Fill(P, dt_proton);
  deltat_electron->Fill(P, dt_electron);

  if (pid != PROTON) deltat_proton_antiID->Fill(P, dt_proton);

  if (pid == PROTON) {
    deltat_proton_withID->Fill(P, dt_proton);
  } else if (pid == PIP) {
    deltat_pion_withID->Fill(P, dt_pion);
  } else if (pid == ELECTRON) {
    deltat_electron_withID->Fill(P, dt_electron);
  }
}
void Histogram::Fill_mom_vs_beta_0th(int pid, double P, double dt) {
  deltat_electron_0th->Fill(P, dt);
  if (pid == ELECTRON) {
    deltat_electron_0th_ID->Fill(P, dt);
  }
}

void Histogram::Write_MomVsBeta() {
  momentum->SetXTitle("Momentum (GeV)");
  momentum->Write();
  mom_vs_beta->SetXTitle("Momentum (GeV)");
  mom_vs_beta->SetYTitle("#beta");
  mom_vs_beta->SetOption("COLZ");
  mom_vs_beta->Write();
  mom_vs_beta_pos->SetXTitle("Momentum (GeV)");
  mom_vs_beta_pos->SetYTitle("#beta");
  mom_vs_beta_pos->SetOption("COLZ");
  mom_vs_beta_pos->Write();
  mom_vs_beta_neg->SetXTitle("Momentum (GeV)");
  mom_vs_beta_neg->SetYTitle("#beta");
  mom_vs_beta_neg->SetOption("COLZ");
  mom_vs_beta_neg->Write();

  mom_vs_beta_proton->SetXTitle("Momentum (GeV)");
  mom_vs_beta_proton->SetYTitle("#beta");
  mom_vs_beta_proton->SetOption("COLZ");
  mom_vs_beta_proton->Write();

  mom_vs_beta_pion->SetXTitle("Momentum (GeV)");
  mom_vs_beta_pion->SetYTitle("#beta");
  mom_vs_beta_pion->SetOption("COLZ");
  mom_vs_beta_pion->Write();

  mom_vs_beta_electron->SetXTitle("Momentum (GeV)");
  mom_vs_beta_electron->SetYTitle("#beta");
  mom_vs_beta_electron->SetOption("COLZ");
  mom_vs_beta_electron->Write();
}

void Histogram::Write_deltat() {
  deltat_proton->SetXTitle("Momentum (GeV)");
  deltat_proton->SetYTitle("#delta t");
  deltat_proton->SetOption("COLZ");
  deltat_proton->Write();
  deltat_pion->SetXTitle("Momentum (GeV)");
  deltat_pion->SetYTitle("#delta t");
  deltat_pion->SetOption("COLZ");
  deltat_pion->Write();
  deltat_electron->SetXTitle("Momentum (GeV)");
  deltat_electron->SetYTitle("#delta t");
  deltat_electron->SetOption("COLZ");
  deltat_electron->Write();
  deltat_proton_withID->SetXTitle("Momentum (GeV)");
  deltat_proton_withID->SetYTitle("#delta t");
  deltat_proton_withID->SetOption("COLZ");
  deltat_proton_withID->Write();
  deltat_proton_antiID->SetXTitle("Momentum (GeV)");
  deltat_proton_antiID->SetYTitle("#delta t");
  deltat_proton_antiID->SetOption("COLZ");
  deltat_proton_antiID->Write();
  deltat_pion_withID->SetXTitle("Momentum (GeV)");
  deltat_pion_withID->SetYTitle("#delta t");
  deltat_pion_withID->SetOption("COLZ");
  deltat_pion_withID->Write();
  deltat_electron_withID->SetXTitle("Momentum (GeV)");
  deltat_electron_withID->SetYTitle("#delta t");
  deltat_electron_withID->SetOption("COLZ");
  deltat_electron_withID->Write();
  deltat_electron_0th->SetXTitle("Momentum (GeV)");
  deltat_electron_0th->SetYTitle("#delta t");
  deltat_electron_0th->SetOption("COLZ");
  deltat_electron_0th->Write();
  deltat_electron_0th_ID->SetXTitle("Momentum (GeV)");
  deltat_electron_0th_ID->SetYTitle("#delta t");
  deltat_electron_0th_ID->SetOption("COLZ");
  deltat_electron_0th_ID->Write();
}
