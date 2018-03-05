/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "histogram.hpp"
#include "constants.hpp"

Histogram::Histogram() { makeHists_deltat(); }

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

void Histogram::makeHists_deltat() {
  for (size_t p = 0; p < particle_num; p++) {
    for (size_t c = 0; c < charge_num; c++) {
      for (size_t i = 0; i < with_id_num; i++) {
        hname.append("delta_t_");
        htitle.append("#Deltat ");
        hname.append(particle_name[p]);
        htitle.append(particle_name[p]);
        hname.append("_");
        htitle.append(" ");
        hname.append(charge_name[c]);
        htitle.append(charge_name[c]);
        hname.append("_");
        htitle.append(" ");
        hname.append(id_name[i]);
        htitle.append(id_name[i]);
        delta_t_hist[p][c][i] =
            new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, p_max, bins, Dt_min, Dt_max);
        hname.clear();
        htitle.clear();
      }
    }
  }
}

void Histogram::Fill_deltat(int pid, int charge, double P, Delta_T *dt) {
  double deltaT = -99;
  int good_ID = 0;

  for (size_t p = 0; p < particle_num; p++) {
    switch (p) {
      case 0:
        good_ID = ELECTRON;
        deltaT = dt->Get_dt_E();
        break;
      case 1:
        good_ID = PIP;
        deltaT = dt->Get_dt_Pi();
        break;
      case 2:
        good_ID = PROTON;
        deltaT = dt->Get_dt_P();
        break;
      case 3:
        good_ID = KP;
        deltaT = dt->Get_dt_K();
        break;
    }

    delta_t_hist[p][0][0]->Fill(P, deltaT);
    if (good_ID == abs(pid)) {
      delta_t_hist[p][0][1]->Fill(P, deltaT);
    } else {
      delta_t_hist[p][0][2]->Fill(P, deltaT);
    }

    if (charge == -1) {
      delta_t_hist[p][2][0]->Fill(P, deltaT);
      if (-good_ID == pid) {
        delta_t_hist[p][2][1]->Fill(P, deltaT);
      } else {
        delta_t_hist[p][2][2]->Fill(P, deltaT);
      }
    } else if (charge == 1) {
      delta_t_hist[p][1][0]->Fill(P, deltaT);
      if (good_ID == pid) {
        delta_t_hist[p][1][1]->Fill(P, deltaT);
      } else {
        delta_t_hist[p][1][2]->Fill(P, deltaT);
      }
    }
  }
}

void Histogram::Write_deltat() {
  for (size_t p = 0; p < particle_num; p++) {
    for (size_t c = 0; c < charge_num; c++) {
      for (size_t i = 0; i < with_id_num; i++) {
        delta_t_hist[p][c][i]->SetXTitle("Momentum (GeV)");
        delta_t_hist[p][c][i]->SetYTitle("#Deltat");
        delta_t_hist[p][c][i]->SetOption("COLZ");
        delta_t_hist[p][c][i]->Write();
      }
    }
  }
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
