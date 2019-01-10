/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "histogram.hpp"

Histogram::Histogram() {
  if (getenv("CLAS12_E") != NULL) {
    if (atof(getenv("CLAS12_E")) < 3) {
      q2_max = 1.0;
      w_max = 3.5;
      p_max = 3.0;
    } else if (atof(getenv("CLAS12_E")) < 7) {
      q2_max = 7.0;
      w_max = 7.0;
      p_max = 7.0;
    }
  }
  // Kinematics
  momentum = new TH1D("mom", "mom", bins, p_min, p_max);
  W_hist = new TH1D("W", "W", bins, zero, w_max);
  Q2_hist = new TH1D("Q2", "Q2", bins, zero, q2_max);
  W_vs_q2 = new TH2D("W_vs_q2", "W_vs_q2", bins, zero, w_max, bins, zero, q2_max);

  MM_neutron = new TH1D("missMass", "missMass", bins, zero, 4.0);

  W_hist_singlePi = new TH1D("W_singlePi", "W_singlePi", bins, zero, w_max);
  Q2_hist_singlePi = new TH1D("Q2_singlePi", "Q2_singlePi", bins, zero, q2_max);
  W_vs_q2_singlePi = new TH2D("W_vs_q2_singlePi", "W_vs_q2_singlePi", bins, zero, w_max, bins, zero, q2_max);

  EC_sampling_fraction = new TH2D("EC_sampling_fraction", "EC_sampling_fraction", bins, p_min, p_max, bins, zero, 1.0);

  makeHists_deltat();
  makeHists_MomVsBeta();
}

Histogram::~Histogram() {}

// W and Q^2
void Histogram::Fill_WvsQ2(double W, double Q2) {
  W_vs_q2->Fill(W, Q2);
  W_hist->Fill(W);
  Q2_hist->Fill(Q2);
}

// W and Q^2
void Histogram::Fill_WvsQ2_singlePi(double W, double Q2, double mm) {
  W_vs_q2_singlePi->Fill(W, Q2);
  W_hist_singlePi->Fill(W);
  Q2_hist_singlePi->Fill(Q2);
  MM_neutron->Fill(mm);
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
  /*
    W_vs_q2_lower->SetXTitle("W (GeV)");
    W_vs_q2_lower->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_q2_lower->SetOption("COLZ");
    W_vs_q2_lower->Write();

    W_hist_lower->SetXTitle("W (GeV)");
    W_hist_lower->Write();

    Q2_hist_lower->SetXTitle("Q^{2} (GeV^{2})");
    Q2_hist_lower->Write();

    W_vs_q2_upper->SetXTitle("W (GeV)");
    W_vs_q2_upper->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_q2_upper->SetOption("COLZ");
    W_vs_q2_upper->Write();

    W_hist_upper->SetXTitle("W (GeV)");
    W_hist_upper->Write();

    Q2_hist_upper->SetXTitle("Q^{2} (GeV^{2})");
    Q2_hist_upper->Write();
  */
  W_vs_q2_singlePi->SetXTitle("W (GeV)");
  W_vs_q2_singlePi->SetYTitle("Q^{2} (GeV^{2})");
  W_vs_q2_singlePi->SetOption("COLZ");
  W_vs_q2_singlePi->Write();

  W_hist_singlePi->SetXTitle("W (GeV)");
  W_hist_singlePi->Write();

  Q2_hist_singlePi->SetXTitle("Q^{2} (GeV^{2})");
  Q2_hist_singlePi->Write();
  MM_neutron->Write();

  /*
    W_vs_q2_lower_singlePi->SetXTitle("W (GeV)");
    W_vs_q2_lower_singlePi->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_q2_lower_singlePi->SetOption("COLZ");
    W_vs_q2_lower_singlePi->Write();

    W_hist_lower_singlePi->SetXTitle("W (GeV)");
    W_hist_lower_singlePi->Write();

    Q2_hist_lower_singlePi->SetXTitle("Q^{2} (GeV^{2})");
    Q2_hist_lower_singlePi->Write();

    W_vs_q2_upper_singlePi->SetXTitle("W (GeV)");
    W_vs_q2_upper_singlePi->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_q2_upper_singlePi->SetOption("COLZ");
    W_vs_q2_upper_singlePi->Write();

    W_hist_upper_singlePi->SetXTitle("W (GeV)");
    W_hist_upper_singlePi->Write();

    Q2_hist_upper_singlePi->SetXTitle("Q^{2} (GeV^{2})");
    Q2_hist_upper_singlePi->Write();
  */
}

void Histogram::makeHists_deltat() {
  for (size_t i = 0; i < with_id_num; i++) {
    hname.append("delta_t_vertex");
    htitle.append("#Deltat vertex particle");
    hname.append("_");
    htitle.append(" ");
    hname.append(id_name[i]);
    htitle.append(id_name[i]);
    delta_t_vertex[i] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, p_max, bins, Dt_min, Dt_max);
    hname.clear();
    htitle.clear();
  }

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
        delta_t_hist[p][c][i] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, p_max, bins, Dt_min, Dt_max);
        hname.clear();
        htitle.clear();
      }
    }
  }
}

void Histogram::Fill_deltat_vertex(int pid, int charge, float dt, float momentum) {
  delta_t_vertex[0]->Fill(momentum, dt);
  if (pid == ELECTRON) {
    delta_t_vertex[1]->Fill(momentum, dt);
  } else {
    delta_t_vertex[2]->Fill(momentum, dt);
  }
}

void Histogram::Fill_deltat_pip(int pid, int charge, float dt, float momentum) {
  delta_t_hist[1][1][0]->Fill(momentum, dt);
}

void Histogram::Write_deltat() {
  for (size_t i = 0; i < with_id_num; i++) {
    delta_t_vertex[i]->SetXTitle("Momentum (GeV)");
    delta_t_vertex[i]->SetYTitle("#Deltat");
    delta_t_vertex[i]->SetOption("COLZ");
    delta_t_vertex[i]->Write();
    delete delta_t_vertex[i];
  }

  for (size_t p = 0; p < particle_num; p++) {
    for (size_t c = 0; c < charge_num; c++) {
      for (size_t i = 0; i < with_id_num; i++) {
        delta_t_hist[p][c][i]->SetXTitle("Momentum (GeV)");
        delta_t_hist[p][c][i]->SetYTitle("#Deltat");
        delta_t_hist[p][c][i]->SetOption("COLZ");
        delta_t_hist[p][c][i]->Write();
        delete delta_t_hist[p][c][i];
      }
    }
  }
}

void Histogram::makeHists_MomVsBeta() {
  for (size_t i = 0; i < with_id_num; i++) {
    hname.append("mom_vs_beta_vertex");
    htitle.append("Momentum vs #beta vertex");
    hname.append("_");
    htitle.append(" ");
    hname.append(id_name[i]);
    htitle.append(id_name[i]);
    momvsbeta_vertex[i] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, p_max, bins, zero, 1.2);
    hname.clear();
    htitle.clear();
  }

  for (size_t p = 0; p < particle_num; p++) {
    for (size_t c = 0; c < charge_num; c++) {
      for (size_t i = 0; i < with_id_num; i++) {
        hname.append("mom_vs_beta_");
        htitle.append("Momentum vs #beta ");
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
        momvsbeta_hist[p][c][i] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, p_max, bins, zero, 1.2);
        hname.clear();
        htitle.clear();
      }
    }
  }
}

void Histogram::Fill_MomVsBeta_vertex(int pid, int charge, double P, double beta) {
  if (beta != 0) {
    momvsbeta_vertex[0]->Fill(P, beta);
    if (pid == ELECTRON) {
      momvsbeta_vertex[1]->Fill(P, beta);
    } else {
      momvsbeta_vertex[2]->Fill(P, beta);
    }
  }
}

void Histogram::Fill_MomVsBeta(int pid, int charge, double P, double beta) {
  int good_ID = 0;
  if (beta != 0) {
    momentum->Fill(P);
    for (size_t p = 0; p < particle_num; p++) {
      switch (p) {
        case 0:
          good_ID = ELECTRON;
          break;
        case 1:
          good_ID = PIP;
          break;
        case 2:
          good_ID = PROTON;
          break;
        case 3:
          good_ID = KP;
          break;
      }

      momvsbeta_hist[p][0][0]->Fill(P, beta);
      if (good_ID == abs(pid)) {
        momvsbeta_hist[p][0][1]->Fill(P, beta);
      } else {
        momvsbeta_hist[p][0][2]->Fill(P, beta);
      }

      if (charge == -1) {
        momvsbeta_hist[p][2][0]->Fill(P, beta);
        if (-good_ID == pid) {
          momvsbeta_hist[p][2][1]->Fill(P, beta);
        } else {
          momvsbeta_hist[p][2][2]->Fill(P, beta);
        }
      } else if (charge == 1) {
        momvsbeta_hist[p][1][0]->Fill(P, beta);
        if (good_ID == pid) {
          momvsbeta_hist[p][1][1]->Fill(P, beta);
        } else {
          momvsbeta_hist[p][1][2]->Fill(P, beta);
        }
      }
    }
  }
}

void Histogram::Write_MomVsBeta() {
  for (size_t i = 0; i < with_id_num; i++) {
    momvsbeta_vertex[i]->SetXTitle("Momentum (GeV)");
    momvsbeta_vertex[i]->SetYTitle("#beta");
    momvsbeta_vertex[i]->SetOption("COLZ");
    momvsbeta_vertex[i]->Write();
    delete momvsbeta_vertex[i];
  }

  momentum->SetXTitle("Momentum (GeV)");
  momentum->Write();
  for (size_t p = 0; p < particle_num; p++) {
    for (size_t c = 0; c < charge_num; c++) {
      for (size_t i = 0; i < with_id_num; i++) {
        momvsbeta_hist[p][c][i]->SetXTitle("Momentum (GeV)");
        momvsbeta_hist[p][c][i]->SetYTitle("#beta");
        momvsbeta_hist[p][c][i]->SetOption("COLZ");
        momvsbeta_hist[p][c][i]->Write();
        delete momvsbeta_hist[p][c][i];
      }
    }
  }
}

void Histogram::Fill_EC(double sf, double momentum) { EC_sampling_fraction->Fill(momentum, sf); }
void Histogram::Write_EC() {
  EC_sampling_fraction->SetXTitle("Momentum (GeV)");
  EC_sampling_fraction->SetYTitle("Sampling Fraction");
  EC_sampling_fraction->SetOption("COLZ");
  EC_sampling_fraction->Write();
}
