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
    } else if (atof(getenv("CLAS12_E")) < 8) {
      q2_max = 4.0;
      w_max = 4.0;
      p_max = 4.0;
    } else if (atof(getenv("CLAS12_E")) < 9) {
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
  makeHists_sector();
  makeHists_deltat();
  makeHists_MomVsBeta();
}

Histogram::~Histogram() {}

// W and Q^2
void Histogram::Fill_WvsQ2(double W, double Q2, int sec) {
  W_vs_q2->Fill(W, Q2);
  W_hist->Fill(W);
  Q2_hist->Fill(Q2);

  if (sec > 0 && sec <= 6) {
    W_vs_q2_sec[sec - 1]->Fill(W, Q2);
    W_sec[sec - 1]->Fill(W);
  }
}

void Histogram::Fill_WvsQ2_det(double W, double Q2, int det) {
  if (det == 1) W_det[0]->Fill(W);
  if (det == 2) W_det[1]->Fill(W);
  if (det == 4) W_det[2]->Fill(W);
}

// W and Q^2
void Histogram::Fill_WvsQ2_singlePi(double W, double Q2, double mm, int sec) {
  W_vs_MM_singlePi[sec - 1]->Fill(W, mm);
  W_vs_q2_singlePi->Fill(W, Q2);
  W_hist_singlePi->Fill(W);
  Q2_hist_singlePi->Fill(Q2);
  MM_neutron->Fill(mm);
  if (sec > 0 && sec <= 6) {
    W_vs_q2_singlePi_sec[sec - 1]->Fill(W, Q2);
    W_singlePi_sec[sec - 1]->Fill(W);
    MM_neutron_sec[sec - 1]->Fill(mm);
  }
}

// W and Q^2
void Histogram::Fill_WvsQ2_Npip(double W, double Q2, double mm, int sec) {
  if (sec > 0 && sec <= 6) {
    W_vs_q2_Npip_sec[sec - 1]->Fill(W, Q2);
    W_Npip_sec[sec - 1]->Fill(W);
    MM_Npip_sec[sec - 1]->Fill(mm);
  }
}

void Histogram::Write_WvsQ2(TFile *out) {
  for (size_t i = 0; i < 3; i++) {
    W_det[i]->SetXTitle("W (GeV)");
    W_det[i]->Write();
  }
  auto WvsQ2_can = std::make_unique<TCanvas>("WvsQ2_can", "W vs Q2 sectors", 1920, 1080);
  WvsQ2_can->Divide(3, 2);
  for (size_t i = 0; i < num_sectors; i++) {
    W_vs_q2_sec[i]->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_q2_sec[i]->SetXTitle("W (GeV)");
    W_vs_q2_sec[i]->SetOption("COLZ");
    WvsQ2_can->cd(i + 1);
    W_vs_q2_sec[i]->Draw("same");
  }
  WvsQ2_can->Write();

  auto W_can = std::make_unique<TCanvas>("W_can", "W sectors", 1920, 1080);
  W_can->Divide(3, 2);
  for (size_t i = 0; i < num_sectors; i++) {
    W_sec[i]->SetXTitle("W (GeV)");
    W_can->cd(i + 1);
    W_sec[i]->Fit("gaus", "QMR+", "QMR+", 0.5, 1.2);
    W_sec[i]->Draw("same");
  }
  W_can->Write();

  W_vs_q2->SetXTitle("W (GeV)");
  W_vs_q2->SetYTitle("Q^{2} (GeV^{2})");
  W_vs_q2->SetOption("COLZ");
  W_vs_q2->Write();

  W_hist->SetXTitle("W (GeV)");
  W_hist->Write();

  Q2_hist->SetXTitle("Q^{2} (GeV^{2})");
  Q2_hist->Write();

  W_vs_q2_singlePi->SetXTitle("W (GeV)");
  W_vs_q2_singlePi->SetYTitle("Q^{2} (GeV^{2})");
  W_vs_q2_singlePi->SetOption("COLZ");
  W_vs_q2_singlePi->Write();

  W_hist_singlePi->SetXTitle("W (GeV)");
  W_hist_singlePi->Write();

  Q2_hist_singlePi->SetXTitle("Q^{2} (GeV^{2})");
  Q2_hist_singlePi->Write();
  MM_neutron->Write();

  auto wvsq2_sec = out->mkdir("wvsq2_sec");
  wvsq2_sec->cd();
  for (size_t i = 0; i < num_sectors; i++) {
    W_vs_q2_sec[i]->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_q2_sec[i]->SetXTitle("W (GeV)");
    W_vs_q2_sec[i]->SetOption("COLZ");
    W_vs_q2_sec[i]->Write();
  }
  auto w_sec = out->mkdir("w_sec");
  w_sec->cd();
  for (size_t i = 0; i < num_sectors; i++) {
    W_sec[i]->SetXTitle("W (GeV)");
    W_sec[i]->Write();
  }
  auto singlePi_sec = out->mkdir("singlePi_sec");
  singlePi_sec->cd();
  for (size_t i = 0; i < num_sectors; i++) {
    W_vs_q2_singlePi_sec[i]->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_q2_singlePi_sec[i]->SetXTitle("W (GeV)");
    W_vs_q2_singlePi_sec[i]->SetOption("COLZ");
    W_vs_q2_singlePi_sec[i]->Write();
  }
  for (size_t i = 0; i < num_sectors; i++) {
    W_vs_MM_singlePi[i]->SetOption("COLZ");
    W_vs_MM_singlePi[i]->SetYTitle("MM (GeV)");
    W_vs_MM_singlePi[i]->SetXTitle("W (GeV)");
    W_vs_MM_singlePi[i]->Write();
  }
  for (size_t i = 0; i < num_sectors; i++) {
    W_singlePi_sec[i]->SetXTitle("W (GeV)");
    W_singlePi_sec[i]->Write();
  }
  for (size_t i = 0; i < num_sectors; i++) {
    MM_neutron_sec[i]->Fit("gaus", "", "", 0.7, 1.1);
    MM_neutron_sec[i]->SetXTitle("Mass (GeV)");
    MM_neutron_sec[i]->Write();
  }

  auto Npip_sec = out->mkdir("Npip_sec");
  Npip_sec->cd();

  for (size_t i = 0; i < num_sectors; i++) {
    W_vs_q2_Npip_sec[i]->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_q2_Npip_sec[i]->SetXTitle("W (GeV)");
    W_vs_q2_Npip_sec[i]->SetOption("COLZ");
    W_vs_q2_Npip_sec[i]->Write();
  }
  for (size_t i = 0; i < num_sectors; i++) {
    W_Npip_sec[i]->SetXTitle("W (GeV)");
    W_Npip_sec[i]->Write();
  }
  for (size_t i = 0; i < num_sectors; i++) {
    MM_Npip_sec[i]->SetXTitle("Mass (GeV)");
    MM_Npip_sec[i]->Write();
  }
}

void Histogram::makeHists_sector() {
  for (size_t i = 0; i < 3; i++) {
    hname.append("W_det_");
    htitle.append("W detector: ");
    hname.append(std::to_string(i + 1));
    htitle.append(std::to_string(i + 1));
    W_det[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, zero, w_max);
    hname.clear();
    htitle.clear();
  }

  for (size_t i = 0; i < num_sectors; i++) {
    hname.clear();
    htitle.clear();
    hname.append("wvsq2_sec_");
    htitle.append("W vs Q^{2} Sector: ");
    hname.append(std::to_string(i + 1));
    htitle.append(std::to_string(i + 1));
    W_vs_q2_sec[i] = new TH2D(hname.c_str(), htitle.c_str(), bins, zero, w_max, bins, zero, q2_max);
    hname.clear();
    htitle.clear();

    hname.append("w_sec_");
    htitle.append("W Sector: ");
    hname.append(std::to_string(i + 1));
    htitle.append(std::to_string(i + 1));
    W_sec[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, zero, w_max);
    hname.clear();
    htitle.clear();

    hname.clear();
    htitle.clear();
    hname.append("wvsq2_sec_singlePi_");
    htitle.append("W vs Q^{2} W_singlePi Sector: ");
    hname.append(std::to_string(i + 1));
    htitle.append(std::to_string(i + 1));
    W_vs_q2_singlePi_sec[i] = new TH2D(hname.c_str(), htitle.c_str(), bins, zero, w_max, bins, zero, q2_max);
    hname.clear();
    htitle.clear();

    hname.append("w_sec_singlePi_");
    htitle.append("W singlePi Sector: ");
    hname.append(std::to_string(i + 1));
    htitle.append(std::to_string(i + 1));
    W_singlePi_sec[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, zero, w_max);
    hname.clear();
    htitle.clear();

    hname.clear();
    htitle.clear();
    hname.append("wvsq2_sec_Npip_");
    htitle.append("W vs Q^{2} W_Npip Sector: ");
    hname.append(std::to_string(i + 1));
    htitle.append(std::to_string(i + 1));
    W_vs_q2_Npip_sec[i] = new TH2D(hname.c_str(), htitle.c_str(), bins, zero, w_max, bins, zero, q2_max);
    hname.clear();
    htitle.clear();

    hname.append("w_sec_Npip_");
    htitle.append("W Npip Sector: ");
    hname.append(std::to_string(i + 1));
    htitle.append(std::to_string(i + 1));
    W_Npip_sec[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, zero, w_max);
    hname.clear();
    htitle.clear();

    hname.append("MM_Sec_");
    htitle.append("MM neutron Sector: ");
    hname.append(std::to_string(i + 1));
    htitle.append(std::to_string(i + 1));
    MM_neutron_sec[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, zero, 4.0);
    hname.clear();
    htitle.clear();

    hname.append("MM_Npip_Sec_");
    htitle.append("MM^{2} neutron pip Sector: ");
    hname.append(std::to_string(i + 1));
    htitle.append(std::to_string(i + 1));
    MM_Npip_sec[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, -2.0, 2.0);
    hname.clear();
    htitle.clear();

    hname.append("W_vs_MM_singlePi_");
    htitle.append("W_vs_MM_singlePi_");
    hname.append(std::to_string(i + 1));
    htitle.append(std::to_string(i + 1));
    W_vs_MM_singlePi[i] = new TH2D(hname.c_str(), htitle.c_str(), bins, zero, w_max, bins, -q2_max, q2_max);
    hname.clear();
    htitle.clear();
  }
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

void Histogram::Fill_deltat_pi(int pid, int charge, float dt, float momentum) {
  if (charge == 1) {
    delta_t_hist[1][0][0]->Fill(momentum, dt);
    if (pid == PIP)
      delta_t_hist[1][0][1]->Fill(momentum, dt);
    else
      delta_t_hist[1][0][2]->Fill(momentum, dt);
  } else if (charge == -1) {
    delta_t_hist[1][1][0]->Fill(momentum, dt);
    if (pid == -PIP)
      delta_t_hist[1][1][1]->Fill(momentum, dt);
    else
      delta_t_hist[1][1][2]->Fill(momentum, dt);
  }
}

void Histogram::Fill_deltat_prot(int pid, int charge, float dt, float momentum) {
  if (charge == 1) {
    delta_t_hist[2][0][0]->Fill(momentum, dt);
    if (pid == PROTON)
      delta_t_hist[2][0][1]->Fill(momentum, dt);
    else
      delta_t_hist[2][0][2]->Fill(momentum, dt);
  } else if (charge == -1) {
    delta_t_hist[2][1][0]->Fill(momentum, dt);
    if (pid == -PROTON)
      delta_t_hist[2][1][1]->Fill(momentum, dt);
    else
      delta_t_hist[2][1][2]->Fill(momentum, dt);
  }
}

void Histogram::Write_deltat() {
  for (size_t i = 0; i < with_id_num; i++) {
    delta_t_vertex[i]->SetXTitle("Momentum (GeV)");
    delta_t_vertex[i]->SetYTitle("#Deltat");
    delta_t_vertex[i]->SetOption("COLZ");
    if (delta_t_vertex[i]->GetEntries() > 1) delta_t_vertex[i]->Write();
    delete delta_t_vertex[i];
  }

  for (size_t p = 0; p < particle_num; p++) {
    for (size_t c = 0; c < charge_num; c++) {
      for (size_t i = 0; i < with_id_num; i++) {
        delta_t_hist[p][c][i]->SetXTitle("Momentum (GeV)");
        delta_t_hist[p][c][i]->SetYTitle("#Deltat");
        delta_t_hist[p][c][i]->SetOption("COLZ");
        if (delta_t_hist[p][c][i]->GetEntries() > 1) delta_t_hist[p][c][i]->Write();
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
