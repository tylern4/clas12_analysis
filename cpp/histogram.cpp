/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "histogram.hpp"

Histogram::Histogram(const std::string& output_file) {
  RootOutputFile = std::make_shared<TFile>(output_file.c_str(), "RECREATE");
  def = new TCanvas("def");
  if (getenv("BEAM_E") != NULL) {
    if (atof(getenv("BEAM_E")) < 3) {
      q2_max = 1.0;
      w_max = 3.5;
      p_max = 3.0;
    } else if (atof(getenv("BEAM_E")) < 8) {
      q2_max = 4.0;
      w_max = 4.0;
      p_max = 4.0;
    } else if (atof(getenv("BEAM_E")) < 9) {
      q2_max = 7.0;
      w_max = 7.0;
      p_max = 7.0;
    }
  }
  // Kinematics
  momentum = std::make_shared<TH1D>("mom", "mom", bins, p_min, p_max);
  W_hist = std::make_shared<TH1D>("W", "W", bins, zero, w_max);
  Q2_hist = std::make_shared<TH1D>("Q2", "Q2", bins, zero, q2_max);
  W_vs_q2 = std::make_shared<TH2D>("W_vs_q2", "W_vs_q2", bins, zero, w_max, bins, zero, q2_max);

  MM_neutron = std::make_shared<TH1D>("missMass", "missMass", bins, zero, 4.0);

  W_hist_singlePi = std::make_shared<TH1D>("W_singlePi", "W_singlePi", bins, zero, w_max);
  Q2_hist_singlePi = std::make_shared<TH1D>("Q2_singlePi", "Q2_singlePi", bins, zero, q2_max);
  W_vs_q2_singlePi =
      std::make_shared<TH2D>("W_vs_q2_singlePi", "W_vs_q2_singlePi", bins, zero, w_max, bins, zero, q2_max);

  EC_sampling_fraction =
      std::make_shared<TH2D>("EC_sampling_fraction", "EC_sampling_fraction", bins, p_min, p_max, bins, zero, 1.0);
  makeHists_sector();
  makeHists_deltat();
  makeHists_MomVsBeta();
}

Histogram::~Histogram() { this->Write(); }

void Histogram::Write() {
  std::cout << GREEN << "Writting" << DEF << std::endl;

  std::cerr << BOLDBLUE << "WvsQ2()" << DEF << std::endl;
  TDirectory* WvsQ2_folder = RootOutputFile->mkdir("W vs Q2");
  WvsQ2_folder->cd();
  Write_WvsQ2();

  std::cerr << BOLDBLUE << "Write_MomVsBeta()" << DEF << std::endl;
  TDirectory* Write_MomVsBeta_folder = RootOutputFile->mkdir("Mom Vs Beta");
  Write_MomVsBeta_folder->cd();
  Write_MomVsBeta();

  std::cerr << BOLDBLUE << "Write_deltat()" << DEF << std::endl;
  // TDirectory* Write_deltat_folder = RootOutputFile->mkdir("Delta_t");
  // Write_deltat_folder->cd();
  Write_deltat();

  std::cerr << BOLDBLUE << "Done Writing!!!" << DEF << std::endl;
}

// W and Q^2
void Histogram::Fill_WvsQ2(double W, double Q2, int sec, float weight) {
  W_vs_q2->Fill(W, Q2);
  W_hist->Fill(W);
  Q2_hist->Fill(Q2);

  if (sec > 0 && sec <= 6) {
    W_vs_q2_sec[sec - 1]->Fill(W, Q2, weight);
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

void Histogram::Write_WvsQ2() {
  for (short i = 0; i < 3; i++) {
    W_det[i]->SetXTitle("W (GeV)");
    W_det[i]->Write();
  }
  auto WvsQ2_can = std::make_unique<TCanvas>("WvsQ2_can", "W vs Q2 sectors", 1920, 1080);
  WvsQ2_can->Divide(3, 2);
  for (short i = 0; i < num_sectors; i++) {
    W_vs_q2_sec[i]->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_q2_sec[i]->SetXTitle("W (GeV)");
    W_vs_q2_sec[i]->SetOption("COLZ");
    WvsQ2_can->cd(i + 1);
    W_vs_q2_sec[i]->Draw("same");
  }
  WvsQ2_can->Write();

  auto W_can = std::make_unique<TCanvas>("W_can", "W sectors", 1920, 1080);
  W_can->Divide(3, 2);
  for (short i = 0; i < num_sectors; i++) {
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

  auto wvsq2_sec = RootOutputFile->mkdir("wvsq2_sec");
  wvsq2_sec->cd();
  for (short i = 0; i < num_sectors; i++) {
    W_vs_q2_sec[i]->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_q2_sec[i]->SetXTitle("W (GeV)");
    W_vs_q2_sec[i]->SetOption("COLZ");
    W_vs_q2_sec[i]->Write();
  }
  auto w_sec = RootOutputFile->mkdir("w_sec");
  w_sec->cd();
  for (short i = 0; i < num_sectors; i++) {
    W_sec[i]->SetXTitle("W (GeV)");
    W_sec[i]->Write();
  }
  auto singlePi_sec = RootOutputFile->mkdir("singlePi_sec");
  singlePi_sec->cd();
  for (short i = 0; i < num_sectors; i++) {
    W_vs_q2_singlePi_sec[i]->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_q2_singlePi_sec[i]->SetXTitle("W (GeV)");
    W_vs_q2_singlePi_sec[i]->SetOption("COLZ");
    W_vs_q2_singlePi_sec[i]->Write();
  }
  for (short i = 0; i < num_sectors; i++) {
    W_vs_MM_singlePi[i]->SetOption("COLZ");
    W_vs_MM_singlePi[i]->SetYTitle("MM (GeV)");
    W_vs_MM_singlePi[i]->SetXTitle("W (GeV)");
    W_vs_MM_singlePi[i]->Write();
  }
  for (short i = 0; i < num_sectors; i++) {
    W_singlePi_sec[i]->SetXTitle("W (GeV)");
    W_singlePi_sec[i]->Write();
  }
  for (short i = 0; i < num_sectors; i++) {
    MM_neutron_sec[i]->Fit("gaus", "QMR+", "QMR+", 0.7, 1.1);
    MM_neutron_sec[i]->SetXTitle("Mass (GeV)");
    MM_neutron_sec[i]->Write();
  }

  auto Npip_sec = RootOutputFile->mkdir("Npip_sec");
  Npip_sec->cd();

  for (short i = 0; i < num_sectors; i++) {
    W_vs_q2_Npip_sec[i]->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_q2_Npip_sec[i]->SetXTitle("W (GeV)");
    W_vs_q2_Npip_sec[i]->SetOption("COLZ");
    W_vs_q2_Npip_sec[i]->Write();
  }
  for (short i = 0; i < num_sectors; i++) {
    W_Npip_sec[i]->SetXTitle("W (GeV)");
    W_Npip_sec[i]->Write();
  }
  for (short i = 0; i < num_sectors; i++) {
    MM_Npip_sec[i]->SetXTitle("Mass (GeV)");
    MM_Npip_sec[i]->Write();
  }
}

void Histogram::makeHists_sector() {
  for (short i = 0; i < 3; i++) {
    W_det[i] = std::make_shared<TH1D>(Form("W_det_%d", i + 1), Form("W detector: %d", i + 1), bins, zero, w_max);
  }

  for (short i = 0; i < num_sectors; i++) {
    W_vs_q2_sec[i] = std::make_shared<TH2D>(Form("wvsq2_sec_%d", i + 1), Form("W vs Q^{2} Sector: %d", i + 1), bins,
                                            zero, w_max, bins, zero, q2_max);

    W_sec[i] = std::make_shared<TH1D>(Form("w_sec_%d", i + 1), Form("W Sector: %d", i + 1), bins, zero, w_max);

    W_vs_q2_singlePi_sec[i] =
        std::make_shared<TH2D>(Form("wvsq2_sec_singlePi_%d", i + 1), Form("W vs Q^{2} W_singlePi Sector: %d", i + 1),
                               bins, zero, w_max, bins, zero, q2_max);

    W_singlePi_sec[i] = std::make_shared<TH1D>(Form("w_sec_singlePi_%d", i + 1), Form("W singlePi Sector: %d", i + 1),
                                               bins, zero, w_max);

    W_vs_q2_Npip_sec[i] =
        std::make_shared<TH2D>(Form("wvsq2_sec_Npip_%d", i + 1), Form("W vs Q^{2} W_Npip Sector: %d", i + 1), bins,
                               zero, w_max, bins, zero, q2_max);

    W_Npip_sec[i] =
        std::make_shared<TH1D>(Form("w_sec_Npip_%d", i + 1), Form("W Npip Sector: %d", i + 1), bins, zero, w_max);

    MM_neutron_sec[i] =
        std::make_shared<TH1D>(Form("MM_Sec_%d", i + 1), Form("MM neutron Sector: %d", i + 1), bins, zero, 4.0);

    MM_Npip_sec[i] = std::make_shared<TH1D>(Form("MM_Npip_Sec_%d", i + 1), Form("MM^{2} neutron pip Sector: %d", i + 1),
                                            bins, -2.0, 2.0);

    W_vs_MM_singlePi[i] = std::make_shared<TH2D>(Form("W_vs_MM_singlePi_%d", i + 1), Form("W_vs_MM_singlePi_%d", i + 1),
                                                 bins, zero, w_max, bins, -q2_max, q2_max);
  }
}

void Histogram::makeHists_deltat() {
  std::string tof = "";
  for (short p = 0; p < particle_num; p++) {
    for (short c = 0; c < charge_num; c++) {
      for (short i = 0; i < with_id_num; i++) {
        tof = "ftof";
        delta_t_hist[p][c][i][0] =
            std::make_shared<TH2D>(Form("delta_t_%s_%s_%s_%s", tof.c_str(), particle_name[p].c_str(),
                                        charge_name[c].c_str(), id_name[i].c_str()),
                                   Form("#Deltat %s %s %s %s", tof.c_str(), particle_name[p].c_str(),
                                        charge_name[c].c_str(), id_name[i].c_str()),
                                   bins, p_min, p_max, bins, Dt_min, Dt_max);

        tof = "ctof";
        delta_t_hist[p][c][i][1] =
            std::make_shared<TH2D>(Form("delta_t_%s_%s_%s_%s", tof.c_str(), particle_name[p].c_str(),
                                        charge_name[c].c_str(), id_name[i].c_str()),
                                   Form("#Deltat %s %s %s %s", tof.c_str(), particle_name[p].c_str(),
                                        charge_name[c].c_str(), id_name[i].c_str()),
                                   bins, 0, 3.0, bins, -6.0, 6.0);
      }
    }
  }
}

void Histogram::Fill_deltat_pi(int pid, int charge, float dt, float momentum, bool fc) {
  if (charge == 1) {
    delta_t_hist[1][0][0][fc]->Fill(momentum, dt);
    if (pid == PIP)
      delta_t_hist[1][0][1][fc]->Fill(momentum, dt);
    else
      delta_t_hist[1][0][2][fc]->Fill(momentum, dt);
  } else if (charge == -1) {
    delta_t_hist[1][1][0][fc]->Fill(momentum, dt);
    if (pid == -PIP)
      delta_t_hist[1][1][1][fc]->Fill(momentum, dt);
    else
      delta_t_hist[1][1][2][fc]->Fill(momentum, dt);
  }
}

void Histogram::Fill_deltat_prot(int pid, int charge, float dt, float momentum, bool fc) {
  if (charge == 1) {
    delta_t_hist[2][0][0][fc]->Fill(momentum, dt);
    if (pid == PROTON)
      delta_t_hist[2][0][1][fc]->Fill(momentum, dt);
    else
      delta_t_hist[2][0][2][fc]->Fill(momentum, dt);
  } else if (charge == -1) {
    delta_t_hist[2][1][0][fc]->Fill(momentum, dt);
    if (pid == -PROTON)
      delta_t_hist[2][1][1][fc]->Fill(momentum, dt);
    else
      delta_t_hist[2][1][2][fc]->Fill(momentum, dt);
  }
}

void Histogram::Write_deltat() {
  TDirectory* ftof_folder = RootOutputFile->mkdir("ftof");
  ftof_folder->cd();
  for (short p = 0; p < particle_num; p++) {
    for (short c = 0; c < charge_num; c++) {
      for (short i = 0; i < with_id_num; i++) {
        delta_t_hist[p][c][i][0]->SetXTitle("Momentum (GeV)");
        delta_t_hist[p][c][i][0]->SetYTitle("#Deltat");
        delta_t_hist[p][c][i][0]->SetOption("COLZ");
        if (delta_t_hist[p][c][i][0]->GetEntries() > 1) delta_t_hist[p][c][i][0]->Write();
      }
    }
  }
  TDirectory* ctof_folder = RootOutputFile->mkdir("ctof");
  ctof_folder->cd();
  for (short p = 0; p < particle_num; p++) {
    for (short c = 0; c < charge_num; c++) {
      for (short i = 0; i < with_id_num; i++) {
        delta_t_hist[p][c][i][1]->SetXTitle("Momentum (GeV)");
        delta_t_hist[p][c][i][1]->SetYTitle("#Deltat");
        delta_t_hist[p][c][i][1]->SetOption("COLZ");
        if (delta_t_hist[p][c][i][1]->GetEntries() > 1) delta_t_hist[p][c][i][1]->Write();
      }
    }
  }
}

void Histogram::makeHists_MomVsBeta() {
  for (short p = 0; p < particle_num; p++) {
    for (short c = 0; c < charge_num; c++) {
      for (short i = 0; i < with_id_num; i++) {
        momvsbeta_hist[p][c][i] = std::make_shared<TH2D>(
            Form("mom_vs_beta_%s_%s_%s", particle_name[p].c_str(), charge_name[c].c_str(), id_name[i].c_str()),
            Form("Momentum vs #beta %s %s %s", particle_name[p].c_str(), charge_name[c].c_str(), id_name[i].c_str()),
            bins, p_min, p_max, bins, zero, 1.2);
      }
    }
  }
}

void Histogram::Fill_MomVsBeta(int pid, int charge, double P, double beta) {
  int good_ID = 0;
  if (beta != 0) {
    momentum->Fill(P);
    for (short p = 0; p < particle_num; p++) {
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
  momentum->SetXTitle("Momentum (GeV)");
  momentum->Write();
  for (short p = 0; p < particle_num; p++) {
    for (short c = 0; c < charge_num; c++) {
      for (short i = 0; i < with_id_num; i++) {
        momvsbeta_hist[p][c][i]->SetXTitle("Momentum (GeV)");
        momvsbeta_hist[p][c][i]->SetYTitle("#beta");
        momvsbeta_hist[p][c][i]->SetOption("COLZ");
        momvsbeta_hist[p][c][i]->Write();
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
