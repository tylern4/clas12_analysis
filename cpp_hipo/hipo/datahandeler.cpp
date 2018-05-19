/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#include "datahandeler.hpp"

DataHandeler::DataHandeler(std::vector<std::string> fin, std::string fout) {
  input_files = fin;
  if (getenv("CLAS12_E") != NULL) energy = atof(getenv("CLAS12_E"));
  e_mu = new TLorentzVector(0.0, 0.0, energy, energy);
  out = new TFile(fout.c_str(), "RECREATE");
  hist = new Histogram();
}

DataHandeler::~DataHandeler() {
  out->cd();
  TDirectory *wvsq2 = out->mkdir("wvsq2");
  wvsq2->cd();
  hist->Write_WvsQ2();

  TDirectory *mom_vs_beta = out->mkdir("mom_vs_beta");
  mom_vs_beta->cd();
  hist->Write_MomVsBeta();

  TDirectory *deltat_ftof = out->mkdir("deltat_ftof");
  deltat_ftof->cd();
  hist->Write_deltat();

  out->Close();
  std::cerr << "\n" << total << "\t" << std::endl;
}

void DataHandeler::run() {
#pragma omp parallel for
  for (std::vector<std::string>::const_iterator file_name = input_files.begin(); file_name != input_files.end();
       ++file_name) {
    std::cout << "Processing: " << *file_name << std::endl;
    file_handeler(*file_name);
  }
}

void DataHandeler::file_handeler(std::string fin) {
  // Load chain from branch h10
  hipo::reader reader;
  reader.open(fin.c_str());

  hipo::node<float> *beta = reader.getBranch<float>("REC::Particle", "beta");
  hipo::node<int8_t> *charge = reader.getBranch<int8_t>("REC::Particle", "charge");
  hipo::node<float> *chi2pid = reader.getBranch<float>("REC::Particle", "chi2pid");
  hipo::node<int32_t> *pid = reader.getBranch<int32_t>("REC::Particle", "pid");
  hipo::node<float> *px = reader.getBranch<float>("REC::Particle", "px");
  hipo::node<float> *py = reader.getBranch<float>("REC::Particle", "py");
  hipo::node<float> *pz = reader.getBranch<float>("REC::Particle", "pz");
  hipo::node<int16_t> *status = reader.getBranch<int16_t>("REC::Particle", "status");
  hipo::node<float> *vx = reader.getBranch<float>("REC::Particle", "vx");
  hipo::node<float> *vy = reader.getBranch<float>("REC::Particle", "vy");
  hipo::node<float> *vz = reader.getBranch<float>("REC::Particle", "vz");
  hipo::node<float> *evnt_BCG = reader.getBranch<float>("REC::Event", "BCG");
  hipo::node<float> *evnt_EVNTime = reader.getBranch<float>("REC::Event", "EVNTime");
  hipo::node<int16_t> *evnt_EvCAT = reader.getBranch<int16_t>("REC::Event", "EvCAT");
  hipo::node<int8_t> *evnt_Helic = reader.getBranch<int8_t>("REC::Event", "Helic");
  hipo::node<double> *evnt_LT = reader.getBranch<double>("REC::Event", "LT");
  hipo::node<int32_t> *evnt_NEVENT = reader.getBranch<int32_t>("REC::Event", "NEVENT");
  hipo::node<int16_t> *evnt_NPGP = reader.getBranch<int16_t>("REC::Event", "NPGP");
  hipo::node<int32_t> *evnt_NRUN = reader.getBranch<int32_t>("REC::Event", "NRUN");
  hipo::node<float> *evnt_PTIME = reader.getBranch<float>("REC::Event", "PTIME");
  hipo::node<float> *evnt_RFTime = reader.getBranch<float>("REC::Event", "RFTime");
  hipo::node<float> *evnt_STTime = reader.getBranch<float>("REC::Event", "STTime");
  hipo::node<int64_t> *evnt_TRG = reader.getBranch<int64_t>("REC::Event", "TRG");
  hipo::node<int8_t> *evnt_TYPE = reader.getBranch<int8_t>("REC::Event", "TYPE");
  hipo::node<float> *ec_chi2 = reader.getBranch<float>("REC::Calorimeter", "chi2");
  hipo::node<int8_t> *ec_detector = reader.getBranch<int8_t>("REC::Calorimeter", "detector");
  hipo::node<float> *ec_du = reader.getBranch<float>("REC::Calorimeter", "du");
  hipo::node<float> *ec_dv = reader.getBranch<float>("REC::Calorimeter", "dv");
  hipo::node<float> *ec_dw = reader.getBranch<float>("REC::Calorimeter", "dw");
  hipo::node<float> *ec_energy = reader.getBranch<float>("REC::Calorimeter", "energy");
  hipo::node<float> *ec_hx = reader.getBranch<float>("REC::Calorimeter", "hx");
  hipo::node<float> *ec_hy = reader.getBranch<float>("REC::Calorimeter", "hy");
  hipo::node<float> *ec_hz = reader.getBranch<float>("REC::Calorimeter", "hz");
  hipo::node<int16_t> *ec_index = reader.getBranch<int16_t>("REC::Calorimeter", "index");
  hipo::node<int8_t> *ec_layer = reader.getBranch<int8_t>("REC::Calorimeter", "layer");
  hipo::node<float> *ec_lu = reader.getBranch<float>("REC::Calorimeter", "lu");
  hipo::node<float> *ec_lv = reader.getBranch<float>("REC::Calorimeter", "lv");
  hipo::node<float> *ec_lw = reader.getBranch<float>("REC::Calorimeter", "lw");
  hipo::node<float> *ec_m2u = reader.getBranch<float>("REC::Calorimeter", "m2u");
  hipo::node<float> *ec_m2v = reader.getBranch<float>("REC::Calorimeter", "m2v");
  hipo::node<float> *ec_m2w = reader.getBranch<float>("REC::Calorimeter", "m2w");
  hipo::node<float> *ec_path = reader.getBranch<float>("REC::Calorimeter", "path");
  hipo::node<int16_t> *ec_pindex = reader.getBranch<int16_t>("REC::Calorimeter", "pindex");
  hipo::node<int8_t> *ec_sector = reader.getBranch<int8_t>("REC::Calorimeter", "sector");
  hipo::node<int16_t> *ec_status = reader.getBranch<int16_t>("REC::Calorimeter", "status");
  hipo::node<float> *ec_time = reader.getBranch<float>("REC::Calorimeter", "time");
  hipo::node<float> *ec_x = reader.getBranch<float>("REC::Calorimeter", "x");
  hipo::node<float> *ec_y = reader.getBranch<float>("REC::Calorimeter", "y");
  hipo::node<float> *ec_z = reader.getBranch<float>("REC::Calorimeter", "z");
  hipo::node<float> *sc_chi2 = reader.getBranch<float>("REC::Scintillator", "chi2");
  hipo::node<int16_t> *sc_component = reader.getBranch<int16_t>("REC::Scintillator", "component");
  hipo::node<int8_t> *sc_detector = reader.getBranch<int8_t>("REC::Scintillator", "detector");
  hipo::node<float> *sc_energy = reader.getBranch<float>("REC::Scintillator", "energy");
  hipo::node<float> *sc_hx = reader.getBranch<float>("REC::Scintillator", "hx");
  hipo::node<float> *sc_hy = reader.getBranch<float>("REC::Scintillator", "hy");
  hipo::node<float> *sc_hz = reader.getBranch<float>("REC::Scintillator", "hz");
  hipo::node<int16_t> *sc_index = reader.getBranch<int16_t>("REC::Scintillator", "index");
  hipo::node<int8_t> *sc_layer = reader.getBranch<int8_t>("REC::Scintillator", "layer");
  hipo::node<float> *sc_path = reader.getBranch<float>("REC::Scintillator", "path");
  hipo::node<int16_t> *sc_pindex = reader.getBranch<int16_t>("REC::Scintillator", "pindex");
  hipo::node<int8_t> *sc_sector = reader.getBranch<int8_t>("REC::Scintillator", "sector");
  hipo::node<int16_t> *sc_status = reader.getBranch<int16_t>("REC::Scintillator", "status");
  hipo::node<float> *sc_time = reader.getBranch<float>("REC::Scintillator", "time");
  hipo::node<float> *sc_x = reader.getBranch<float>("REC::Scintillator", "x");
  hipo::node<float> *sc_y = reader.getBranch<float>("REC::Scintillator", "y");
  hipo::node<float> *sc_z = reader.getBranch<float>("REC::Scintillator", "z");
  hipo::node<float> *cc_chi2 = reader.getBranch<float>("REC::Cherenkov", "chi2");
  hipo::node<int8_t> *cc_detector = reader.getBranch<int8_t>("REC::Cherenkov", "detector");
  hipo::node<float> *cc_dphi = reader.getBranch<float>("REC::Cherenkov", "dphi");
  hipo::node<float> *cc_dtheta = reader.getBranch<float>("REC::Cherenkov", "dtheta");
  hipo::node<int16_t> *cc_index = reader.getBranch<int16_t>("REC::Cherenkov", "index");
  hipo::node<float> *cc_nphe = reader.getBranch<float>("REC::Cherenkov", "nphe");
  hipo::node<float> *cc_path = reader.getBranch<float>("REC::Cherenkov", "path");
  hipo::node<float> *cc_phi = reader.getBranch<float>("REC::Cherenkov", "phi");
  hipo::node<int16_t> *cc_pindex = reader.getBranch<int16_t>("REC::Cherenkov", "pindex");
  hipo::node<int8_t> *cc_sector = reader.getBranch<int8_t>("REC::Cherenkov", "sector");
  hipo::node<int16_t> *cc_status = reader.getBranch<int16_t>("REC::Cherenkov", "status");
  hipo::node<float> *cc_theta = reader.getBranch<float>("REC::Cherenkov", "theta");
  hipo::node<float> *cc_time = reader.getBranch<float>("REC::Cherenkov", "time");
  hipo::node<float> *cc_x = reader.getBranch<float>("REC::Cherenkov", "x");
  hipo::node<float> *cc_y = reader.getBranch<float>("REC::Cherenkov", "y");
  hipo::node<float> *cc_z = reader.getBranch<float>("REC::Cherenkov", "z");

  double P;
  bool electron_cuts;
  int total = 0;
  int sc_d = 0;
  double W = 0;
  double Q2 = 0;
  double sf = 0;
  double P_x = 0;
  double P_y = 0;
  double P_z = 0;
  double per = 0;
  int index = 0;
  int num_pip = 0;
  TVector3 e_mu_prime_3;
  TLorentzVector e_mu_prime;
  bool good_e = false;

  int current_event = 0;
  while (reader.next() == true) {
    current_event++;
    if (pid->getLength() == 0) continue;

    if (!std::floor(current_event % 1000)) std::cerr << "\t\t" << std::floor(current_event) << "\r\r" << std::flush;
    num_pip = 0;
    good_e = false;
    for (int j = 0; j < ec_pindex->getLength(); j++) {
      if (ec_pindex->getLength() == 0) continue;
      try {
        index = ec_pindex->getValue(j);
        if (pid->getValue(index) == ELECTRON) {
          e_mu_prime_3.SetXYZ(px->getValue(index), py->getValue(index), pz->getValue(index));
          P = e_mu_prime_3.Mag();
          e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
          sf = ec_energy->getValue(j) / e_mu_prime.P();
          good_e = true;
        }
      } catch (std::exception &e) {
        total++;
      }
    }
    if (!good_e) continue;
    good_e = false;
    for (int j = 0; j < sc_time->getLength(); j++) {
      if (sc_time->getLength() == 0) continue;
      try {
        Delta_T *dt = new Delta_T(sc_time->getValue(0), sc_path->getValue(0));
        index = sc_pindex->getValue(j);
        sc_d = sc_detector->getValue(j);
        // I think 12 is FTOF
        if (sc_d == 12) {
          P_x = px->getValue(index) * px->getValue(index);
          P_y = py->getValue(index) * py->getValue(index);
          P_z = pz->getValue(index) * pz->getValue(index);
          P = TMath::Sqrt(P_x + P_y + P_z);

          dt->deltat(P, sc_time->getValue(j), sc_path->getValue(j));

          if (index == 0) {
            hist->Fill_MomVsBeta_vertex(pid->getValue(index), charge->getValue(index), P, beta->getValue(index));
            hist->Fill_deltat_vertex(pid->getValue(index), charge->getValue(index), P, dt);
          } else {
            hist->Fill_MomVsBeta(pid->getValue(index), charge->getValue(index), P, beta->getValue(index));
            hist->Fill_deltat(pid->getValue(index), charge->getValue(index), P, dt);
          }
        }
        if (pid->getValue(sc_pindex->getValue(j)) == PIP && abs(dt->Get_dt_Pi()) < 0.5) num_pip++;
        if (pid->getValue(sc_pindex->getValue(j)) == ELECTRON && sc_detector->getValue(sc_pindex->getValue(j)) == 12)
          good_e = true;
        delete dt;
      } catch (std::exception &e) {
        total++;
      }
    }
    if (!good_e) continue;
    // && sf >= 0.07 && sf <= 0.26
    if (good_e && e_mu_prime.P() > 1.5) {
      hist->Fill_EC(sf, e_mu_prime.P());
      W = physics::W_calc(*e_mu, e_mu_prime);
      Q2 = physics::Q2_calc(*e_mu, e_mu_prime);
      hist->Fill_WvsQ2(W, Q2);
      if (num_pip == 1 && pid->getLength() == 2) hist->Fill_WvsQ2_singlePi(W, Q2);
    }
  }
}
