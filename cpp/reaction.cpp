/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "reaction.hpp"

Reaction::Reaction() {
  double energy = CLAS12_E;
  if (getenv("CLAS12_E") != NULL) energy = atof(getenv("CLAS12_E"));
  _beam = std::make_unique<TLorentzVector>();
  _beam->SetPxPyPzE(0.0, 0.0, sqrt(energy * energy - MASS_E * MASS_E), energy);

  _target = std::make_unique<TLorentzVector>(0.0, 0.0, 0.0, MASS_P);
  _elec = std::make_unique<TLorentzVector>();
  _prot = std::make_unique<TLorentzVector>();
  _pip = std::make_unique<TLorentzVector>();
  _pim = std::make_unique<TLorentzVector>();
  _neutron = std::make_unique<TLorentzVector>();
  _other = std::make_unique<TLorentzVector>();
}

Reaction::~Reaction() {}

void Reaction::SetElec(float px, float py, float pz) {
  _hasE = true;
  _elec->SetXYZM(px, py, pz, MASS_E);

  _W = physics::W_calc(*_beam, *_elec);
  _Q2 = physics::Q2_calc(*_beam, *_elec);
}

void Reaction::SetProton(float px, float py, float pz) {
  _numProt++;
  _numPos++;
  _hasP = true;
  _prot->SetXYZM(px, py, pz, MASS_P);
}
void Reaction::SetPip(float px, float py, float pz) {
  _numPip++;
  _numPos++;
  _hasPip = true;
  _pip->SetXYZM(px, py, pz, MASS_PIP);
}
void Reaction::SetPim(float px, float py, float pz) {
  _numPim++;
  _numNeg++;
  _hasPim = true;
  _pim->SetXYZM(px, py, pz, MASS_PIM);
}

void Reaction::SetOther(float px, float py, float pz, int pid) {
  if (pid == NEUTRON) {
    _hasNeutron = true;
    _numNeutral++;
    _neutron->SetXYZM(px, py, pz, _mass_map[pid]);
  } else {
    _hasOther = true;
    _numOther++;
    _other->SetXYZM(px, py, pz, _mass_map[pid]);
  }
}

void Reaction::CalcMissMass() {
  auto mm = std::make_unique<TLorentzVector>();
  *mm += (*_beam - *_elec + *_target);
  if (SinglePip() || NeutronPip()) {
    *mm -= *_pip;
    _MM = mm->M();
    _MM2 = mm->M2();
  } else if (TwoPion()) {
    *mm -= *_prot;
    *mm -= *_pip;
    *mm -= *_pim;
    _MM = mm->M();
    _MM2 = mm->M2();
  } else if (ProtonPim()) {
    *mm -= *_prot;
    *mm -= *_pim;
    _MM = mm->M();
    _MM2 = mm->M2();
  } else if (SingleP()) {
    *mm -= *_prot;
    _MM = mm->M();
    _MM2 = mm->M2();
  }
}

float Reaction::MM() {
  CalcMissMass();
  return _MM;
}
float Reaction::MM2() {
  CalcMissMass();
  return _MM2;
}

float Reaction::W() { return _W; }
float Reaction::Q2() { return _Q2; }
