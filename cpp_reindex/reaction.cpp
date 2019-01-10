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

  _hasE = false;
  _hasP = false;
  _hasPip = false;
  _hasPim = false;

  _MM = std::nan("-99");
  _MM2 = std::nan("-99");

  _W = std::nan("-99");
  _Q2 = std::nan("-99");
}

Reaction::~Reaction() {}

void Reaction::SetElec(float px, float py, float pz, float mass) {
  _hasE = true;
  _elec->SetXYZM(px, py, pz, mass);

  // Can calculate W and Q2 here
  _W = physics::W_calc(*_beam, *_elec);
  _Q2 = physics::Q2_calc(*_beam, *_elec);
}

void Reaction::SetProton(float px, float py, float pz, float mass) {
  _hasP = true;
  _prot->SetXYZM(px, py, pz, mass);
}
void Reaction::SetPip(float px, float py, float pz, float mass) {
  _hasPip = true;
  _pip->SetXYZM(px, py, pz, mass);
}
void Reaction::SetPim(float px, float py, float pz, float mass) {
  _hasPim = true;
  _pim->SetXYZM(px, py, pz, mass);
}

void Reaction::CalcMissMass() {
  TLorentzVector mm;
  mm = (*_beam - *_elec);
  mm += *_target;
  if (NeutronPip()) {
    mm -= *_pip;
    _MM = mm.M();
    _MM2 = mm.M2();
  } else if (twoPionEvent()) {
    mm -= *_prot;
    mm -= *_pip;
    mm -= *_pim;
    _MM = mm.M();
    _MM2 = mm.M2();
  } else if (ProtonPimEvent()) {
    mm -= *_prot;
    mm -= *_pim;
    _MM = mm.M();
    _MM2 = mm.M2();
  }
}

float Reaction::MM() {
  if (_MM != _MM) CalcMissMass();
  return _MM;
}
float Reaction::MM2() {
  if (_MM2 != _MM2) CalcMissMass();
  return _MM2;
}

float Reaction::W() { return _W; }
float Reaction::Q2() { return _Q2; }

bool Reaction::twoPionEvent() { return (_hasE && _hasP && _hasPip && _hasPim); }
bool Reaction::ProtonPimEvent() { return (_hasE && _hasP && _hasPim && !_hasPip); }
bool Reaction::NeutronPip() { return (_hasE && !_hasP && _hasPip && !_hasPim); }
