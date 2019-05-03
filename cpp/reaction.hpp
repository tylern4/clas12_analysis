/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef REACTION_H_GUARD
#define REACTION_H_GUARD

#include <iostream>
#include "TLorentzVector.h"
#include "constants.hpp"
#include "physics.hpp"

class Reaction {
 private:
  std::unique_ptr<TLorentzVector> _beam;
  std::unique_ptr<TLorentzVector> _elec;
  std::unique_ptr<TLorentzVector> _gamma;
  std::unique_ptr<TLorentzVector> _target;
  std::unique_ptr<TLorentzVector> _prot;
  std::unique_ptr<TLorentzVector> _pip;
  std::unique_ptr<TLorentzVector> _pim;
  std::unique_ptr<TLorentzVector> _neutron;
  std::unique_ptr<TLorentzVector> _other;

  bool _hasE = false;
  bool _hasP = false;
  bool _hasPip = false;
  bool _hasPim = false;
  bool _hasOther = false;
  bool _hasNeutron = false;

  bool _boosted = false;

  short _numProt = 0;
  short _numPip = 0;
  short _numPim = 0;
  short _numPos = 0;
  short _numNeg = 0;
  short _numNeutral = 0;
  short _numOther = 0;

  float _MM = std::nan("-99");
  float _MM2 = std::nan("-99");

  float _W = std::nan("-99");
  float _Q2 = std::nan("-99");

  std::map<int, double> _mass_map = {{PROTON, MASS_P}, {-PROTON, MASS_P},  {NEUTRON, MASS_N},  {PIP, MASS_PIP},
                                     {PIM, MASS_PIM},  {PI0, MASS_PI0},    {KP, MASS_KP},      {KM, MASS_KM},
                                     {PHOTON, MASS_G}, {ELECTRON, MASS_E}, {-ELECTRON, MASS_E}};

 public:
  Reaction();
  ~Reaction();

  void SetElec(float px, float py, float pz);
  void SetProton(float px, float py, float pz);
  void SetPip(float px, float py, float pz);
  void SetPim(float px, float py, float pz);
  void SetOther(float px, float py, float pz, int pid);
  void CalcMissMass();

  float MM();
  float MM2();
  float W();
  float Q2();

  bool TwoPion() { return ((_numPip == 1 && _numPim == 1) && (_hasE && !_hasP && _hasPip && _hasPim)); }
  bool ProtonPim() { return ((_numProt == 1 && _numPim == 1) && (_hasE && _hasP && !_hasPip && _hasPim)); }
  bool SinglePip() { return ((_numPip == 1) && (_hasE && !_hasP && _hasPip && !_hasPim)); }
  bool SingleP() { return ((_numProt == 1) && (_hasE && _hasP && !_hasPip && !_hasPim)); }
  bool NeutronPip() { return ((_numPip == 1 && _numNeutral == 1) && (_hasE && !_hasP && _hasPip && !_hasPim)); }
};

#endif
