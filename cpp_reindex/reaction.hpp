/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef REACTION_H_GUARD
#define REACTION_H_GUARD

#include "TLorentzVector.h"
#include "constants.hpp"
#include "physics.hpp"

class Reaction {
 private:
  TLorentzVector *_beam;
  TLorentzVector *_elec;
  TLorentzVector *_target;
  TLorentzVector *_prot;
  TLorentzVector *_pip;
  TLorentzVector *_pim;

  bool _hasE;
  bool _hasP;
  bool _hasPip;
  bool _hasPim;

  float _MM;
  float _MM2;

  float _W;
  float _Q2;

 public:
  Reaction();
  Reaction(TLorentzVector *beam);
  ~Reaction();

  void SetElec(float px, float py, float pz, float mass);
  void SetProton(float px, float py, float pz, float mass);
  void SetPip(float px, float py, float pz, float mass);
  void SetPim(float px, float py, float pz, float mass);
  void CalcMissMass();

  float MM();
  float MM2();
  float W();
  float Q2();

  bool twoPionEvent();
  bool ProtonPimEvent();
};

#endif
