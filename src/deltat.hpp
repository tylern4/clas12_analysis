/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef DT_H_GUARD
#define DT_H_GUARD
#include <iostream>
#include "branches.hpp"

class Delta_T {
 private:
  std::shared_ptr<Branches12> _data;
  float _sc_t_v = NAN;
  float _sc_r_v = NAN;
  float _vertex = NAN;
  float _sc_t = NAN;
  float _sc_r = NAN;
  float _beta = NAN;
  float _momentum = NAN;
  bool _ctof = false;

  float _vertex_time(float sc_time, float sc_pathlength, float relatavistic_beta);
  float _deltat(int num);

 public:
  Delta_T(std::shared_ptr<Branches12> data);
  ~Delta_T();

  void dt_calc(float momentum, float time_1b, float path_1b, float time_1a, float path_1a, float time_2, float path_2,
               float time_ctof, float path_ctof);
  void dt_calc(int i);
  float dt_E();
  float dt_P();
  float dt_Pi();
  float dt_K();
  float dt(int pid);
  float momentum();
  bool ctof();
};

#endif
