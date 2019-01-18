/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef DT_H_GUARD
#define DT_H_GUARD
#include <iostream>
#include <map>
#include "constants.hpp"

class Delta_T {
 private:
  std::map<int, double> _mass_map = {{PROTON, MASS_P}, {-PROTON, MASS_P},  {NEUTRON, MASS_N},  {PIP, MASS_PIP},
                                     {PIM, MASS_PIM},  {PI0, MASS_PI0},    {KP, MASS_KP},      {KM, MASS_KM},
                                     {PHOTON, MASS_G}, {ELECTRON, MASS_E}, {-ELECTRON, MASS_E}};
  float _sc_t_v = std::nanf("-99");
  float _sc_r_v = std::nanf("-99");
  float _vertex = std::nanf("-99");
  float _sc_t = std::nanf("-99");
  float _sc_r = std::nanf("-99");
  float _beta = std::nanf("-99");
  float _momentum = std::nanf("-99");

  bool _ctof = false;

  float _vertex_time(float sc_time, float sc_pathlength, float relatavistic_beta);
  float _deltat(int num);

 public:
  Delta_T(float time_1b, float path_1b, float time_1a, float path_1a, float time_2, float path_2);
  ~Delta_T();

  void dt_calc(float momentum, float time_1b, float path_1b, float time_1a, float path_1a, float time_2, float path_2,
               float time_ctof, float path_ctof);
  float dt_E();
  float dt_P();
  float dt_Pi();
  float dt_K();
  float dt(int pid);
  float momentum();
};

#endif
