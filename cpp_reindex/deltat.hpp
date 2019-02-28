/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef DT_H_GUARD
#define DT_H_GUARD
#include <iostream>
#include <map>
#include "clas12branches.hpp"
#include "constants.hpp"

class Delta_T {
 private:
  std::map<int, double> _mass_map = {{PROTON, MASS_P}, {-PROTON, MASS_P},  {NEUTRON, MASS_N},  {PIP, MASS_PIP},
                                     {PIM, MASS_PIM},  {PI0, MASS_PI0},    {KP, MASS_KP},      {KM, MASS_KM},
                                     {PHOTON, MASS_G}, {ELECTRON, MASS_E}, {-ELECTRON, MASS_E}};
  float _sc_t_v = NAN;
  float _sc_r_v = NAN;
  float _vertex = NAN;
  float _sc_t = NAN;
  float _sc_r = NAN;
  float _beta = NAN;
  float _momentum = NAN;

  float _time_1b_vert = NAN;
  float _path_1b_vert = NAN;
  float _time_1a_vert = NAN;
  float _path_1a_vert = NAN;
  float _time_2_vert = NAN;
  float _path_2_vert = NAN;

  std::shared_ptr<Clas12Branches> _data;

  float _vertex_time(float sc_time, float sc_pathlength, float relatavistic_beta);
  float _deltat(int num);

 public:
  Delta_T(std::shared_ptr<Clas12Branches> data);
  ~Delta_T();

  void get_det_info(int i);
  float dt_E(int i);
  float dt_P(int i);
  float dt_Pi(int i);
  float dt_K(int i);
  float dt(int i, int pid);
  float momentum(int i);
};

#endif
