/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "deltat.hpp"

Delta_T::Delta_T(std::shared_ptr<Branches12> data) {
  _data = data;
  if (!std::isnan(_data->sc_ftof_1b_time(0))) {
    _sc_t_v = _data->sc_ftof_1b_time(0);
    _sc_r_v = _data->sc_ftof_1b_path(0);
  } else if (!std::isnan(_data->sc_ftof_1a_time(0))) {
    _sc_t_v = _data->sc_ftof_1a_time(0);
    _sc_r_v = _data->sc_ftof_1a_path(0);
  } else if (!std::isnan(_data->sc_ftof_2_time(0))) {
    _sc_t_v = _data->sc_ftof_2_time(0);
    _sc_r_v = _data->sc_ftof_2_path(0);
  }
  _vertex = _vertex_time(_sc_t_v, _sc_r_v, 1.0);
}

Delta_T::~Delta_T() {}

float Delta_T::_vertex_time(float sc_time, float sc_pathlength, float relatavistic_beta) {
  return sc_time - sc_pathlength / (relatavistic_beta * c_special_units);
}

float Delta_T::_deltat(int pid) {
  _beta = 1.0 / sqrt(1.0 + (mass[pid] / _momentum) * (mass[pid] / _momentum));
  if (_sc_t == _sc_t && _sc_r == _sc_r) {
    return _vertex - _vertex_time(_sc_t, _sc_r, _beta);
  } else {
    return std::nanf("-99");
  }
}

void Delta_T::dt_calc(float momentum, float time_1b, float path_1b, float time_1a, float path_1a, float time_2,
                      float path_2, float time_ctof, float path_ctof) {
  _momentum = momentum;
  if (time_1b == time_1b) {
    _sc_t = time_1b;
    _sc_r = path_1b;
  } else if (time_1a == time_1a) {
    _sc_t = time_1a;
    _sc_r = path_1a;
  } else if (time_ctof == time_ctof) {
    _ctof = true;
    _sc_t = time_ctof;
    _sc_r = path_ctof;
  } else if (time_2 == time_2) {
    _sc_t = time_2;
    _sc_r = path_2;
  }
}

void Delta_T::dt_calc(int i) {
  _momentum = _data->p(i);
  if (!std::isnan(_data->sc_ftof_1b_time(i))) {
    _sc_t = _data->sc_ftof_1b_time(i);
    _sc_r = _data->sc_ftof_1b_path(i);
  } else if (!std::isnan(_data->sc_ftof_1a_time(i))) {
    _sc_t = _data->sc_ftof_1a_time(i);
    _sc_r = _data->sc_ftof_1a_path(i);
  } else if (!std::isnan(_data->sc_ctof_time(i))) {
    _ctof = true;
    _sc_t = _data->sc_ctof_time(i);
    _sc_r = _data->sc_ctof_path(i);
  } else if (!std::isnan(_data->sc_ftof_2_time(i))) {
    _sc_t = _data->sc_ftof_2_time(i);
    _sc_r = _data->sc_ftof_2_path(i);
  }
}

float Delta_T::dt_E() { return _deltat(ELECTRON); }
float Delta_T::dt_P() { return _deltat(PROTON); }
float Delta_T::dt_Pi() { return _deltat(PIP); }
float Delta_T::dt_K() { return _deltat(KP); }
float Delta_T::dt(int pid) { return _deltat(pid); }
float Delta_T::momentum() { return _momentum; }
bool Delta_T::ctof() { return _ctof; }
